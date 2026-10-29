"""
Finding and reading the inputs.

`ieeg` (the BIDS-iEEG derivative of ieeg-prep, or a plain BIDS-iEEG dataset) gives the
recordings, their good channels and the expert labels; `methods` (a fan-in port) gives the
tables of seizure-dynasd, seizure-ei and seizure-fragility. With links each port arrives in
INPUT_DIR/<port>/; without links everything sits side by side in INPUT_DIR, so files are found
by their names anywhere under the port (symlinks followed).
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("consensus")

# The method tables, by file name (conventions.md, "Tables written by the analysis steps").
METHOD_FILES = {
    "seizure_onsets.csv": "dynasd",
    "channel_onsets.csv": "dynasd",
    "ei.csv": "ei",
    "fragility.csv": "fragility",
}
METHODS = ("dynasd", "ei", "fragility")

TRUE = {"true", "1", "yes", "y", "t"}
FALSE = {"false", "0", "no", "n", "f"}


def walk(root: Path) -> list[Path]:
    """Every regular file under root (symlinked files and folders followed, hidden skipped)."""
    out: list[Path] = []
    if not root.is_dir():
        return out
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:  # a symlink loop
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        out += [Path(dirpath) / f for f in sorted(filenames) if not f.startswith(".")]
    return [p for p in out if p.is_file()]


def port_root(input_dir: Path, port: str) -> Path:
    """INPUT_DIR/<port>/ when the port was linked, else INPUT_DIR (everything side by side)."""
    sub = input_dir / port
    return sub if sub.is_dir() and walk(sub) else input_dir


def to_bool(value) -> bool | None:
    """'true'/'false' (and yes/no, 1/0) to a bool; empty or unknown to None."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    s = str(value).strip().lower()
    if s in TRUE:
        return True
    if s in FALSE:
        return False
    return None


def to_float(value) -> float | None:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(v) else v


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix == ".tsv" else ","
    # BIDS files may carry a BOM; "n/a" is BIDS for missing.
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, encoding="utf-8-sig",
                       na_values=[]).replace({"n/a": ""})


# ---- the ieeg port ----------------------------------------------------------------------

@dataclass
class Recording:
    recording_id: str
    channels: list[str] = field(default_factory=list)  # good channels, in file order
    onset_labels: dict[str, bool] = field(default_factory=dict)  # channel -> onset, when known
    onset_s: float | None = None  # expert onset (recordings.csv, else the events.tsv seizure row)
    seizure_end_s: float | None = None
    duration_s: float | None = None
    stim_induced: bool | None = None
    lvfa: bool | None = None
    channels_file: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def has_labels(self) -> bool:
        return any(self.onset_labels.values())


def stem_of(name: str, suffix: str) -> str:
    """Recording id from a BIDS file name: the entity chain without desc-preproc and the suffix."""
    s = name[: -len(suffix)] if name.endswith(suffix) else name
    return re.sub(r"_desc-preproc$", "", s)


@dataclass
class IeegIndex:
    recordings: dict[str, Recording]
    index_path: Path | None
    columns: list[str]  # recordings.csv columns (to tell whether stim_induced / lvfa exist)


def read_ieeg(root: Path) -> IeegIndex:
    files = walk(root)
    channels = {stem_of(p.name, "_channels.tsv"): p for p in files if p.name.endswith("_channels.tsv")}
    events = {stem_of(p.name, "_events.tsv"): p for p in files if p.name.endswith("_events.tsv")}
    index = sorted((p for p in files if p.name == "recordings.csv"), key=lambda p: (len(p.parts), str(p)))
    if len(index) > 1:
        log.info("several recordings.csv found, using %s", index[0])
    recs: dict[str, Recording] = {}
    columns: list[str] = []
    rows: list[dict] = []
    if index:
        df = read_table(index[0])
        columns = list(df.columns)
        rows = df.to_dict("records")
    else:
        rows = [{"recording_id": rid} for rid in sorted(channels)]

    for row in rows:
        rid = str(row.get("recording_id", "")).strip()
        if not rid:
            continue
        rec = Recording(rid)
        rec.onset_s = to_float(row.get("onset_s"))
        rec.duration_s = to_float(row.get("duration_s"))
        rec.stim_induced = to_bool(row.get("stim_induced")) if "stim_induced" in row else None
        rec.lvfa = to_bool(row.get("lvfa")) if "lvfa" in row else None
        ch = channels.get(rid)
        if ch is None and row.get("path"):
            # recordings.csv path -> <stem>_desc-preproc_ieeg.edf: look beside it.
            edf = Path(str(row["path"])).name
            ch = channels.get(stem_of(re.sub(r"_ieeg\.[A-Za-z0-9]+$", "", edf), ""))
        if ch is not None:
            read_channels(rec, ch, root)
        else:
            rec.notes.append("no _channels.tsv")
        ev = events.get(rid)
        if ev is not None:
            read_events(rec, ev)
        recs[rid] = rec
    return IeegIndex(recs, index[0] if index else None, columns)


def read_channels(rec: Recording, path: Path, root: Path) -> None:
    df = read_table(path)
    if "name" not in df.columns:
        rec.notes.append(f"{path.name} has no name column")
        return
    status = df["status"].str.strip().str.lower() if "status" in df.columns else pd.Series([""] * len(df))
    good = df[status != "bad"]
    rec.channels = [str(c).strip() for c in good["name"]]
    rec.channels_file = os.path.relpath(path, root)
    if "onset" in df.columns:
        labels = {str(n).strip(): to_bool(v) for n, v in zip(good["name"], good["onset"])}
        # A label column that is empty for every channel means "unknown", not "no onset".
        if any(v is not None for v in labels.values()):
            rec.onset_labels = {k: bool(v) for k, v in labels.items() if v is not None}


def read_events(rec: Recording, path: Path) -> None:
    df = read_table(path)
    if not {"onset", "trial_type"} <= set(df.columns):
        return
    sz = df[df["trial_type"].str.strip().str.lower() == "seizure"]
    if sz.empty:
        return
    on = to_float(sz.iloc[0]["onset"])
    dur = to_float(sz.iloc[0].get("duration")) if "duration" in sz.columns else None
    if rec.onset_s is None and on is not None:
        rec.onset_s = on
        rec.notes.append("onset_s from events.tsv")
    if on is not None and dur is not None and dur > 0:
        rec.seizure_end_s = on + dur


# ---- the methods port -------------------------------------------------------------------

@dataclass
class MethodTables:
    seizure_onsets: pd.DataFrame | None = None
    channel_onsets: pd.DataFrame | None = None
    ei: pd.DataFrame | None = None
    fragility: pd.DataFrame | None = None
    paths: dict[str, str] = field(default_factory=dict)
    ignored: list[str] = field(default_factory=list)

    @property
    def present(self) -> list[str]:
        out = []
        if self.channel_onsets is not None:
            out.append("dynasd")
        if self.ei is not None:
            out.append("ei")
        if self.fragility is not None:
            out.append("fragility")
        return out


REQUIRED = {
    "seizure_onsets.csv": {"recording_id", "onset_s"},
    "channel_onsets.csv": {"recording_id", "channel"},
    "ei.csv": {"recording_id", "channel", "ei"},
    "fragility.csv": {"recording_id", "channel", "fragility_post"},
}


def read_methods(root: Path, input_dir: Path) -> MethodTables:
    found: dict[str, list[Path]] = {}
    for p in walk(root):
        if p.name in METHOD_FILES:
            found.setdefault(p.name, []).append(p)
    out = MethodTables()
    for name, paths in sorted(found.items()):
        # Shallowest first; a second table of the same kind (e.g. a step run twice under
        # different LABEL folders) is ignored and logged.
        paths.sort(key=lambda p: (len(p.parts), str(p)))
        for extra in paths[1:]:
            out.ignored.append(os.path.relpath(extra, input_dir))
            log.info("ignoring %s: already using %s", extra, paths[0])
        df = read_table(paths[0])
        missing = REQUIRED[name] - set(df.columns)
        if missing:
            log.info("ignoring %s: lacks column(s) %s", paths[0], ", ".join(sorted(missing)))
            out.ignored.append(os.path.relpath(paths[0], input_dir))
            continue
        df["recording_id"] = df["recording_id"].str.strip()
        if "channel" in df.columns:
            df["channel"] = df["channel"].str.strip()
        setattr(out, name[:-4], df)
        out.paths[name] = os.path.relpath(paths[0], input_dir)
    return out
