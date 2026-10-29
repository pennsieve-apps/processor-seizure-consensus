"""
A tiny synthetic input for processor-seizure-consensus, following conventions.md:

- the `ieeg` port of ieeg-prep: recordings.csv, per recording _channels.tsv (with `onset`
  labels) and _events.tsv (`seizure` row). The consensus step never reads the signals, so no
  EDF is written (recordings.csv still names its path);
- the `methods` fan-in: seizure-dynasd's seizure_onsets.csv + channel_onsets.csv, seizure-ei's
  ei.csv and seizure-fragility's fragility.csv, with the onset channels scored higher (plus
  noise) so the expected answer is known.

Recordings: P1 seizure 1 (labelled, spontaneous), P1 seizure 2 (labelled, stimulation-induced,
LVFA; ei.csv lacks two channels), P2 seizure 1 (no labels: `onset` column empty), and a channel
marked bad in P2.

Layouts: `linked` (INPUT_DIR/ieeg/..., INPUT_DIR/methods/<step>/... with each step's
provenance.json, as the methods apps write them) and `flat` (everything
side by side, nested folders). `flags=False` leaves out the proposed stim_induced/lvfa columns.

    python -m tests.make_fixture data/fixture [--no-flags]
"""
from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

import numpy as np

CHANNELS = [f"{s}{i:02d}-{s}{i + 1:02d}" for s in ("LA", "LB") for i in range(1, 9)]  # 16 bipolar
RECORDINGS = [
    # recording id, onset_s, onset channel indices (None = unlabelled), stim, lvfa
    ("sub-P1_ses-implant_task-ictal1000_run-01", 120.0, [2, 3, 4], False, False),
    ("sub-P1_ses-implant_task-ictal2000_run-01", 118.5, [9, 10, 11], True, True),
    ("sub-P2_ses-implant_task-ictal3000_run-01", 121.0, None, False, False),
]
BAD = ("sub-P2_ses-implant_task-ictal3000_run-01", "LB08-LB09")


def write_tsv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_ieeg(root: Path, flags: bool) -> None:
    rec_rows = []
    for rid, onset, onset_ch, stim, lvfa in RECORDINGS:
        sub, ses = rid.split("_")[0], rid.split("_")[1]
        d = root / sub / ses / "ieeg"
        rows = []
        for i, ch in enumerate(CHANNELS):
            a, b = ch.split("-")
            bad = (rid, ch) == BAD
            label = "" if onset_ch is None else str(i in onset_ch).lower()
            rows.append([ch, "SEEG", "µV", "n/a", "n/a", 512, "bad" if bad else "good",
                         "flat" if bad else "n/a", f"{a}|{b}", "", "", label, ""])
        write_tsv(d / f"{rid}_desc-preproc_channels.tsv",
                  ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status",
                   "status_description", "contacts", "soz", "resected", "onset", "spread"], rows)
        write_tsv(d / f"{rid}_desc-preproc_events.tsv", ["onset", "duration", "trial_type", "source"],
                  [[onset, 45.0, "seizure", "expert"], [onset - 100, 60.0, "baseline", "prep"]])
        row = [rid, sub[4:], ses[4:], rid.split("_")[2][5:], "01", f"{sub}/{ses}/ieeg/{rid}_desc-preproc_ieeg.edf",
               512, 16, 15 if rid == BAD[0] else 16, 300.0, "bipolar", "ictal", onset, onset - 100, onset - 40, "bids-edf"]
        if flags:
            row += [str(stim).lower(), str(lvfa).lower()]
        rec_rows.append(row)
    header = ["recording_id", "subject", "session", "task", "run", "path", "sfreq", "n_channels", "n_good",
              "duration_s", "montage", "kind", "onset_s", "baseline_start_s", "baseline_end_s", "source_format"]
    if flags:
        header += ["stim_induced", "lvfa"]
    write_csv(root / "recordings.csv", header, rec_rows)
    (root / "dataset_description.json").write_text(
        '{"Name": "synthetic", "BIDSVersion": "1.9.0", "DatasetType": "derivative", "GeneratedBy": [{"Name": "ieeg-prep"}]}\n')


def method_tables(seed: int = 42) -> dict[str, tuple[list[str], list[list]]]:
    rng = np.random.default_rng(seed)
    so, co, ei, fr = [], [], [], []
    for rid, onset, onset_ch, _, _ in RECORDINGS:
        truth = onset_ch if onset_ch is not None else [5, 6]
        so.append([rid, round(onset + 1.5, 3), "dynasd"])
        for i, ch in enumerate(CHANNELS):
            hot = i in truth
            if hot or i % 3 == 0:  # onset channels and a few others are recruited
                t = onset + (0.5 + rng.uniform(0, 1) if hot else 8 + rng.uniform(0, 10))
                co.append([rid, ch, round(t, 3), "true", round(rng.uniform(0.6, 1.0) if hot else rng.uniform(0.2, 0.6), 4), "dynasd"])
            else:
                co.append([rid, ch, "", "false", round(rng.uniform(0, 0.3), 4), "dynasd"])
            if not (rid.endswith("2000_run-01") and i in (14, 15)):  # ei lacks two channels here
                e = (0.7 + 0.3 * rng.uniform()) if hot else 0.4 * rng.uniform()
                ei.append([rid, ch, round(e, 4), round(onset + rng.uniform(0, 5), 3), round(rng.uniform(1, 5), 3)])
            pre = 0.4 + 0.1 * rng.uniform()
            post = (0.75 + 0.2 * rng.uniform()) if hot else (0.3 + 0.35 * rng.uniform())
            fr.append([rid, ch, round(pre, 4), round(post, 4), 0])
        # fragility rank: 1 = highest fragility_post within the recording
        rows = [r for r in fr if r[0] == rid]
        for k, r in enumerate(sorted(rows, key=lambda r: -r[3])):
            r[4] = k + 1
    return {
        "seizure_onsets.csv": (["recording_id", "onset_s", "method"], so),
        "channel_onsets.csv": (["recording_id", "channel", "onset_s", "recruited", "score", "method"], co),
        "ei.csv": (["recording_id", "channel", "ei", "detection_s", "ratio_peak"], ei),
        "fragility.csv": (["recording_id", "channel", "fragility_pre", "fragility_post", "rank"], fr),
    }


STEP = {"seizure_onsets.csv": "dynasd", "channel_onsets.csv": "dynasd", "ei.csv": "ei", "fragility.csv": "fragility"}


def make_fixture(out: Path, layout: str = "both", flags: bool = True) -> dict[str, Path]:
    out = Path(out)
    roots = {}
    tables = method_tables()
    for lay in (("linked", "flat") if layout == "both" else (layout,)):
        root = out / lay
        if root.exists():
            shutil.rmtree(root)
        if lay == "linked":
            write_ieeg(root / "ieeg", flags)
            for name, (header, rows) in tables.items():
                write_csv(root / "methods" / STEP[name] / name, header, rows)
            # as the platform merges the fan-in: each method's folder, provenance.json included
            for step in sorted(set(STEP.values())):
                (root / "methods" / step / "provenance.json").write_text(
                    '{"app": "processor-seizure-%s"}\n' % step)
        else:
            write_ieeg(root / "prep", flags)
            for name, (header, rows) in tables.items():
                write_csv(root / "results" / STEP[name] / name, header, rows)
        roots[lay] = root
    return roots


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--layout", default="both", choices=["both", "linked", "flat"])
    ap.add_argument("--no-flags", action="store_true")
    a = ap.parse_args()
    for k, v in make_fixture(Path(a.out), a.layout, not a.no_flags).items():
        print(k, v)
