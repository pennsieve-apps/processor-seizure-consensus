"""
processor-seizure-consensus: fuses the channel rankings of seizure-dynasd, seizure-ei and
seizure-fragility into one seizure-onset-zone estimate per recording, and scores every method
and the consensus against the expert onset labels when the recording has them.

Reads INPUT_DIR/ieeg/ and INPUT_DIR/methods/ (or everything side by side in INPUT_DIR), writes
OUTPUT_DIR/{onset_zone.csv, scores.csv, provenance.json} (the single `onset-zone` output).
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from processor import inputs, rank, scoring

log = logging.getLogger("consensus")
APP = "processor-seizure-consensus"

PARAMETERS = ("FUSION", "FRAGILITY_SCORE", "MISSING_CHANNELS", "ZONE_MODE", "ZONE_K",
              "ZONE_THRESHOLD", "EVENT_TOLERANCE_START_S", "EVENT_TOLERANCE_END_S")

ZONE_COLUMNS = ["recording_id", "channel", "dynasd_onset_s", "dynasd_score", "dynasd_recruited", "ei",
                "fragility", "rank_dynasd", "rank_ei", "rank_fragility", "n_methods", "consensus_score",
                "rank", "in_onset_zone", "expert_onset"]


def version() -> str:
    for p in (Path(__file__).resolve().parent.parent / "VERSION", Path("/processor/VERSION")):
        if p.is_file():
            return p.read_text().strip()
    return "unknown"


def env(name: str, default: str) -> str:
    v = os.environ.get(name, "")
    return v.strip() if v.strip() else default


def params() -> dict:
    p = {
        "FUSION": env("FUSION", "mean-rank").lower(),
        "FRAGILITY_SCORE": env("FRAGILITY_SCORE", "post").lower(),
        "MISSING_CHANNELS": env("MISSING_CHANNELS", "worst").lower(),
        "ZONE_MODE": env("ZONE_MODE", "top-k").lower(),
        "ZONE_K": env("ZONE_K", "10"),
        "ZONE_THRESHOLD": env("ZONE_THRESHOLD", "0.9"),
        "EVENT_TOLERANCE_START_S": env("EVENT_TOLERANCE_START_S", "30"),
        "EVENT_TOLERANCE_END_S": env("EVENT_TOLERANCE_END_S", "60"),
    }
    choices = {"FUSION": rank.FUSIONS, "FRAGILITY_SCORE": ("post", "difference"),
               "MISSING_CHANNELS": rank.MISSING, "ZONE_MODE": rank.ZONE_MODES}
    for k, allowed in choices.items():
        if p[k] not in allowed:
            raise ValueError(f"{k} must be one of {', '.join(allowed)}, got '{p[k]}'")
    try:
        p["ZONE_K"] = int(p["ZONE_K"])
        p["ZONE_THRESHOLD"] = float(p["ZONE_THRESHOLD"])
        p["EVENT_TOLERANCE_START_S"] = float(p["EVENT_TOLERANCE_START_S"])
        p["EVENT_TOLERANCE_END_S"] = float(p["EVENT_TOLERANCE_END_S"])
    except ValueError as e:
        raise ValueError(f"bad numeric parameter: {e}") from None
    if p["ZONE_K"] < 1:
        raise ValueError("ZONE_K must be >= 1")
    if not 0 <= p["ZONE_THRESHOLD"] <= 1:
        raise ValueError("ZONE_THRESHOLD must be in [0, 1]")
    return p


def detected_onset(tables: inputs.MethodTables, rid: str) -> tuple[float | None, str]:
    """seizure-dynasd's seizure onset: the earliest onset_s in seizure_onsets.csv, else the
    earliest recruited channel in channel_onsets.csv."""
    if tables.seizure_onsets is not None:
        t = tables.seizure_onsets[tables.seizure_onsets["recording_id"] == rid]
        v = pd.to_numeric(t["onset_s"], errors="coerce").dropna()
        if len(v):
            return float(v.min()), "seizure_onsets.csv"
    if tables.channel_onsets is not None:
        vals = rank.method_values("dynasd", tables.channel_onsets, rid, "post")
        v = vals.loc[vals["dynasd_recruited"], "dynasd_onset_s"].dropna()
        if len(v):
            return float(v.min()), "channel_onsets.csv (earliest recruited channel)"
    return None, ""


def process_recording(rid: str, rec: inputs.Recording | None, tables: inputs.MethodTables, p: dict):
    """onset_zone rows and score rows of one recording."""
    method_tables = {"dynasd": tables.channel_onsets, "ei": tables.ei, "fragility": tables.fragility}
    values = {m: rank.method_values(m, t, rid, p["FRAGILITY_SCORE"])
              for m, t in method_tables.items() if t is not None}
    values = {m: v for m, v in values.items() if v["_usable"].any()}
    notes: list[str] = list(rec.notes) if rec else ["not in the ieeg input"]
    channels = list(rec.channels) if rec and rec.channels else []
    if not channels:
        seen: dict[str, None] = {}
        for v in values.values():
            seen.update(dict.fromkeys(v["channel"]))
        channels = list(seen)
        if channels:
            notes.append("channel set from the method tables")
    if not values:
        raise RuntimeError("no method has scores for this recording")
    if not channels:
        raise RuntimeError("no channels")
    for m, v in values.items():
        extra = set(v["channel"]) - set(channels)
        if extra:
            notes.append(f"{m}: {len(extra)} channel(s) not among the good channels, ignored")

    ranks = pd.DataFrame({m: rank.rank_method(v, channels, p["MISSING_CHANNELS"]) for m, v in values.items()})
    consensus = rank.fuse(ranks, p["FUSION"])
    cons_rank = pd.Series(rank.rank_keys([(-s if np.isfinite(s) else np.inf,) for s in consensus]), index=channels)

    labels = rec.onset_labels if rec else {}
    has_labels = bool(rec and rec.has_labels)
    n_expert = sum(1 for c in channels if labels.get(c))
    mode, k = p["ZONE_MODE"], p["ZONE_K"]
    if mode == "expert-size":
        if has_labels:
            k = n_expert
        else:
            notes.append("ZONE_MODE=expert-size without labels: top-k with ZONE_K")
            mode = "top-k"
    zone = rank.select(consensus, mode, k, p["ZONE_THRESHOLD"])

    tab = pd.DataFrame({"recording_id": rid, "channel": channels})
    for m, v in values.items():
        vv = v.drop_duplicates("channel").set_index("channel")
        for col in [c for c in vv.columns if not c.startswith("_")]:
            tab[col] = tab["channel"].map(vv[col])
        tab[f"rank_{m}"] = ranks[m].to_numpy()
    tab["n_methods"] = ranks.notna().sum(axis=1).to_numpy()
    tab["consensus_score"] = consensus.to_numpy().round(6)
    tab["rank"] = cons_rank.to_numpy()
    tab["in_onset_zone"] = zone.to_numpy()
    tab["expert_onset"] = [labels.get(c) if labels else None for c in channels]
    for c in ZONE_COLUMNS:
        if c not in tab.columns:
            tab[c] = None
    tab = tab[ZONE_COLUMNS].sort_values(["rank", "channel"], kind="mergesort")

    scores: list[dict] = []

    def add(method, metric, value):
        scores.append({"recording_id": rid, "method": method, "metric": metric,
                       "value": None if value is None or (isinstance(value, float) and not np.isfinite(value)) else value})

    if has_labels:
        known = [c for c in channels if c in labels]
        lab = np.array([labels[c] for c in known], bool)
        all_scores = {m: rank.score_from_rank(ranks[m]) for m in ranks.columns}
        all_scores["consensus"] = consensus
        for m, s in all_scores.items():
            s_known = s.reindex(known)
            add(m, "auc", scoring.auc(s_known.to_numpy(), lab))
            if m == "consensus":
                sel = zone.reindex(known).to_numpy()
            else:
                sel = rank.select(s, mode, k, p["ZONE_THRESHOLD"]).reindex(known).to_numpy()
            pr, rc, f1 = scoring.prf(sel, lab)
            add(m, "precision", pr)
            add(m, "recall", rc)
            add(m, "f1", f1)
            top = rank.select(s, "top-k", n_expert, 0).reindex(known).to_numpy()
            add(m, "r_precision", scoring.prf(top, lab)[0])
            add(m, "n_selected", int(np.asarray(sel, bool).sum()))
        add("expert", "n_onset_channels", int(lab.sum()))
        add("expert", "n_channels_labelled", len(known))

    det, det_src = detected_onset(tables, rid)
    if det is not None:
        add("dynasd", "detected_onset_s", round(det, 3))
        if rec and rec.onset_s is not None:
            add("dynasd", "onset_error_s", round(det - rec.onset_s, 3))
            add("dynasd", "abs_onset_error_s", round(abs(det - rec.onset_s), 3))
            ev = scoring.onset_event(det, rec.onset_s, rec.seizure_end_s, rec.duration_s,
                                     p["EVENT_TOLERANCE_START_S"], p["EVENT_TOLERANCE_END_S"])
            add("dynasd", "event_tp", ev["event_tp"])
            add("dynasd", "event_fp", ev["event_fp"])
        if det_src != "seizure_onsets.csv":
            notes.append(f"detected onset from {det_src}")

    info = {"recording_id": rid, "status": "ok", "channels": len(channels), "methods": list(values),
            "expert_onset_channels": n_expert if has_labels else None,
            "zone_size": int(zone.sum()), "notes": "; ".join(notes)}
    return tab, scores, info


def summaries(scores: pd.DataFrame, recs: dict, splits: dict[str, bool]) -> list[dict]:
    """Rows over recordings: ALL, and per stim_induced / lvfa value when available."""
    groups = {"ALL": list(scores["recording_id"].unique())}
    for flag, available in splits.items():
        if not available:
            continue
        for val in (True, False):
            groups[f"{flag}={str(val).lower()}"] = [r for r, rec in recs.items() if getattr(rec, flag) is val]
    out = []
    for gname, rids in groups.items():
        g = scores[scores["recording_id"].isin(rids)]
        out.append({"recording_id": gname, "method": "all", "metric": "n_recordings",
                    "value": int(g["recording_id"].nunique())})
        for (method, metric), sub in g.groupby(["method", "metric"], sort=True):
            if metric in scoring.SUMMARISED:
                for stat, v in scoring.summarise(list(pd.to_numeric(sub["value"], errors="coerce"))).items():
                    out.append({"recording_id": gname, "method": method, "metric": f"{metric}_{stat}", "value": v})
            elif metric == "event_tp":
                tp = pd.to_numeric(sub["value"], errors="coerce").dropna()
                if len(tp):
                    out.append({"recording_id": gname, "method": method, "metric": "event_sensitivity",
                                "value": float(tp.sum() / len(tp))})
                    out.append({"recording_id": gname, "method": method, "metric": "event_n", "value": int(len(tp))})
    for flag, available in splits.items():
        if not available:
            out.append({"recording_id": f"{flag}=unavailable", "method": "all", "metric": "split_unavailable", "value": 1})
    return out


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S", stream=sys.stdout, force=True)
    started = time.time()
    input_dir = Path(os.environ.get("INPUT_DIR", "/data/input"))
    output_dir = Path(os.environ.get("OUTPUT_DIR", "/data/output"))
    try:
        p = params()
    except ValueError as e:
        log.error("ERROR: %s", e)
        return 1
    log.info("%s %s: %s", APP, version(), " ".join(f"{k}={v}" for k, v in p.items()))

    ieeg_root = inputs.port_root(input_dir, "ieeg")
    methods_root = inputs.port_root(input_dir, "methods")
    index = inputs.read_ieeg(ieeg_root)
    tables = inputs.read_methods(methods_root, input_dir)
    log.info("ieeg: %d recording(s)%s; methods: %s", len(index.recordings),
             f" from {os.path.relpath(index.index_path, input_dir)}" if index.index_path else " (no recordings.csv)",
             ", ".join(f"{k} ({v})" for k, v in tables.paths.items()) or "none")
    if not tables.present:
        log.error("ERROR: no channel_onsets.csv, ei.csv or fragility.csv under %s", methods_root)
        return 1

    rids: dict[str, None] = dict.fromkeys(index.recordings)
    for t in (tables.channel_onsets, tables.ei, tables.fragility, tables.seizure_onsets):
        if t is not None:
            rids.update(dict.fromkeys(sorted(set(t["recording_id"]) - set(rids))))
    if not rids:
        log.error("ERROR: no recordings found")
        return 1

    splits = {"stim_induced": "stim_induced" in index.columns, "lvfa": "lvfa" in index.columns}
    for flag, ok in splits.items():
        if not ok:
            log.info("note: recordings.csv has no %s column: the %s split is unavailable", flag, flag)

    zones, score_rows, infos = [], [], []
    for rid in rids:
        t0 = time.time()
        rec = index.recordings.get(rid)
        try:
            tab, sc, info = process_recording(rid, rec, tables, p)
        except Exception as e:  # one recording failing must not stop the others
            log.info("%s: FAILED: %s", rid, e)
            infos.append({"recording_id": rid, "status": "failed", "notes": str(e)})
            continue
        info["seconds"] = round(time.time() - t0, 3)
        zones.append(tab)
        score_rows += sc
        infos.append(info)
        auc = next((s["value"] for s in sc if s["method"] == "consensus" and s["metric"] == "auc"), None)
        log.info("%s: %d channels, methods %s, zone %d channel(s)%s%s %.2fs", rid, info["channels"],
                 "+".join(info["methods"]), info["zone_size"],
                 f", expert onset {info['expert_onset_channels']}, consensus AUC {auc:.3f}" if auc is not None else "",
                 f" [{info['notes']}]" if info["notes"] else "", info["seconds"])

    ok = [i for i in infos if i["status"] == "ok"]
    output_dir.mkdir(parents=True, exist_ok=True)
    if zones:
        z = pd.concat(zones, ignore_index=True)
        for col in ("in_onset_zone", "dynasd_recruited", "expert_onset"):
            z[col] = z[col].map(lambda v: "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(bool(v)).lower())
        z.to_csv(output_dir / "onset_zone.csv", index=False, float_format="%.6g")
    scores = pd.DataFrame(score_rows, columns=["recording_id", "method", "metric", "value"])
    labelled = {r: index.recordings[r] for r in scores["recording_id"].unique() if r in index.recordings}
    summary = summaries(scores, labelled, splits) if len(scores) else []
    pd.concat([scores, pd.DataFrame(summary, columns=scores.columns)], ignore_index=True).to_csv(
        output_dir / "scores.csv", index=False, float_format="%.6g")

    prov = {
        "app": APP, "version": version(),
        "tools": {"python": sys.version.split()[0],
                  **{pkg: _v(pkg) for pkg in ("numpy", "pandas", "timescoring")}},
        "parameters": p,
        "inputs": {"ieeg": os.path.relpath(index.index_path, input_dir) if index.index_path else None,
                   "channels": sorted(r.channels_file for r in index.recordings.values() if r.channels_file),
                   "methods": tables.paths, "ignored": tables.ignored},
        "splits": {k: ("available" if v else "unavailable: recordings.csv has no " + k + " column")
                   for k, v in splits.items()},
        "started": datetime.fromtimestamp(started, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seconds": round(time.time() - started, 2),
        "recordings": infos,
    }
    (output_dir / "provenance.json").write_text(json.dumps(prov, indent=2, default=str))

    cons_auc = scores[(scores["method"] == "consensus") & (scores["metric"] == "auc")]["value"]
    cons_auc = pd.to_numeric(cons_auc, errors="coerce").dropna()
    log.info("done: %d ok, %d failed%s; %.1fs", len(ok), len(infos) - len(ok),
             f"; consensus AUC median {cons_auc.median():.3f} over {len(cons_auc)} labelled recording(s)" if len(cons_auc) else "",
             time.time() - started)
    if not ok:
        log.error("ERROR: every recording failed")
        return 1
    return 0


def _v(pkg: str) -> str:
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return "unknown"


if __name__ == "__main__":
    sys.exit(main())
