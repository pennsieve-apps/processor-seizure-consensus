"""
Scores against the expert labels: channel-level (AUC, precision/recall/F1, R-precision) for
every method and the consensus, onset-time error of seizure-dynasd, and summaries over
recordings (all, and split by stimulation-induced vs spontaneous and by LVFA when
recordings.csv says which is which).
"""
from __future__ import annotations

import numpy as np
from timescoring.annotations import Annotation
from timescoring.scoring import EventScoring

CHANNEL_METRICS = ("auc", "precision", "recall", "f1", "r_precision")
SUMMARISED = CHANNEL_METRICS + ("onset_error_s", "abs_onset_error_s")


def auc(score: np.ndarray, label: np.ndarray) -> float:
    """Area under the ROC curve: the probability that a random onset channel scores higher
    than a random other channel, ties counting one half (Mann-Whitney U / (n1 n0))."""
    score = np.asarray(score, float)
    label = np.asarray(label, bool)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n1, n0 = int(label.sum()), int((~label).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    s = score[order]
    ranks = np.empty(len(s))
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2 + 1
        i = j + 1
    r = np.empty(len(s))
    r[order] = ranks
    u = r[label].sum() - n1 * (n1 + 1) / 2
    return float(u / (n1 * n0))


def prf(selected: np.ndarray, label: np.ndarray) -> tuple[float, float, float]:
    selected = np.asarray(selected, bool)
    label = np.asarray(label, bool)
    tp = int((selected & label).sum())
    p = tp / selected.sum() if selected.sum() else float("nan")
    r = tp / label.sum() if label.sum() else float("nan")
    if np.isfinite(p) and np.isfinite(r):
        f = 0.0 if p + r == 0 else 2 * p * r / (p + r)
    else:
        f = float("nan")
    return float(p), float(r), float(f)


def onset_event(detected_s: float, ref_onset_s: float, ref_end_s: float | None,
                duration_s: float | None, tol_start: float, tol_end: float) -> dict:
    """SzCORE-style event scoring of one detected onset with timescoring's EventScoring.

    Reference event: [expert onset, seizure end] (end from the events.tsv seizure row; onset
    + 10 s when unknown). Hypothesis: a 1-s event at the detected onset, since seizure-dynasd
    reports an onset, not an end. A detection counts (tp = 1) when it overlaps the reference
    widened by tol_start before and tol_end after (SzCORE: 30 s, 60 s)."""
    ref_end = ref_end_s if (ref_end_s is not None and ref_end_s > ref_onset_s) else ref_onset_s + 10.0
    span = max(duration_s or 0.0, ref_end + tol_end + 1.0, detected_s + 2.0)
    # Events must lie inside the recording; clip at 0 (a detection before the clip start is
    # not possible anyway).
    ref = Annotation([(max(0.0, ref_onset_s), ref_end)], 10, int(np.ceil(span * 10)))
    hyp = Annotation([(max(0.0, detected_s), max(0.0, detected_s) + 1.0)], 10, int(np.ceil(span * 10)))
    params = EventScoring.Parameters(toleranceStart=tol_start, toleranceEnd=tol_end, minOverlap=0,
                                     maxEventDuration=max(5 * 60, ref_end - ref_onset_s + 1),
                                     minDurationBetweenEvents=90)
    sc = EventScoring(ref, hyp, params)
    return {"event_tp": int(sc.tp), "event_fp": int(sc.fp)}


def summarise(values: list[float]) -> dict[str, float]:
    v = np.asarray([x for x in values if x is not None and np.isfinite(x)], float)
    if not len(v):
        return {"n": 0}
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return {"n": int(len(v)), "median": float(med), "q1": float(q1), "q3": float(q3),
            "min": float(v.min()), "max": float(v.max()), "mean": float(v.mean())}
