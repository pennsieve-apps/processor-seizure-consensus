"""
Per-method channel ranking, rank fusion and the onset-zone selection.

Ranks run 1 = most likely onset channel; ties share the average rank. Every method is ranked
over the same channel set (the recording's good channels), so ranks are comparable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from processor.inputs import to_bool

FUSIONS = ("mean-rank", "borda")
ZONE_MODES = ("top-k", "threshold", "expert-size")
MISSING = ("worst", "ignore")


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def method_values(method: str, table: pd.DataFrame, recording_id: str, fragility_score: str) -> pd.DataFrame:
    """The method's rows for one recording: channel, the values shown in onset_zone.csv and a
    sort key (smaller = more likely onset), NaN where the method gave nothing usable."""
    t = table[table["recording_id"] == recording_id].drop_duplicates("channel", keep="first")
    out = pd.DataFrame({"channel": t["channel"].to_numpy()})
    if method == "dynasd":
        onset = _num(t["onset_s"]) if "onset_s" in t.columns else pd.Series(np.nan, index=t.index)
        score = _num(t["score"]) if "score" in t.columns else pd.Series(np.nan, index=t.index)
        if "recruited" in t.columns:
            recruited = t["recruited"].map(to_bool)
            recruited = recruited.where(recruited.notna(), onset.notna())
        else:
            recruited = onset.notna()
        recruited = recruited.astype(bool).to_numpy()
        out["dynasd_onset_s"] = onset.to_numpy()
        out["dynasd_score"] = score.to_numpy()
        out["dynasd_recruited"] = recruited
        # Recruited channels first, by onset time (earliest first), ties by score (highest
        # first); then the rest by score. Encoded as one lexicographic key.
        on = np.where(recruited & np.isfinite(onset.to_numpy()), onset.to_numpy(), np.inf)
        sc = np.where(np.isfinite(score.to_numpy()), -score.to_numpy(), np.inf)
        out["_key"] = list(zip((~recruited).astype(int), on, sc))
        out["_usable"] = recruited | np.isfinite(score.to_numpy())
    elif method == "ei":
        v = _num(t["ei"]).to_numpy()
        out["ei"] = v
        out["_key"] = list(zip(np.where(np.isfinite(v), -v, np.inf)))
        out["_usable"] = np.isfinite(v)
    elif method == "fragility":
        post = _num(t["fragility_post"]).to_numpy()
        pre = _num(t["fragility_pre"]).to_numpy() if "fragility_pre" in t.columns else np.full(len(t), np.nan)
        v = post - pre if fragility_score == "difference" else post
        out["fragility"] = v
        out["_key"] = list(zip(np.where(np.isfinite(v), -v, np.inf)))
        out["_usable"] = np.isfinite(v)
    else:  # pragma: no cover
        raise ValueError(method)
    return out


def rank_keys(keys: list[tuple]) -> np.ndarray:
    """Average ranks (1 = smallest key) of tuple keys; equal keys share their mean rank."""
    order = sorted(range(len(keys)), key=lambda i: keys[i])
    ranks = np.empty(len(keys))
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and keys[order[j + 1]] == keys[order[i]]:
            j += 1
        ranks[[order[k] for k in range(i, j + 1)]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def rank_method(values: pd.DataFrame, channels: list[str], missing: str) -> pd.Series:
    """Rank of every channel (index) for one method. Channels the method did not score:
    `worst` gives them the shared last places; `ignore` leaves them NaN."""
    v = values[values["channel"].isin(channels) & values["_usable"]]
    r = pd.Series(np.nan, index=channels, dtype=float)
    if v.empty:
        return r
    r.loc[v["channel"].to_numpy()] = rank_keys(list(v["_key"]))
    lacking = r.isna()
    if missing == "worst" and lacking.any():
        n_have = int((~lacking).sum())
        r[lacking] = n_have + (int(lacking.sum()) + 1) / 2
    return r


def fuse(ranks: pd.DataFrame, fusion: str) -> pd.Series:
    """Consensus score in [0, 1] (1 = most likely onset) from per-method ranks (columns).

    mean-rank: 1 - mean over methods of the normalised rank (r - 1) / (n - 1).
    borda: each method gives n - r points; score = points / (methods x (n - 1)). A channel
    a method did not rank gets 0 points from it (with MISSING=ignore) or the points of
    the shared last place (MISSING=worst).
    """
    n = len(ranks)
    if n <= 1:
        return pd.Series(1.0, index=ranks.index)
    if fusion == "mean-rank":
        norm = (ranks - 1) / (n - 1)
        return 1 - norm.mean(axis=1, skipna=True)
    if fusion == "borda":
        points = (n - ranks).fillna(0)
        return points.sum(axis=1) / (ranks.shape[1] * (n - 1))
    raise ValueError(fusion)


def score_from_rank(rank: pd.Series) -> pd.Series:
    n = rank.notna().sum()
    return 1 - (rank - 1) / max(n - 1, 1)


def select(score: pd.Series, mode: str, k: int, threshold: float) -> pd.Series:
    """Which channels form the onset zone: the k best (ties broken by channel order, so
    exactly k) or every channel with score >= threshold. NaN scores are never selected."""
    s = score.dropna()
    sel = pd.Series(False, index=score.index)
    if mode == "threshold":
        sel[s.index[s >= threshold - 1e-12]] = True
        return sel
    k = max(0, min(int(k), len(s)))
    order = sorted(range(len(s)), key=lambda i: -s.iloc[i])  # stable: file order on ties
    sel[s.index[order[:k]]] = True
    return sel
