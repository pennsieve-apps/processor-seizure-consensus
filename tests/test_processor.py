"""The synthetic fixture through the processor (both input layouts), the ranking and scoring
helpers, and the Lambda handler."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from processor import handler, main as app, rank, scoring
from tests.make_fixture import BAD, CHANNELS, RECORDINGS, make_fixture

LABELLED = [r[0] for r in RECORDINGS if r[2] is not None]
UNLABELLED = [r[0] for r in RECORDINGS if r[2] is None]


@pytest.fixture(scope="session")
def roots(tmp_path_factory) -> dict[str, Path]:
    return make_fixture(tmp_path_factory.mktemp("fixture"), "both")


def run(input_dir: Path, out: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(out))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    return app.main()


def zone(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / "onset_zone.csv", dtype={"in_onset_zone": str, "expert_onset": str})


def scores(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / "scores.csv")


def value(sc: pd.DataFrame, rid: str, method: str, metric: str):
    v = sc[(sc.recording_id == rid) & (sc.method == method) & (sc.metric == metric)]["value"]
    assert len(v) == 1, (rid, method, metric)
    return v.iloc[0]


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_both_layouts(roots, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(roots[layout], out, monkeypatch) == 0
    z = zone(out)
    assert list(z.columns) == app.ZONE_COLUMNS
    assert {"recording_id", "channel", "rank", "in_onset_zone"} <= set(z.columns)  # conventions.md
    assert set(z.recording_id) == {r[0] for r in RECORDINGS}
    # the bad channel is left out
    assert not ((z.recording_id == BAD[0]) & (z.channel == BAD[1])).any()
    for rid, _, onset_ch, _, _ in RECORDINGS:
        zr = z[z.recording_id == rid]
        assert (zr.in_onset_zone == "true").sum() == 10  # ZONE_K default
        if onset_ch is not None:
            top3 = set(zr.nsmallest(3, "rank").channel)
            assert top3 == {CHANNELS[i] for i in onset_ch}
    sc = scores(out)
    assert list(sc.columns) == ["recording_id", "method", "metric", "value"]
    for rid in LABELLED:
        for m in ("dynasd", "ei", "fragility", "consensus"):
            assert value(sc, rid, m, "auc") > 0.9
        assert value(sc, rid, "consensus", "recall") == 1
        assert value(sc, rid, "consensus", "precision") == pytest.approx(0.3)
        assert value(sc, rid, "consensus", "r_precision") == 1
        assert value(sc, rid, "dynasd", "onset_error_s") == pytest.approx(1.5)
        assert value(sc, rid, "dynasd", "event_tp") == 1
    for rid in UNLABELLED:
        assert sc[(sc.recording_id == rid) & (sc.metric == "auc")].empty
    assert value(sc, "ALL", "consensus", "auc_n") == 2
    assert value(sc, "ALL", "dynasd", "event_sensitivity") == 1
    # stim/LVFA split present in this fixture
    assert value(sc, "stim_induced=true", "all", "n_recordings") == 1
    assert value(sc, "lvfa=false", "consensus", "auc_n") == 1
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["app"] == "processor-seizure-consensus" and prov["version"] == "0.1.0"
    assert prov["splits"] == {"stim_induced": "available", "lvfa": "available"}
    assert len(prov["inputs"]["methods"]) == 4
    assert all(r["status"] == "ok" for r in prov["recordings"])


def test_no_flags_reports_gap(tmp_path, monkeypatch):
    root = make_fixture(tmp_path / "fx", "linked", flags=False)["linked"]
    out = tmp_path / "out"
    assert run(root, out, monkeypatch) == 0
    sc = scores(out)
    assert value(sc, "stim_induced=unavailable", "all", "split_unavailable") == 1
    assert value(sc, "lvfa=unavailable", "all", "split_unavailable") == 1
    assert not sc.recording_id.str.startswith("stim_induced=true").any()
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["splits"]["lvfa"].startswith("unavailable")


def test_single_method_and_none(roots, tmp_path, monkeypatch):
    src = roots["linked"]
    one = tmp_path / "one"
    shutil.copytree(src / "ieeg", one / "ieeg")
    shutil.copytree(src / "methods" / "fragility", one / "methods" / "fragility")
    assert run(one, tmp_path / "out1", monkeypatch) == 0
    z = zone(tmp_path / "out1")
    assert z.rank_dynasd.isna().all() and z.rank_fragility.notna().all()
    assert (z.n_methods == 1).all()
    none = tmp_path / "none"
    shutil.copytree(src / "ieeg", none / "ieeg")
    assert run(none, tmp_path / "out2", monkeypatch) == 1


def test_missing_channels_and_fusion(roots, tmp_path, monkeypatch):
    rid = RECORDINGS[1][0]  # ei.csv lacks LB07-LB08, LB08-LB09 here
    assert run(roots["linked"], tmp_path / "w", monkeypatch) == 0
    w = zone(tmp_path / "w").query("recording_id == @rid").set_index("channel")
    assert w.loc["LB07-LB08", "rank_ei"] == pytest.approx(15.5)  # shared last places
    assert run(roots["linked"], tmp_path / "i", monkeypatch, MISSING_CHANNELS="ignore") == 0
    i = zone(tmp_path / "i").query("recording_id == @rid").set_index("channel")
    assert np.isnan(i.loc["LB07-LB08", "rank_ei"]) and i.loc["LB07-LB08", "n_methods"] == 2
    assert run(roots["linked"], tmp_path / "b", monkeypatch, FUSION="borda") == 0
    b = zone(tmp_path / "b")
    m = zone(tmp_path / "w")
    # with every channel ranked by every method (MISSING_CHANNELS=worst) Borda and mean rank agree
    assert (b["rank"].to_numpy() == m["rank"].to_numpy()).all()


def test_zone_modes(roots, tmp_path, monkeypatch):
    assert run(roots["linked"], tmp_path / "e", monkeypatch, ZONE_MODE="expert-size", ZONE_K=4) == 0
    z = zone(tmp_path / "e")
    for rid in LABELLED:
        assert (z[z.recording_id == rid].in_onset_zone == "true").sum() == 3
    assert (z[z.recording_id == UNLABELLED[0]].in_onset_zone == "true").sum() == 4
    sc = scores(tmp_path / "e")
    assert value(sc, LABELLED[0], "consensus", "f1") == 1
    assert run(roots["linked"], tmp_path / "t", monkeypatch, ZONE_MODE="threshold", ZONE_THRESHOLD=0.85) == 0
    z = zone(tmp_path / "t")
    assert ((z.in_onset_zone == "true") == (z.consensus_score >= 0.85)).all()
    assert run(roots["linked"], tmp_path / "x", monkeypatch, ZONE_MODE="nonsense") == 1


def test_fragility_difference(roots, tmp_path, monkeypatch):
    assert run(roots["linked"], tmp_path / "d", monkeypatch, FRAGILITY_SCORE="difference") == 0
    z = zone(tmp_path / "d")
    fr = pd.read_csv(roots["linked"] / "methods" / "fragility" / "fragility.csv")
    m = z.merge(fr, on=["recording_id", "channel"])
    assert np.allclose(m.fragility, m.fragility_post - m.fragility_pre, atol=1e-5)


def test_duplicate_tables_ignored(roots, tmp_path, monkeypatch):
    root = tmp_path / "dup"
    shutil.copytree(roots["linked"], root)
    (root / "methods" / "ei" / "rerun").mkdir()
    shutil.copy(root / "methods" / "ei" / "ei.csv", root / "methods" / "ei" / "rerun" / "ei.csv")
    assert run(root, tmp_path / "out", monkeypatch) == 0
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    assert prov["inputs"]["methods"]["ei.csv"] == "methods/ei/ei.csv"
    assert prov["inputs"]["ignored"] == ["methods/ei/rerun/ei.csv"]


def test_dynasd_ranking():
    t = pd.DataFrame({"recording_id": "r", "channel": ["a", "b", "c", "d", "e"],
                      "onset_s": ["12", "", "10", "10", ""], "recruited": ["true", "false", "true", "true", "false"],
                      "score": ["0.9", "0.5", "0.1", "0.8", "0.7"]})
    v = rank.method_values("dynasd", t, "r", "post")
    r = rank.rank_method(v, list("abcde"), "worst")
    assert list(r) == [3, 5, 2, 1, 4]  # recruited by onset (ties by score), then by score


def test_auc_matches_sklearn():
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(1)
    s = rng.integers(0, 5, 60).astype(float)  # many ties
    y = rng.random(60) < 0.3
    assert scoring.auc(s, y) == pytest.approx(roc_auc_score(y, s))
    assert np.isnan(scoring.auc(s, np.zeros(60, bool)))


def test_onset_event_tolerances():
    assert scoring.onset_event(95.0, 120.0, 160.0, 300.0, 30, 60)["event_tp"] == 1   # 25 s early
    assert scoring.onset_event(80.0, 120.0, 160.0, 300.0, 30, 60)["event_tp"] == 0   # 40 s early
    assert scoring.onset_event(215.0, 120.0, 160.0, 300.0, 30, 60)["event_tp"] == 1  # within end + 60 s
    late = scoring.onset_event(230.0, 120.0, 160.0, 300.0, 30, 60)
    assert late == {"event_tp": 0, "event_fp": 1}


def test_lambda_handler(roots, tmp_path, monkeypatch):
    monkeypatch.setenv("ZONE_K", "3")  # stale value from a previous warm invocation
    out = tmp_path / "lam"
    res = handler.handler({"inputDir": str(roots["linked"]), "outputDir": str(out), "FUSION": "borda"}, None)
    assert res["status"] == "ok"
    z = zone(out)
    assert (z.groupby("recording_id").in_onset_zone.apply(lambda s: (s == "true").sum()) == 10).all()
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["parameters"]["FUSION"] == "borda" and prov["parameters"]["ZONE_K"] == 10
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "empty"), "outputDir": str(tmp_path / "o2")}, None)
