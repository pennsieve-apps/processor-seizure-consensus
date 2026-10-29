# processor-seizure-consensus

The fusion step of the **Seizure Onset Localisation (SEEG)** workflow:
`ieeg-prep → seizure-dynasd → {seizure-ei, seizure-fragility} → seizure-consensus → ieeg-report`

Each onset-localisation method looks at a seizure from a different angle:

- **seizure-dynasd** (DynaSD; Ojemann et al. 2026): when each channel is recruited into the
  seizure, from how far its dynamics diverge from the pre-ictal baseline;
- **seizure-ei** (Epileptogenicity Index; Bartolomei et al. 2008): how early and how strongly
  a channel shows fast (beta/gamma) activity at onset;
- **seizure-fragility** (neural fragility; Li et al. 2021): how close a channel is to
  destabilising the network, from linear models fitted to short windows.

None of them is right on every seizure. This step puts them on one scale, so a channel ranked
high by all three methods stands out from one ranked high by only one.

- **For an epileptologist:** `onset_zone.csv` lists, for every seizure recording, the
  channels ranked by the consensus, with each method's own value and rank beside it, and
  marks the top channels as the estimated onset zone. It is a summary of quantitative
  markers to set beside visual review, not a diagnosis.
- **For a researcher:** when the recordings carry expert onset labels (the `onset` column of
  `_channels.tsv`), `scores.csv` scores every method and the consensus against them: AUC,
  precision/recall/F1 of the onset zone, R-precision, the error of the detected seizure onset
  time and SzCORE event scoring (Dan et al. 2024) of that onset, per recording and summarised
  over recordings, also split by stimulation-induced vs spontaneous seizures and by low-voltage
  fast activity (LVFA) when recordings.csv says which is which.

## Ports

| Port | Direction | Placement | What |
|---|---|---|---|
| `ieeg` | input (required) | `INPUT_DIR/ieeg/` | ieeg-prep's BIDS-iEEG derivative (or a plain BIDS-iEEG dataset): `recordings.csv`, `<stem>_desc-preproc_channels.tsv` (good channels, `onset` labels), `<stem>_desc-preproc_events.tsv` (`seizure` row). Signals are not read. |
| `methods` | input (required, fan-in) | `INPUT_DIR/methods/` | seizure-dynasd (`seizure_onsets.csv`, `channel_onsets.csv`), seizure-ei (`ei.csv`), seizure-fragility (`fragility.csv`). Found by file name anywhere under the port; any one of `channel_onsets.csv`, `ei.csv`, `fragility.csv` is enough. |
| `onset-zone` | output | `OUTPUT_DIR/` | `onset_zone.csv`, `scores.csv`, `provenance.json`. |

Without links (everything side by side in `INPUT_DIR`) the files are found by name anywhere
under `INPUT_DIR`; symlinks are followed. When two tables of the same kind are found (e.g. a
step run twice under different `LABEL` folders), the shallowest one is used and the others are
listed in `provenance.json → inputs.ignored`.

## Parameters

| Name | Default | What it does |
|---|---|---|
| `FUSION` | `mean-rank` | `mean-rank`: consensus score = 1 − mean over methods of the normalised rank (r − 1)/(n − 1). `borda`: each method gives n − r points; score = points / (methods × (n − 1)). With every channel ranked by every method (`MISSING_CHANNELS=worst`) both give the same order; they differ with `MISSING_CHANNELS=ignore`, where Borda gives a missing channel 0 points from that method and mean rank averages over the methods that ranked it. |
| `FRAGILITY_SCORE` | `post` | Ranks fragility by `fragility_post` (mean fragility after onset) or by `difference` (`fragility_post − fragility_pre`), which favours channels that became fragile at onset over channels that always are. |
| `MISSING_CHANNELS` | `worst` | A good channel a method did not score: `worst` = it shares that method's last places (average rank); `ignore` = that method is left out of its consensus (`n_methods` shows how many counted). |
| `ZONE_MODE` | `top-k` | `top-k`: the `ZONE_K` best channels (exactly k; ties broken by channel order). `threshold`: every channel with `consensus_score ≥ ZONE_THRESHOLD`. `expert-size`: as many channels as the expert marked as onset — **this uses the labels, so it is for scoring methods only**, never for unlabelled clinical data; a recording without labels falls back to `top-k`. |
| `ZONE_K` | `10` | Zone size for `top-k` (~5 % of the ~200 bipolar channels of a demo clip). |
| `ZONE_THRESHOLD` | `0.9` | Score threshold for `threshold` (0–1; 1 = ranked first by every method). |
| `EVENT_TOLERANCE_START_S` | `30` | SzCORE event scoring: a detected onset up to this long before the expert onset still counts. |
| `EVENT_TOLERANCE_END_S` | `60` | SzCORE event scoring: a detected onset up to this long after the end of the seizure still counts. |

## How each method ranks channels

Ranks run from 1 (most likely onset) over the recording's **good** channels (`status` not
`bad` in `_channels.tsv`; without a `_channels.tsv`, the channels the methods report). Ties
share their average rank. Channels a method reports but that are not good channels are ignored
(logged).

| Method | Table | Ranked by |
|---|---|---|
| dynasd | `channel_onsets.csv` | Recruited channels first (`recruited` true; a channel with an `onset_s` counts as recruited when the column is missing), by `onset_s` (earliest first), ties by `score` (highest first); then the channels that were not recruited, by `score`. |
| ei | `ei.csv` | `ei`, highest first. |
| fragility | `fragility.csv` | `fragility_post` (or post − pre, see `FRAGILITY_SCORE`), highest first. |

The step's own `rank` column is the rank of the consensus score.

## Outputs

**`onset_zone.csv`**: one row per recording and good channel, sorted by rank.

| Column | |
|---|---|
| `recording_id`, `channel` | Keys (conventions.md) |
| `dynasd_onset_s`, `dynasd_score`, `dynasd_recruited` | seizure-dynasd's values (empty when that method was not given) |
| `ei` | seizure-ei's EI |
| `fragility` | the fragility value used for ranking (`fragility_post` or the difference) |
| `rank_dynasd`, `rank_ei`, `rank_fragility` | per-method ranks (1 = best) |
| `n_methods` | methods that ranked this channel |
| `consensus_score` | 0–1, higher = more likely onset |
| `rank` | consensus rank |
| `in_onset_zone` | `true` for the selected channels (`ZONE_MODE`) |
| `expert_onset` | the expert label from `_channels.tsv` (empty when unknown), for convenience |

**`scores.csv`** (`recording_id, method, metric, value`): written always; empty of channel
metrics when no recording has labels. `method` is `dynasd`, `ei`, `fragility`, `consensus`,
`expert` (label counts) or `all`.

Per recording with labels (a recording counts as labelled when at least one good channel has
`onset=true`; channels with an empty label are left out of its scores):

| Metric | |
|---|---|
| `auc` | ROC AUC of the method's ranking against the onset labels (Mann-Whitney, ties ½) |
| `precision`, `recall`, `f1` | of the method's own selection under `ZONE_MODE` (each method selects by its own rank; `consensus` uses `in_onset_zone`) |
| `r_precision` | precision of the top-N channels, N = number of expert onset channels (= recall = F1 at that size): the "expert-size" score, independent of `ZONE_K` |
| `n_selected` | channels selected |
| `expert/n_onset_channels`, `expert/n_channels_labelled` | |

Per recording with a dynasd onset (`seizure_onsets.csv`, earliest row; else the earliest
recruited channel of `channel_onsets.csv`, noted in provenance) and an expert onset
(recordings.csv `onset_s`, else the `seizure` event):

| Metric | |
|---|---|
| `dynasd/detected_onset_s` | |
| `dynasd/onset_error_s`, `abs_onset_error_s` | detected − expert (negative = early) |
| `dynasd/event_tp`, `event_fp` | SzCORE event scoring with timescoring's `EventScoring`: reference = [expert onset, seizure end] (end from the `seizure` event's duration, onset + 10 s when unknown), hypothesis = a 1-s event at the detected onset (dynasd reports an onset, not an end). `event_tp = 1` when the detection falls between onset − `EVENT_TOLERANCE_START_S` and end + `EVENT_TOLERANCE_END_S`. |

Summaries, with `recording_id` = `ALL`, and `stim_induced=true|false`, `lvfa=true|false` when
recordings.csv has those columns: `n_recordings`, and for `auc, precision, recall, f1,
r_precision, onset_error_s, abs_onset_error_s` the `<metric>_n, _median, _q1, _q3, _min,
_max, _mean` over recordings; `event_sensitivity` (fraction of recordings whose onset was
detected within tolerance) and `event_n`. A false-alarm rate is not reported: each clip holds
one seizure and dynasd reports one onset, so it would be meaningless.

**Stimulation / LVFA split: a gap in the contract.** conventions.md does not define where the
seizure type comes from. This app reads two optional recordings.csv columns, `stim_induced`
and `lvfa` (`true`/`false`), which ieeg-prep would fill from `labels.csv` recording-level rows
(the demo's `annotations.tsv` has `stim` and `lvfa` per seizure). When a column is missing,
the split is not made, `scores.csv` gets a row `stim_induced=unavailable, all,
split_unavailable, 1` (likewise for `lvfa`), the log says so and `provenance.json → splits`
records it. The app never reads `annotations.tsv` itself.

**`provenance.json`**: app and version, tool versions (Python, NumPy, pandas, timescoring),
parameters, inputs (recordings.csv, channels files, method tables used and ignored), splits,
timings and per-recording status/notes.

A recording that fails is logged and skipped; the step fails when every recording failed, or
when no `channel_onsets.csv`, `ei.csv` or `fragility.csv` was found.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2 GiB, 15 min timeout. The step reads only small
tables; time and memory grow linearly with recordings × channels.

Measured in Docker (linux/amd64 image under emulation on an M-series Mac, `--memory 1g`) on a
stand-in for ieeg-prep's output built from the demo data (3 clips of HUP224, 208 bipolar
channels each, expert onset labels) with label-blind *random* method tables of real size
(the upstream method steps had not produced output yet): **0.3 s** for the processing, 2.6 s
wall time including container start; cgroup `memory.peak` **122 MB**. The image is 294 MB
unpacked, ~69 MB compressed. 2 GiB (and the Lambda's 2048 MB) leaves ample room; 24 or
several hundred recordings stay within seconds.

## Running it locally

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python        # unit tests on the synthetic fixture (no Docker)
make run PYTHON=.venv/bin/python         # Docker over the fixture -> data/output/
make run-ieeg IN=/path/with/ieeg-and-methods   # Docker over real upstream outputs
```

The tests run the fixture through both input layouts (linked ports and everything side by
side), check the ranking, fusion, missing-channel policies, zone modes, the AUC against
scikit-learn, the SzCORE tolerances, the stim/LVFA split and its absence, and the Lambda
handler (parameters from the event, stale environment cleared, failure raised).

## Limitations

- The consensus is only as good as its inputs: three methods that share a failure mode (e.g.
  an artefact-laden channel that looks "early" and "fragile") will agree on it.
- Ranks discard how much better one channel is than the next; a method that is confident on
  one recording and guessing on another weighs the same in both.
- `in_onset_zone` with `top-k` always selects k channels, even when the seizure has a focal
  onset in two channels or a diffuse one in thirty; `threshold` adapts to agreement between
  methods, not to the true zone size. `expert-size` needs labels and is for evaluation only.
- AUC and F1 are computed against per-seizure expert onset channels (the unequivocal
  electrographic onset in the demo), not against resection or outcome; agreement with
  experts is not proof of the epileptogenic zone.
- Bipolar channels carry a label when either contact has it (conventions.md), which inflates
  the expert set slightly at its edges.
- Onset-time error uses dynasd only; EI and fragility do not report a seizure onset here.

## Notes on the conventions

- **Fan-in collision of `provenance.json`.** conventions.md asks every step to write
  `provenance.json` into each output folder, and says fan-in entries land side by side and two
  entries with the same name fail the run. Linking dynasd's `onsets`, ei's and fragility's
  outputs into `methods` would then collide on `provenance.json`. This app finds the tables
  anywhere under the port, so it works whether the platform lands each link as its own folder
  or the method steps write into a subfolder; the contract needs one of those.
- recordings.csv has no `stim_induced` / `lvfa` columns (see above).

## Licences

Apache-2.0 for this app's code (`LICENSE`). The image bundles timescoring (MIT), NumPy and
pandas (BSD-3-Clause), awslambdaric (Apache-2.0) and their dependencies; no copyleft code.
See `THIRD_PARTY.md`.

## Citations

Scoring approach (SzCORE event scoring, with tolerances of 30 s before and 60 s after):

> Dan J, Pale U, Amirshahi A, Cappelletti W, Ingolfsson TM, Wang X, Cossettini A, Bernini A,
> Benini L, Beniczky S, Atienza D, Ryvlin P. **SzCORE: Seizure Community Open-Source Research
> Evaluation framework for the validation of electroencephalography-based automated seizure
> detection algorithms.** Epilepsia 66(S3):14–24 (2025; online 2024). doi:10.1111/epi.18113,
> PMID 39292446.

> Dan J, Pale U. **timescoring** 0.0.7: library for measuring performance of time series
> classification. https://github.com/esl-epfl/timescoring (MIT). No DOI.

The methods whose rankings are fused:

> Ojemann WKS, Xu Z, Shi H, Walsh K, Pattnaik A, Sinha N, Lavelle S, Aguila C, Gallagher R,
> Revell A, LaRocque JJ, Korzun J, Kulick-Soper CV, Zhou DJ, Galer PD, Sinha SR, Shinohara RT,
> Davis KA, Litt B, Conrad EC. **Unsupervised seizure annotation and detection with neural
> dynamic divergence.** medRxiv (2026), preprint. doi:10.64898/2026.02.15.26346325.

> Bartolomei F, Chauvel P, Wendling F. **Epileptogenicity of brain structures in human
> temporal lobe epilepsy: a quantified study from intracerebral EEG.** Brain 131(7):1818–1830
> (2008). doi:10.1093/brain/awn111, PMID 18556663.

> Li A, Huynh C, Fitzgerald Z, Cajigas I, Brusko D, Jagid J, Claudio AO, Kanner AM, Hopp J,
> Chen S, Haagensen J, Johnson E, Anderson W, Crone N, Inati S, Zaghloul KA, Bulacio J,
> Gonzalez-Martinez J, Sarma SV. **Neural fragility as an EEG marker of the seizure onset
> zone.** Nat Neurosci 24(10):1465–1474 (2021). doi:10.1038/s41593-021-00901-w, PMID 34354282.

Demo data:

> Ojemann WK, Zhou D, Armstrong C, Sinha N, Litt B, Conrad E. **Multi-Expert Seizure
> Annotation.** Pennsieve Discover (2026). doi:10.26275/dely-fzwx (CC BY 4.0).

Journal DOIs were checked against Crossref (title, authors, year), PMIDs against PubMed
E-utilities, the dataset DOI against DataCite.
