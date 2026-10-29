# sample-input

This app consumes the *outputs* of other steps of the Seizure Onset Localisation workflow:
ieeg-prep's `ieeg` derivative (for the channels and the expert labels) and the tables of
seizure-dynasd, seizure-ei and seizure-fragility. There is nothing to download for it directly.

For a local run, generate the synthetic fixture:

```sh
make fixture      # writes data/fixture/linked/ and data/fixture/flat/ (gitignored)
make run          # runs the Docker image over data/fixture/linked/
```

`tests/make_fixture.py` writes three recordings with 16 bipolar SEEG channels: recordings.csv,
`_channels.tsv` with `onset` labels (one recording unlabelled, one channel bad), `_events.tsv`
with the `seizure` row, and the four method tables, with the expert onset channels scored
higher plus noise.

For real data, run the upstream steps on the Multi-Expert Seizure Annotation dataset
(Pennsieve Discover 530, doi:10.26275/dely-fzwx, CC BY 4.0; `scripts/fetch-sample-data.sh` in
workflow-epilepsy-ieeg downloads it), then put their outputs in one folder:

```
data/real/ieeg/      <- ieeg-prep's ieeg output (recordings.csv, sub-*/…)
data/real/methods/   <- seizure-dynasd's onsets, seizure-ei's and seizure-fragility's outputs
```

and `make run-ieeg IN=data/real`.
