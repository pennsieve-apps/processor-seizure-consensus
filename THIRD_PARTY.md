# Third-party software

The code in this repository (`processor/`, `tests/`, `Dockerfile`, …) is Apache-2.0 (see
`LICENSE`). The image bundles the components below, installed unmodified from PyPI as binary
wheels (`requirements.txt`, exact versions). None of them is copyleft; each package's licence
file is inside the image under `/usr/local/lib/python3.12/site-packages/<package>-<version>.dist-info/`.

| Component | Version | What it does here | Licence | Source |
|---|---|---|---|---|
| timescoring | 0.0.7 | SzCORE event scoring of the detected seizure onset (`EventScoring`) | MIT (LICENSE in the repository, © 2023 ESL-EPFL; the wheel's metadata leaves the field empty) | https://github.com/esl-epfl/timescoring |
| NumPy | 2.5.4 | Arrays, ranks, AUC | BSD-3-Clause (bundled parts 0BSD, MIT, Zlib, CC0-1.0) | https://numpy.org |
| pandas | 3.1.0 | Reading and writing the tables | BSD-3-Clause | https://pandas.pydata.org |
| python-dateutil | 3.0.0 | pandas dependency | Apache-2.0 OR BSD-3-Clause | PyPI |
| six | 1.17.0 | python-dateutil dependency | MIT | PyPI |
| tzdata | 2026.6 | pandas dependency (time zones) | Apache-2.0 | PyPI |
| awslambdaric | 4.2.3 | AWS Lambda runtime interface client (the `lambda` compute type) | Apache-2.0 | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| simplejson, snapshot-restore-py | 4.2.0, 1.0.0 | awslambdaric dependencies | MIT OR AFL-2.1; Apache-2.0 | PyPI |
| Python | 3.12.15 | Interpreter | PSF-2.0 | https://www.python.org |
| python:3.12.15-slim | sha256:a0f57064… | Base image (Debian slim + Python) | PSF-2.0; Debian packages under their own free licences | https://hub.docker.com/_/python |

Development only (not in the image): pytest (MIT), scikit-learn (BSD-3-Clause; the tests
check the in-house AUC against `roc_auc_score`).

## Data

The demo data the app was tried on (not in the image or the repository): Ojemann WK, Zhou D,
Armstrong C, Sinha N, Litt B, Conrad E. *Multi-Expert Seizure Annotation.* Pennsieve Discover
(2026). doi:10.26275/dely-fzwx. CC BY 4.0.
