# eICU pipeline

These modules provide portable interfaces for authorized local copies of the eICU Collaborative Research Database 2.0. The repository contains no raw records or patient-level derived data.

Copy [`public_config_template.yaml`](public_config_template.yaml) outside the repository and set `EICU_ROOT`, `ONSETBOUNDS_OUTPUT_ROOT`, and the path to locally fitted training constants.

The analysis proceeds through hospital-aware cohort construction, hourly feature construction, independent model and threshold selection, test-cohort evaluation, and bootstrap inference. The public entry points are importable modules:

- [`cohort.py`](cohort.py) for first-unit-stay selection;
- [`features.py`](features.py) for hourly aggregation and externally supplied training medians;
- [`model_selection.py`](model_selection.py) for weighted average precision and the `score > c` threshold rule;
- [`evaluate.py`](evaluate.py) for interval-record evaluation and hospital-then-patient bootstrap resampling.

The primary procedure resamples hospitals and then patients within sampled hospitals. The conditional analysis resamples patients conditional on the 40 test hospitals. These modules are building blocks for an authorized local driver rather than a standalone database extraction command. Data access is described in [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md).
