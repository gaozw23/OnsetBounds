# MIMIC-IV pipeline

These modules provide portable interfaces for authorized local copies of MIMIC-IV 3.1. The repository contains no raw records or patient-level derived data.

Copy [`public_config_template.yaml`](public_config_template.yaml) outside the repository and set `MIMIC_IV_ROOT`, `ONSETBOUNDS_OUTPUT_ROOT`, and the path to locally fitted training constants.

The analysis proceeds through cohort and interval construction, hourly feature construction, model and threshold selection, test-cohort evaluation, and patient bootstrap inference. The public entry points are importable modules:

- [`cohort.py`](cohort.py) for analysis-stay selection;
- [`features.py`](features.py) for hourly aggregation and externally supplied training medians;
- [`model_selection.py`](model_selection.py) for weighted average precision and the `score > c` threshold rule;
- [`evaluate.py`](evaluate.py) for aggregation of externally generated alert and interval records.

These modules are building blocks for an authorized local driver rather than a standalone database extraction command. Data access is described in [`docs/DATA_ACCESS.md`](../../docs/DATA_ACCESS.md).
