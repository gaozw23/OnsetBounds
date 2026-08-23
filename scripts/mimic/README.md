# MIMIC-IV pipeline

This directory exposes database-column interfaces, deterministic cohort/feature helpers, frozen metric rules, and evaluation code. It does not contain MIMIC-IV records, fitted models, preprocessing arrays, thresholds, or patient-level derivatives. Authorized users must provide MIMIC-IV 3.1 locally and external TRAIN-derived constants through configuration.

Stage order and frozen decisions are documented in `docs/MIMIC_PIPELINE.md`. Schema names such as `subject_id` are code interfaces, not embedded identifier values.
