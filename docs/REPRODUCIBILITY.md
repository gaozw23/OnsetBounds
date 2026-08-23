# Reproducibility

## Public, no-data workflow

1. Use Python 3.12 and install `.[simulation,dev]`.
2. Run `python -m pytest -q`.
3. Run the synthetic example.
4. Rebuild the five figures from `paper_repro/aggregate/figure_data`.

These steps do not access MIMIC-IV or eICU and do not reproduce patient-level analyses.

## Authorized real-data workflow

The frozen stage order was: cohort and interval construction; hourly features using TRAIN-only constants; model selection on TRAIN/VALIDATION; threshold selection on VALIDATION non-events using strict `score > c`; TEST opening; patient-level contributions and point bounds; bootstrap inference; aggregate reporting. Exact locks are listed in `ANALYSIS_LOCKS.md`.

The public pipeline modules expose statistical rules and portable interfaces. They intentionally omit fitted models, coefficient vectors, imputation/scaling arrays, thresholds, TEST records, and bootstrap draws. Those objects must not be inferred from public aggregate output.

The frozen environment used Python 3.12.13, NumPy 2.5.2, SciPy 1.18.0, pandas 2.3.3, scikit-learn 1.9.0, XGBoost 3.3.0, Matplotlib 3.11.1, and PyYAML 6.0.3. Final bootstrap counts were MIMIC patient 5000/5000, eICU two-stage hospital/patient 5000/5000, and eICU conditional-on-observed-hospitals 1000/1000.
