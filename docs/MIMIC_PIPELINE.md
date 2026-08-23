# MIMIC-IV pipeline

The MIMIC-IV 3.1 workflow selected the analysis unit, applied prespecified clinical exclusions, constructed a creatinine baseline and interval-observed event brackets, materialized causal hourly features, selected Elastic Net/XGBoost/Trend vehicles without TEST access, and selected first-alert thresholds on VALIDATION non-events under `score > c`.

After TEST opening, the frozen code computed patient-policy contributions, single-policy identified intervals, direct shared-onset paired intervals, patient bootstrap outer sets, and paired max-statistic simultaneous sets. Public helpers retain closed endpoint equality and the original metric denominators. The public code does not fit the frozen models or contain fitted constants.
