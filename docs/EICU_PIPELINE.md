# eICU pipeline

The eICU v2.0 replication used hospital-aware cohort splitting, database-specific feature semantics, and independent model fitting; no MIMIC coefficient transport was reported. Threshold selection used VALIDATION non-event patient maxima and strict `score > c` with the higher threshold chosen when target-distance ties occurred.

Primary sampling inference used a two-stage hospital/patient bootstrap. A separate conditional-on-observed-hospitals bootstrap remained at 1000 replicates. The public implementation includes the two-stage resampling primitive, but no hospital identifiers, frozen patient contributions, draws, models, preprocessing arrays, or thresholds.
