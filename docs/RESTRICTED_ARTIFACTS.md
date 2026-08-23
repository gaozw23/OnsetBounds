# Deliberately excluded artifacts

The public freeze excludes all raw MIMIC-IV/eICU records and these patient-derived objects:

- cohort rows, outcome brackets, trajectory audits, source-row extracts, and exact timestamps;
- canonical hourly matrices and score rows;
- subject, admission, stay, patient-unit-stay, unique-patient, raw-hospital, and study identifiers;
- fitted model files and complete coefficient vectors;
- TRAIN-derived imputation/scaling arrays and patient-specific transformations;
- patient-policy and paired contribution tables;
- bootstrap draw/checkpoint files and TEST-open records containing restricted paths;
- virtual environments, logs, transcripts, caches, temporary renders, and LaTeX auxiliaries.

These exclusions prevent redistribution of controlled data and reduce re-identification and model-inversion risk. Public CSVs are aggregate manuscript sources listed in `paper_repro/aggregate/PUBLIC_AGGREGATE_DATA_MANIFEST.csv`.
