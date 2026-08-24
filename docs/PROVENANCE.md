# Reproducibility provenance

The numerical results reported in the manuscript were generated with the statistical implementation preserved at:

- tag: `v1.0.0-submission`
- commit: `6df6e4b7679a4dcfcd9ddf17b9e5a7f9dab2d097`

Later commits reorganize documentation and repository presentation. They do not alter the statistical implementation used for the reported results.

## Final bootstrap runs

| Analysis | Procedure | Valid replicates |
| --- | --- | ---: |
| MIMIC-IV | Patient bootstrap | 5000 / 5000 |
| eICU primary | Hospital-then-patient bootstrap | 5000 / 5000 |
| eICU conditional sensitivity | Patient bootstrap conditional on the 40 test hospitals | 1000 / 1000 |

Models and alert thresholds were not refitted during bootstrap.

Final bootstrap lock: `a9335f5f768459276a7e2834aefd2b78f7824cd8adab04ed4244a553a6afdedd`

For implementation checks, see [`TEST_CASE_MAPPING.csv`](TEST_CASE_MAPPING.csv).
