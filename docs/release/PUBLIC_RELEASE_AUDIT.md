# OnsetBounds public-release status

STATUS = PUBLIC_REPOSITORY_REMOTE_RELEASE_CURRENT

REPOSITORY_URL = https://github.com/gaozw23/OnsetBounds

REMOTE_CREATED = YES

PUSH_PERFORMED = YES

VISIBILITY = PUBLIC

STATISTICAL_FREEZE_COMMIT = 6df6e4b7679a4dcfcd9ddf17b9e5a7f9dab2d097

STATISTICAL_FREEZE_TAG = v1.0.0-submission

RELEASE_URL = https://github.com/gaozw23/OnsetBounds/releases/tag/v1.0.0-submission

LICENSE_STATUS = MIT

PYTEST_FUNCTIONS = 17

LOGICAL_MANDATORY_CASES = 21

TEST_MAPPING_STATUS = PARTIAL_DIRECT_COVERAGE_DISCLOSED

MIMIC_PRIMARY_BOOTSTRAP = 5000/5000

EICU_PRIMARY_TWO_STAGE_BOOTSTRAP = 5000/5000

EICU_CONDITIONAL_BOOTSTRAP = 1000/1000

FINAL_BOOTSTRAP_LOCK_HASH = a9335f5f768459276a7e2834aefd2b78f7824cd8adab04ed4244a553a6afdedd

STATISTICAL_ANALYSIS_RERUN = NO

NUMERICAL_RESULTS_CHANGED = NO

The 17 pytest functions and 21 logical requirements are not interchangeable counts. `docs/TEST_CASE_MAPPING.csv` records direct, partial, and absent direct coverage without inventing tests. In particular, the simulated cohort-containment requirement is not a direct public unit test.

The statistical implementation used for the frozen numerical results is identified by the statistical freeze commit and tag above. Later documentation, provenance, and presentation commits do not alter that implementation or its reported numerical results.

Historical pre-remote manifest and audit files are retained under `docs/release/history/` with explicit provenance notices. Their historical state fields remain unchanged.
