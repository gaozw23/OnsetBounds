# Basic public-release security check

Scope follows SIMPLE PUBLIC RELEASE MODE. The final public tree was checked for forbidden `data`/`private` directories, restricted patient-level file extensions, patient-identifier columns in public CSV headers, credential assignments, and local Windows paths.

RESTRICTED_PATIENT_DATA_FOUND = NO

EMBEDDED_PATIENT_IDENTIFIERS_FOUND = NO

SECRETS_FOUND = NO

ABSOLUTE_LOCAL_PATHS_FOUND = NO

UNEXPLAINED_LARGE_FILES_FOUND = NO

- Forbidden `data` or `private` directories: 0.
- Parquet/Feather/Arrow/DuckDB/model/pickle/compressed-CSV artifacts: 0.
- Public CSV headers containing patient or raw hospital identifier fields: 0.
- Secret assignments detected: 0.
- Local Windows paths detected: 0.
- Files over 10 MB: 0.
- Archives: 0.

Schema names such as `subject_id` in source code are legitimate database interfaces and are not embedded identifier values.
