# Third-party dependencies and data-access notices

A source-file provenance review classified all tracked public source code as `ORIGINAL`. No source code was copied or adapted from MIT-LCP mimic-code, MIT-LCP eicu-code, another GitHub repository, Stack Overflow, or another third-party source. No third-party copyright or license header is present in the public source tree.

Database schema names, clinical item IDs, and configuration facts are interoperability metadata rather than copied source code. They do not require third-party source-code attribution. The following external packages and controlled databases retain their own licenses or access terms:

| Dependency or data ecosystem | Status | Public repository use |
|---|---|---|
| MIMIC-IV 3.1 | Credentialed restricted database; no data redistributed | Schema interfaces, item IDs, and documentation only |
| eICU Collaborative Research Database 2.0 | Credentialed restricted database; no data redistributed | Schema interfaces and documentation only |
| NumPy, pandas, SciPy, scikit-learn, XGBoost, Matplotlib, PyYAML, pytest | Normal package dependencies under their respective upstream licenses | Imported APIs only; no dependency source is vendored |

The repository's original source code is licensed under the MIT License in `LICENSE`.

LICENSE_STATUS = MIT
