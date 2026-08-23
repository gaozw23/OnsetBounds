# Data access

The final analyses used MIMIC-IV 3.1 and the eICU Collaborative Research Database 2.0. Both are distributed through PhysioNet under credentialing, training, and data-use requirements. This repository redistributes neither database.

Authorized users should download the databases directly from PhysioNet, keep them outside the repository, and expose their local roots through `MIMIC_IV_ROOT` or `EICU_ROOT` (or an external YAML configuration). The pipeline templates deliberately contain no default machine path. Do not commit raw data, source-row excerpts, identifiers, or patient-derived artifacts.
