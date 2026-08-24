# Data access

The clinical analyses used MIMIC-IV 3.1 and the eICU Collaborative Research Database 2.0. Users must obtain both databases through PhysioNet and follow the applicable credentialing, training, and data-use requirements.

This repository does not redistribute raw MIMIC-IV data, raw eICU data, or patient-level derived records. It contains source code and aggregate manuscript data only.

Keep local database copies outside the repository. The pipeline templates accept database roots through `MIMIC_IV_ROOT` or `EICU_ROOT`, or through an external YAML configuration. See the [MIMIC-IV](../scripts/mimic/README.md) and [eICU](../scripts/eicu/README.md) pipeline notes.
