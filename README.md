# OnsetBounds

OnsetBounds is the research software accompanying *Sharp performance bounds for clinical early-warning systems with interval-uncertain event onset: a partial-identification framework* by Zhuowei Gao and Weipeng Sun, School of Mathematics, Jilin University.

## Overview

Clinical event onset may be observed only within an interval \([L,U]\), rather than at an exact time. OnsetBounds computes sharp single-policy performance bounds and direct paired comparison bounds while preserving this interval uncertainty.

The framework separates identification uncertainty due to interval-uncertain onset from sampling uncertainty in the estimated bounds.

## Installation

```bash
git clone https://github.com/gaozw23/OnsetBounds.git
cd OnsetBounds
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[simulation,dev]"
```

## Reproducing the paper

### Synthetic example

```bash
python paper_repro/synthetic/generate_synthetic_example.py --config paper_repro/synthetic/example_config.yaml --output paper_repro/synthetic/example_output
```

The command generates synthetic records only and does not use clinical data.

### Simulation

The public simulation runner and example configuration are under [`scripts/simulation/`](scripts/simulation/). The full manuscript design and aggregate outputs are described in [the reproducibility guide](docs/REPRODUCIBILITY.md).

### Figures and tables

```bash
python paper_repro/build_figures.py --input paper_repro/aggregate/figure_data --output paper_repro/reproduced_figures
```

The figure builder uses the public aggregate sources under [`paper_repro/aggregate/`](paper_repro/aggregate/).

### MIMIC-IV and eICU

Authorized-data interfaces and configuration templates are documented with the [MIMIC-IV](scripts/mimic/README.md) and [eICU](scripts/eicu/README.md) code. Reproduction of the clinical analyses requires local database access.

## Data access

Raw MIMIC-IV and eICU data and patient-level derived records are not distributed. See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md).

## Code version

The reported numerical results correspond to tag `v1.0.0-submission` and commit `6df6e4b7679a4dcfcd9ddf17b9e5a7f9dab2d097`. Later documentation changes do not alter that statistical implementation. See [docs/PROVENANCE.md](docs/PROVENANCE.md) for the bootstrap provenance.

## Tests

The public test suite contains 17 pytest functions covering the implementation checks summarized in Supplementary Table S12. The correspondence is listed in [`docs/TEST_CASE_MAPPING.csv`](docs/TEST_CASE_MAPPING.csv). Cohort containment is verified in the completed simulation outputs rather than by a standalone unit test.

```bash
python -m pytest -q
```

## Citation

Please cite *Sharp performance bounds for clinical early-warning systems with interval-uncertain event onset: a partial-identification framework* by Zhuowei Gao and Weipeng Sun, School of Mathematics, Jilin University. Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).

External Python packages retain their respective upstream licenses; no third-party source code is vendored in this repository.
