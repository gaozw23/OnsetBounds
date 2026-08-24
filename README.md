# OnsetBounds

Code and reproducibility materials accompanying the manuscript
*Sharp performance bounds for clinical early-warning systems with interval-uncertain event onset: a partial-identification framework*
by Zhuowei Gao and Weipeng Sun,
School of Mathematics, Jilin University.

## Overview

Clinical event onset may be observed only within an interval \([L,U]\), rather than at an exact time. OnsetBounds computes sharp single-policy performance bounds and direct paired comparison bounds while preserving this interval uncertainty.

The framework separates identification uncertainty due to interval-uncertain onset from sampling uncertainty in the estimated bounds.

## Repository scope

- Core OnsetBounds implementation.
- Simulation code.
- MIMIC-IV and eICU authorized-data pipelines.
- Synthetic example.
- Aggregate manuscript figure and table data.
- Implementation tests.

Raw MIMIC-IV/eICU data and patient-level derivatives are not distributed in this repository.

## Manuscript code version

The code version frozen for the submission manuscript is:

`6df6e4b7679a4dcfcd9ddf17b9e5a7f9dab2d097`

Tag:

`v1.0.0-submission`

That commit and tag are preserved permanently. Subsequent repository commits may contain documentation, provenance, and presentation updates only; they do not change the statistical implementation used for the reported results.

## Repository structure

```text
.
├── configs/             Public configuration templates.
├── docs/                Data-access, methods, and reproducibility documentation.
├── paper_repro/         Synthetic example and public aggregate figure/table materials.
├── scripts/             Simulation and authorized real-data pipeline entry points.
├── src/onsetbounds/     Core Python package.
├── tests/               Public implementation tests.
├── pyproject.toml       Package metadata and dependency groups.
└── requirements.txt     Frozen reference environment.
```

## Data

The authorized real-data applications use MIMIC-IV 3.1 and the eICU Collaborative Research Database 2.0. Access requires PhysioNet credentialing and acceptance of the applicable data-use agreements. Raw data are not hosted on GitHub. See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md) for access and local configuration guidance.

## Environment and installation

```bash
git clone https://github.com/gaozw23/OnsetBounds.git
cd OnsetBounds
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[simulation,dev]"
```

## Reproduction

### 1. Lightweight verification

```bash
python -m pytest -q
```

### 2. Synthetic example

```bash
python paper_repro/synthetic/generate_synthetic_example.py --config paper_repro/synthetic/example_config.yaml --output paper_repro/synthetic/example_output
```

The command generates synthetic records only and does not use clinical data.

### 3. Simulations

Simulation entry points and the example configuration are under [`scripts/simulation/`](scripts/simulation/). Production simulations are computationally intensive; review the configuration before running them.

### 4. Authorized real-data analysis

The [`scripts/mimic/`](scripts/mimic/) and [`scripts/eicu/`](scripts/eicu/) directories document the authorized-data pipelines. Reproduction requires local access to the corresponding PhysioNet databases and cannot be performed from this repository alone.

### 5. Tables and figures

```bash
python paper_repro/build_figures.py --input paper_repro/aggregate/figure_data --output paper_repro/reproduced_figures
```

Manuscript figures can be rebuilt from the public aggregate data under [`paper_repro/aggregate/`](paper_repro/aggregate/); no patient-level data are required for this step.

## Tests

The public tests cover single-policy bounds, paired bounds, bootstrap and inference utilities, weighted average precision, and the synthetic pipeline.

Table S12 of the accompanying manuscript lists 21 logical implementation requirements. The public suite contains 17 pytest functions, and one function may cover more than one logical requirement. [`docs/TEST_CASE_MAPPING.csv`](docs/TEST_CASE_MAPPING.csv) distinguishes direct, partial, indirect, and not-directly-covered requirements without claiming 21 standalone tests.

## Citation

Please cite *Sharp performance bounds for clinical early-warning systems with interval-uncertain event onset: a partial-identification framework* by Zhuowei Gao and Weipeng Sun, School of Mathematics, Jilin University. Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## License

This project is released under the MIT License. See [`LICENSE`](LICENSE).
