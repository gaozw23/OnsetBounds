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

Version `1.1.0` (`v1.1.0-jbi-revision`) adds the JBI revision method code below. The frozen submission analyses remain at tag `v1.0.0-submission` and commit `6df6e4b7679a4dcfcd9ddf17b9e5a7f9dab2d097`; their core statistical implementation is unchanged. See [docs/PROVENANCE.md](docs/PROVENANCE.md) for the submission bootstrap provenance.

### JBI revision additions

- [`scripts/revision/bootstrap/`](scripts/revision/bootstrap/): A1 raw-scale bootstrap validation, the three-policy contribution layout, independent-unit oracle moments and uncertainty classification, and the conditional-binomial finite-sample guard for iid patient binary targets. The [fixed synthetic configuration](configs/revision/bootstrap.json) retains the final scenarios, pilot policies, seeds and sampling settings. Geometry and confidence sets call the existing `onsetbounds` functions.
- [`scripts/revision/clinical_sensitivity.py`](scripts/revision/clinical_sensitivity.py): C1 shared-onset width decomposition, alert concordance and right-endpoint references; C2 early/late onset-location restriction and exact ranking change points; D3 closed lead-window sensitivity at 0–6, 6–12, 6–24 and 6–48 hours. The original primary-evaluable population, brackets and first alerts stay fixed.
- [`scripts/revision/S1_published_score_core/`](scripts/revision/S1_published_score_core/): Supplementary Software S1, the existing 11-file public method package for the local exploratory application of the published Ugwuowo et al. (2020) AKI score. It includes the formula, retrospective timeline, unit mappings, frozen policy configuration and previously supplied synthetic test sources. It does not establish original hospital deployment equivalence or revalidate the published model's performance.

Install revision dependencies with `python -m pip install -e ".[revision]"`. For A1, add `scripts/revision/bootstrap` to the Python module search path and import `validation.validate_replication`; it accepts the JSON configuration plus independently constructed `(18, 2)` population endpoint and MCSE arrays. `oracle.empty_state`, `add_units` and `estimates` accumulate independent patient units, or hospital-summed contributions for the hospital target. The published files are portable method components; internal production orchestration, checkpoints and result archives are excluded.

The clinical module accepts aligned, authorized local arrays in hours since each patient's monitoring origin, using `+inf` for no alert. `shared_onset`, `lead_windows` and `onset_location` take the unchanged primary-evaluable event subset with one column per policy; `right_endpoint_reference` takes the full cohort and its original event/evaluability masks. `alert_concordance` accepts a chosen frozen population. These functions return aggregates without loading clinical files. Follow the S1 README for its distinct nanosecond timeline interface and landmark policy semantics.

No patient-level clinical data, derived scores, alert times or trajectories are distributed. Complete reproduction of restricted clinical analyses still requires authorized data access and the corresponding local preparation. No tests, bootstrap runs, simulations, model training, threshold selection or clinical TEST analyses were rerun for this release; the included synthetic test files are source artifacts only.

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
