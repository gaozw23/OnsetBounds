# OnsetBounds

Reference code for **Sharp Performance Bounds for Clinical Early-Warning Systems with Interval-Uncertain Event Onset**.

OnsetBounds evaluates an alert against every event time compatible with an observed onset bracket. It returns sharp bounds for valid, premature, late/missed, PPV, and bounded-utility performance. Direct paired bounds preserve the shared latent onset and can be tighter than subtracting separate marginal intervals.

## Repository scope

The core package, simulations, unit tests, synthetic example, aggregate table/figure sources, and figure builder are public. MIMIC-IV and eICU pipelines require separately authorized local database access; restricted source data and patient-level derivatives are not included.

## Quick start

Python 3.12 is the frozen reference version.

```bash
python -m venv .venv
python -m pip install -e ".[simulation,dev]"
python -m pytest -q
python synthetic/generate_synthetic_example.py --config synthetic/example_config.yaml --output synthetic/example_output
python figures/build_figures.py --input results/aggregate/figure_data --output reproduced_figures
```

The synthetic command creates new pseudo-random records and does not use clinical data. Generated output directories are ignored by Git.

## Authorized real-data reproduction

The frozen applications used MIMIC-IV 3.1 and the eICU Collaborative Research Database 2.0. Both require PhysioNet credentialing and applicable data-use agreements. Configure local roots through environment variables or an external YAML file; never put raw data or fitted patient-derived artifacts in this repository. See `docs/DATA_ACCESS.md`, `docs/MIMIC_PIPELINE.md`, and `docs/EICU_PIPELINE.md`.

Aggregate manuscript sources are under `results/aggregate`; all five figures can be rebuilt without private data. Elastic Net, XGBoost, and Trend are evaluation vehicles. Robust ranking remained unresolved. No MIMIC-to-eICU coefficient transport was reported because exact semantic compatibility was unsupported.

## Citation and license review

Use `CITATION.cff` for citation metadata. Repository licensing is intentionally marked `REVIEW_REQUIRED`; see `LICENSE_REVIEW_REQUIRED.md` and `docs/THIRD_PARTY_NOTICES.md` before public release.
