# Reproducibility

## Public workflow

Use Python 3.12 and install the simulation and development dependencies:

```bash
python -m pip install -e ".[simulation,dev]"
```

Run the public test suite:

```bash
python -m pytest -q
```

Run the synthetic example:

```bash
python paper_repro/synthetic/generate_synthetic_example.py --config paper_repro/synthetic/example_config.yaml --output paper_repro/synthetic/example_output
```

The small public simulation runner can be invoked with:

```bash
python scripts/simulation/run_simulation.py --config scripts/simulation/configs/example.yaml --output simulation_example.csv
```

It illustrates the package interface; the completed manuscript simulation results are supplied as aggregate files under `paper_repro/aggregate/`.

Rebuild the manuscript figures from aggregate data:

```bash
python paper_repro/build_figures.py --input paper_repro/aggregate/figure_data --output paper_repro/reproduced_figures
```

## Clinical data workflow

The [MIMIC-IV](../scripts/mimic/README.md) and [eICU](../scripts/eicu/README.md) modules provide portable interfaces for authorized local data. Raw databases, patient-level derived records, fitted objects, and bootstrap draws are not included. See [DATA_ACCESS.md](DATA_ACCESS.md) for access requirements and [PROVENANCE.md](PROVENANCE.md) for the code version associated with the reported results.
