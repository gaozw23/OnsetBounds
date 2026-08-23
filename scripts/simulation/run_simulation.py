"""Small public simulation runner; not the frozen manuscript simulation run."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from onsetbounds.bounds import single_policy_bounds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    rng = np.random.default_rng(int(config["seed"]))
    rows = []
    for replication in range(int(config["replications"])):
        n = int(config["n"])
        event = rng.random(n) < float(config["event_prevalence"])
        onset = rng.uniform(12, 48, n)
        width = float(config["interval_width_hours"])
        lower = np.maximum(0, onset - rng.uniform(0, width, n))
        upper = lower + width
        alert = onset - rng.normal(12, 8, n)
        bounds = [single_policy_bounds(float(l), float(u), float(a)) for l, u, a, d in zip(lower, upper, alert, event) if d]
        rows.append({
            "replication": replication,
            "events": len(bounds),
            "valid_lower": float(np.mean([x.valid_lower for x in bounds])) if bounds else np.nan,
            "valid_upper": float(np.mean([x.valid_upper for x in bounds])) if bounds else np.nan,
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
