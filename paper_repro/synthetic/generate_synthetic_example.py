"""Generate a wholly synthetic interval-onset example."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from onsetbounds.bootstrap import outer_confidence_set, patient_bootstrap
from onsetbounds.bounds import single_policy_bounds
from onsetbounds.paired import paired_valid_bounds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    rng = np.random.default_rng(int(cfg["seed"]))
    n = int(cfg["patients"])
    event = rng.random(n) < float(cfg["event_prevalence"])
    latent = rng.uniform(18, 60, n)
    gap_before = rng.uniform(1, 14, n)
    gap_after = rng.uniform(1, 14, n)
    lower, upper = latent - gap_before, latent + gap_after
    policy_names = [p["name"] for p in cfg["policies"]]
    alerts = {p["name"]: np.where(rng.random(n) < 0.88, latent - float(p["alert_shift_hours"]) + rng.normal(0, 3, n), np.inf) for p in cfg["policies"]}
    patient = []
    summary = []
    for name in policy_names:
        bounds = [single_policy_bounds(float(lower[i]), float(upper[i]), float(alerts[name][i])) for i in np.flatnonzero(event)]
        lo = np.array([b.valid_lower for b in bounds], dtype=float)
        hi = np.array([b.valid_upper for b in bounds], dtype=float)
        records = list(zip(lo, hi))
        boot = patient_bootstrap(records, lambda x: np.array([np.mean([r[0] for r in x]), np.mean([r[1] for r in x])]), int(cfg["bootstrap_replicates"]), int(cfg["seed"]) + len(summary))
        outer = outer_confidence_set(float(lo.mean()), float(hi.mean()), boot[:, 0], boot[:, 1])
        summary.append({"object": name, "type": "single_policy", **outer})
        patient.extend({"synthetic_index": int(i), "policy": name, "lower": lower[i], "upper": upper[i], "alert": alerts[name][i]} for i in np.flatnonzero(event))
    a, b = policy_names[:2]
    paired = [paired_valid_bounds(float(lower[i]), float(upper[i]), float(alerts[a][i]), float(alerts[b][i])) for i in np.flatnonzero(event)]
    summary.append({"object": f"{a}-{b}", "type": "direct_paired", "estimate_lower": np.mean([x.lower for x in paired]), "estimate_upper": np.mean([x.upper for x in paired]), "outer_lower": np.nan, "outer_upper": np.nan, "critical": np.nan})
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(summary)
    result.to_csv(output / "synthetic_bounds.csv", index=False)
    fig, ax = plt.subplots(figsize=(6, 3.4))
    singles = result[result.type == "single_policy"].reset_index(drop=True)
    for i, row in singles.iterrows():
        ax.plot([row.estimate_lower, row.estimate_upper], [i, i], lw=6, solid_capstyle="round")
        ax.plot([row.outer_lower, row.outer_upper], [i, i], color="black", lw=1, linestyle="--")
    ax.set_yticks(range(len(singles)), singles.object)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Valid-sensitivity identified interval")
    ax.set_title("Wholly synthetic OnsetBounds example")
    fig.tight_layout()
    fig.savefig(output / "synthetic_bounds.pdf", metadata={"Creator": "OnsetBounds synthetic example"})
    plt.close(fig)
    (output / "README.txt").write_text("Generated from independent pseudo-random synthetic data; no clinical records were used.\n", encoding="utf-8")
    print(f"SYNTHETIC_EXAMPLE_PASSED patients={n} events={int(event.sum())} policies={len(policy_names)}")


if __name__ == "__main__":
    main()
