"""Rebuild the five manuscript figures using only public aggregate CSV files."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODELS = ["ELASTIC_NET", "XGBOOST", "TREND"]
LABELS = {"ELASTIC_NET": "Elastic net", "XGBOOST": "XGBoost", "TREND": "Trend"}
COLORS = {"ELASTIC_NET": "#0072B2", "XGBOOST": "#D55E00", "TREND": "#009E73"}


def save(fig: plt.Figure, output: Path, name: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / name, bbox_inches="tight", metadata={"Creator": "OnsetBounds public aggregate-data builder", "CreationDate": None, "ModDate": None})
    plt.close(fig)


def figure1(output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(6.25, 4.85), gridspec_kw={"height_ratios": [1.15, 1]})
    ax = axes[0]
    ax.set_title("(A) Single-policy geometry", loc="left", fontweight="bold")
    ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis("off")
    ax.plot([3, 8], [3.3, 3.3], lw=6, color="#90CAF9", solid_capstyle="butt")
    ax.text(5.5, 3.55, r"valid window $W=[A+\ell,A+u]$", ha="center")
    rows = [(2.6, 4, 7, "certain valid", r"$v^-=v^+=1$"), (1.65, 1.8, 5, "ambiguous", r"$[v^-,v^+]=[0,1]$"), (.7, .5, 2.4, "impossible", r"$v^-=v^+=0$")]
    for y, left, right, label, result in rows:
        ax.plot([left, right], [y, y], lw=5, color="#F9A825", solid_capstyle="round")
        ax.text(.05, y, label, va="center")
        ax.text(8.35, y, result, va="center")
    ax = axes[1]
    ax.set_title("(B) Paired comparison with a shared latent onset", loc="left", fontweight="bold")
    ax.set_xlim(0, 10); ax.set_ylim(0, 3.3); ax.axis("off")
    ax.plot([1.1, 6.8], [2.5, 2.5], lw=6, color="#0072B2", label=r"$W_a$")
    ax.plot([3.1, 8.8], [1.85, 1.85], lw=6, color="#D55E00", label=r"$W_b$")
    ax.plot([1.8, 8.2], [.9, .9], lw=3, color="black", label=r"$\Gamma=[L,U]$")
    ax.text(1.2, 2.8, r"$W_a$"); ax.text(8.55, 2.15, r"$W_b$"); ax.text(4.75, .42, r"same $T\in\Gamma$")
    ax.text(8.95, 2.7, r"$[d_{ab,-},d_{ab,+}]=[-1,1]$", ha="right")
    save(fig, output, "figure1_geometry.pdf")


def figure2(data: Path, output: Path) -> None:
    primary = pd.read_csv(data / "primary_summary.csv")
    oracle = pd.read_csv(data / "oracle_summary.csv")
    labels = {"right_endpoint": "Right endpoint", "midpoint": "Midpoint", "uniform_mi": "Uniform MI", "working_weibull_aft_mi": "Weibull AFT MI"}
    colors = {"right_endpoint": "#7f3c8d", "midpoint": "#11a579", "uniform_mi": "#3969ac", "working_weibull_aft_mi": "#e68310"}
    markers = {"right_endpoint": "o", "midpoint": "s", "uniform_mi": "^", "working_weibull_aft_mi": "D"}
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0))
    ax = axes[0, 0]
    points = primary[(primary.target == "a") & primary.method.isin(colors)]
    for method in colors:
        grouped = points[points.method == method].groupby("nominal_width").agg(x=("realized_mean_width", "mean"), y=("bias", "mean")).reset_index()
        ax.plot(grouped.x, grouped.y, marker=markers[method], color=colors[method], lw=1.2, ms=4, label=labels[method])
    ax.axhline(0, color=".45", lw=.7, ls="--"); ax.set(title="A. Point-estimator bias", xlabel="Realized mean interval width (hours)", ylabel="Bias in policy-a valid sensitivity"); ax.legend(frameon=False, ncol=2)
    ax = axes[0, 1]
    bounds = primary[(primary.target == "a") & (primary.method == "OnsetBounds")].merge(oracle[["oracle_key", "a_lower", "a_upper"]], on="oracle_key")
    elo, ehi = bounds.a_lower + bounds.lower_endpoint_bias, bounds.a_upper + bounds.upper_endpoint_bias
    ax.scatter(bounds.a_lower, elo, s=12, facecolors="none", edgecolors="#3969ac", label="Lower endpoint")
    ax.scatter(bounds.a_upper, ehi, s=12, color="#e68310", alpha=.75, label="Upper endpoint")
    limits = [min(bounds.a_lower.min(), elo.min()), max(bounds.a_upper.max(), ehi.max())]
    ax.plot(limits, limits, color=".4", lw=.8, ls="--"); ax.set(xlim=limits, ylim=limits, title="B. Estimated versus oracle endpoints", xlabel="Oracle endpoint", ylabel="Mean estimated endpoint"); ax.legend(frameon=False)
    ax = axes[1, 0]
    coverage = primary[(primary.target == "a") & (primary.method == "OnsetBounds")]
    for width, subset in coverage.groupby("nominal_width"):
        ax.scatter(np.repeat(width, len(subset)) + np.linspace(-.16, .16, len(subset)), subset.complete_identified_set_coverage, s=9, color=".65", alpha=.55)
    means = coverage.groupby("nominal_width")["complete_identified_set_coverage"].mean()
    ax.plot(means.index, means.values, color="#3969ac", marker="o", lw=1.3, label="Mean across prespecified cells"); ax.axhline(.95, color="#c73e1d", lw=.8, ls="--", label="0.95")
    ax.set(ylim=(.88, 1.005), title="C. Outer-set coverage", xlabel="Nominal interval width (hours)", ylabel="Coverage"); ax.legend(frameon=False)
    ax = axes[1, 1]
    paired = primary[(primary.target == "paired") & (primary.method == "direct_paired_OnsetBounds")]
    for scenario, marker in zip(("A", "B", "C"), ("o", "s", "^")):
        subset = paired[paired.scenario == scenario]
        ax.scatter(subset.mean_marginal_subtraction_width, subset.mean_direct_paired_width, s=18, marker=marker, label=f"Scenario {scenario}", alpha=.8)
    maximum = float(paired.mean_marginal_subtraction_width.max())
    ax.plot([0, maximum], [0, maximum], color=".4", lw=.8, ls="--"); ax.set(xlim=(0, maximum*1.03), ylim=(0, maximum*1.03), title="D. Paired identified-set width", xlabel="Marginal-subtraction width", ylabel="Direct shared-onset width"); ax.legend(frameon=False)
    fig.tight_layout(pad=1, w_pad=1.5, h_pad=1.5); save(fig, output, "figure2_simulation_results.pdf")


def operating_panel(ax: plt.Axes, curve: pd.DataFrame, outer: dict[tuple[str, str], tuple[float, float]], database: str) -> None:
    for model in MODELS:
        rows = curve[curve.model == model].sort_values("test_non_event_first_alert_rate")
        x = 100 * rows.test_non_event_first_alert_rate.to_numpy(); lo = rows.valid_lower.to_numpy(); hi = rows.valid_upper.to_numpy(); mid = (lo + hi) / 2
        ax.plot(x, mid, color=COLORS[model], lw=.9, alpha=.8); ax.vlines(x, lo, hi, color=COLORS[model], lw=1.3); ax.scatter(x, mid, s=12, facecolor="white", edgecolor=COLORS[model], lw=.8)
        primary = rows[np.isclose(rows.target, .10)].iloc[0]; xp = 100 * primary.test_non_event_first_alert_rate; ol, ou = outer[(database, model)]
        ax.vlines(xp, ol, ou, color=COLORS[model], lw=.8, linestyles="dashed"); ax.scatter([xp], [(primary.valid_lower + primary.valid_upper)/2], s=24, color=COLORS[model], label=LABELS[model])
    ax.set(xlim=(0, None), ylim=(0, 1), xlabel="TEST non-event first-alert rate (%)", ylabel="Valid-sensitivity bound"); ax.grid(axis="y", color="#D9D9D9", lw=.5); ax.legend(frameon=False, loc="upper left")


def figures3_to_5(data: Path, output: Path) -> None:
    single = pd.read_csv(data / "FINAL_BOOTSTRAP_SINGLE_POLICY_CI_DATA.csv")
    outer = {(r.database, r.model): (r.outer_lower, r.outer_upper) for r in single.itertuples() if r.metric == "valid"}
    mimic_curve = pd.read_csv(data / "mimic_operating_curve_results.csv"); eicu_curve = pd.read_csv(data / "EICU_ALL_OPERATING_POINTS.csv")
    eicu_summary = pd.read_csv(data / "EICU_TEST_COHORT_SUMMARY.csv").iloc[0]; mimic_resolution = pd.read_csv(data / "mimic_interval_resolution_summary.csv").iloc[0]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), constrained_layout=True)
    operating_panel(axes[0], mimic_curve, outer, "MIMIC"); axes[0].set_title("A  MIMIC-IV TEST")
    operating_panel(axes[1], eicu_curve, outer, "eICU"); axes[1].set_title("B  eICU TEST replication"); axes[1].set_ylabel("")
    for idx, (label, values, color) in enumerate([("MIMIC-IV", [mimic_resolution["median"], mimic_resolution.p90], "#0072B2"), ("eICU", [eicu_summary.median_interval_hours, eicu_summary.p90_interval_hours], "#D55E00")]):
        axes[2].plot(values, [idx, idx], color=color, lw=3, solid_capstyle="round"); axes[2].scatter(values, [idx, idx], color=["white", color], edgecolor=color, s=28, lw=1)
        for value in values: axes[2].text(value, idx+.13, f"{value:.2f}", ha="center", va="bottom", fontsize=7)
    axes[2].set(yticks=[0,1], yticklabels=["MIMIC-IV","eICU"], ylim=(-.5,1.5), xlabel="Onset-bracket width (hours)", title="C  Temporal resolution: median to p90"); axes[2].grid(axis="x", color="#D9D9D9", lw=.5)
    save(fig, output, "figure3_real_data_performance.pdf")
    paired = pd.read_csv(data / "FINAL_BOOTSTRAP_PAIRED_CI_DATA.csv"); contrasts = ["ELASTIC_NET-XGBOOST","ELASTIC_NET-TREND","XGBOOST-TREND"]; labels=["Elastic net - XGBoost","Elastic net - Trend","XGBoost - Trend"]
    fig, axes = plt.subplots(1,2,figsize=(7.2,2.55),sharex=True,constrained_layout=True)
    for ax,database,title in zip(axes,["MIMIC","eICU"],["A  MIMIC-IV TEST","B  eICU TEST replication"]):
        rows={r.contrast:r for r in paired[paired.database==database].itertuples()}; y=np.arange(3)[::-1]
        for yi,c in zip(y,contrasts):
            r=rows[c]; ax.plot([r.marginal_lower,r.marginal_upper],[yi,yi],color="#9E9E9E",lw=5,solid_capstyle="round"); ax.plot([r.simultaneous_outer_lower,r.simultaneous_outer_upper],[yi,yi],color="#222",lw=1,ls="--"); ax.plot([r.direct_lower,r.direct_upper],[yi,yi],color="#0072B2",lw=3,solid_capstyle="round")
        ax.axvline(0,color="#444",lw=.7); ax.set(yticks=y,yticklabels=labels,xlabel="Difference in valid sensitivity",title=title); ax.grid(axis="x",color="#E0E0E0",lw=.5)
    axes[1].set_yticklabels([]); save(fig, output, "figure4_paired_comparisons.pdf")
    mimic_primary=pd.read_csv(data/"mimic_primary_r010_performance.csv").set_index("model"); eicu_primary=pd.read_csv(data/"EICU_PRIMARY_R010_RESULTS.csv").set_index("model"); retesting=pd.read_csv(data/"EICU_RETESTING_GAP_RESULTS.csv")
    fig,axes=plt.subplots(1,3,figsize=(7.2,2.5),constrained_layout=True); x=np.arange(3)
    axes[0].bar(x-.18,100*np.array([mimic_resolution.prop_gt6,mimic_resolution.prop_gt12,mimic_resolution.prop_gt24]),.36,label="MIMIC-IV",color="#0072B2"); axes[0].bar(x+.18,100*np.array([eicu_summary.prop_gt6h,eicu_summary.prop_gt12h,eicu_summary.prop_gt24h]),.36,label="eICU",color="#D55E00"); axes[0].set(xticks=x,xticklabels=[">6 h",">12 h",">24 h"],ylim=(0,100),ylabel="Finite brackets (%)",title="A  Onset-resolution tails"); axes[0].legend(frameon=False)
    axes[1].bar(x-.18,[mimic_primary.loc[m].identification_width for m in MODELS],.36,color="#0072B2"); axes[1].bar(x+.18,[eicu_primary.loc[m].identification_width for m in MODELS],.36,color="#D55E00"); axes[1].set(xticks=x,xticklabels=["EN","XGB","Trend"],ylim=(0,.7),ylabel="Valid-sensitivity ID width",title="B  Primary identification width")
    q1=retesting[retesting.group=="Q1"].set_index("model"); q4=retesting[retesting.group=="Q4"].set_index("model"); axes[2].bar(x-.18,[q1.loc[m].valid_interval_width for m in MODELS],.36,label="Q1: median 14.13 h",color="#56B4E9"); axes[2].bar(x+.18,[q4.loc[m].valid_interval_width for m in MODELS],.36,label="Q4: median 24.00 h",color="#CC79A7"); axes[2].set(xticks=x,xticklabels=["EN","XGB","Trend"],ylim=(0,.7),ylabel="Valid-sensitivity ID width",title="C  eICU retesting-gap quartiles"); axes[2].legend(frameon=False,loc="upper left")
    for ax in axes: ax.grid(axis="y",color="#E0E0E0",lw=.5)
    save(fig,output,"figure5_measurement_resolution.pdf")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--input",type=Path,required=True); parser.add_argument("--output",type=Path,required=True); args=parser.parse_args()
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"axes.titlesize":9.5,"axes.labelsize":8.5,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,"axes.linewidth":.7,"pdf.fonttype":42})
    figure1(args.output); figure2(args.input,args.output); figures3_to_5(args.input,args.output)
    print("FIGURE_REPRODUCTION_PASSED=5/5")


if __name__ == "__main__": main()
