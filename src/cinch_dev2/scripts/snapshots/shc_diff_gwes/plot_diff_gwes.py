#!/usr/bin/env python3
"""Plot true-data and toy-benchmark SHC-conditioned Diff-GWES results."""

from pathlib import Path
import textwrap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"
OUT = WF / "output/figures"
INK, MUTED, GRID = "#20252C", "#68717D", "#E2E7EB"
BLUE, ORANGE, GREY = "#3976A8", "#D46A3A", "#B8C1C9"


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B8C0C8")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)


def header(fig, title, subtitle):
    fig.text(0.06, 0.975, title, ha="left", va="top", fontsize=17, weight="bold", color=INK)
    fig.text(0.06, 0.938, textwrap.fill(subtitle, 105), ha="left", va="top", fontsize=9, color=MUTED)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    true = pd.read_csv(WF / "output/true/validation_candidates.tsv", sep="\t")
    toy = pd.read_csv(WF / "output/toy/validation_candidates.tsv", sep="\t")
    summaries = pd.concat([pd.read_csv(WF / f"output/{mode}/summary.tsv", sep="\t").assign(mode=mode) for mode in ("true", "toy")], ignore_index=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.labelsize": 9, "axes.titlesize": 10, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), facecolor="white")
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.08, top=0.87, hspace=0.36, wspace=0.25)
    header(fig, "SHC-conditioned Diff-GWES: true data and injected positive control", "Discovery and validation samples are disjoint within SHC. Validation uses 999 SHC-preserving permutations; toy data contain 20 known strong background-dependent triples.")

    ax = axes[0, 0]
    ax.scatter(true.discovery_delta, true.validation_delta, s=15, alpha=0.55, color=BLUE, edgecolors="none")
    lim = max(abs(true.discovery_delta).max(), abs(true.validation_delta.dropna()).max(), 0.2)
    ax.plot([-lim, lim], [-lim, lim], color=INK, linestyle="--", linewidth=1)
    ax.axhline(0, color=GREY, linewidth=0.8); ax.axvline(0, color=GREY, linewidth=0.8)
    ax.set(xlabel="Discovery Delta", ylabel="Validation Delta", title="A  True data: effect replication")
    style(ax)

    ax = axes[0, 1]
    nontruth = toy[~toy.is_injected_truth]
    truth = toy[toy.is_injected_truth]
    ax.scatter(nontruth.discovery_delta, nontruth.validation_delta, s=12, alpha=0.28, color=GREY, edgecolors="none", label="Non-injected candidate")
    ax.scatter(truth.discovery_delta, truth.validation_delta, s=34, alpha=0.95, color=ORANGE, edgecolors=INK, linewidths=0.3, label="Injected truth")
    lim = max(abs(toy.discovery_delta).max(), abs(toy.validation_delta.dropna()).max(), 1.0)
    ax.plot([-lim, lim], [-lim, lim], color=INK, linestyle="--", linewidth=1)
    ax.axhline(0, color=GREY, linewidth=0.8); ax.axvline(0, color=GREY, linewidth=0.8)
    ax.set(xlabel="Discovery Delta", ylabel="Validation Delta", title="B  Toy data: 20 injected triples")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    style(ax)

    ax = axes[1, 0]
    bins = np.linspace(0, 1, 21)
    true_p = true.loc[true.validation_eligible, "perm_p"].dropna()
    toy_null_p = toy.loc[toy.validation_eligible & toy.is_negative_control, "perm_p"].dropna()
    ax.hist(true_p, bins=bins, histtype="step", linewidth=1.8, color=BLUE, label=f"True candidates (n={len(true_p)})")
    ax.hist(toy_null_p, bins=bins, histtype="step", linewidth=1.8, color=ORANGE, label=f"Predeclared negative controls (n={len(toy_null_p)})")
    ax.axhline(len(true_p) / 20, color=BLUE, linestyle=":", linewidth=1)
    ax.axhline(len(toy_null_p) / 20, color=ORANGE, linestyle=":", linewidth=1)
    ax.set(xlabel="Empirical permutation p", ylabel="Candidate count", title="C  Held-out permutation p distributions")
    ax.legend(frameon=False, fontsize=7)
    style(ax)

    ax = axes[1, 1]
    labels = ["Truth recovered", "Truth missed", "Negative-control hit", "True-data hit"]
    toy_summary = summaries[summaries["mode"] == "toy"].iloc[0]
    true_summary = summaries[summaries["mode"] == "true"].iloc[0]
    values = [int(toy_summary.truth_significant), int(toy_summary.n_injected_truth - toy_summary.truth_significant), int(toy_summary.negative_controls_significant), int(true_summary.n_significant)]
    colors = [ORANGE, GREY, GREY, BLUE]
    bars = ax.bar(labels, values, color=colors, width=0.62)
    ax.bar_label(bars, padding=3, fontsize=9, color=INK)
    ax.set(ylabel="Triples", title="D  Positive-control outcome")
    ax.tick_params(axis="x", rotation=18)
    ax.set_ylim(0, max(values) * 1.22 + 1)
    style(ax)
    fig.savefig(OUT / "shc_diff_gwes_true_toy_validation.png", dpi=220)
    fig.savefig(OUT / "shc_diff_gwes_true_toy_validation.svg")
    plt.close(fig)

    attrition = []
    for row in summaries.itertuples(index=False):
        for stage, value in [("Discovery tests", row.n_discovery_tests), ("Discovery eligible", row.n_discovery_eligible), ("Validation selected", row.n_validation_candidates), ("Validation eligible", row.n_validation_eligible), ("Significant", row.n_significant)]:
            attrition.append({"mode": row.mode, "stage": stage, "n": value})
    attrition = pd.DataFrame(attrition)
    attrition.to_csv(OUT / "pipeline_attrition.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(9.5, 4.8), facecolor="white")
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.17, top=0.78)
    header(fig, "SHC-conditioned Diff-GWES pipeline attrition", "Counts are shown on a log scale; significant toy triples are the 20 injected positive controls, while true Kp data yield no held-out q<0.05 triple.")
    stages = attrition.stage.unique()
    x = np.arange(len(stages)); width = 0.36
    for offset, (mode, color) in zip((-width/2, width/2), (("true", BLUE), ("toy", ORANGE))):
        vals = attrition[attrition["mode"] == mode].set_index("stage").loc[stages, "n"].to_numpy()
        ax.bar(x + offset, np.maximum(vals, 0.8), width, color=color, label=mode.capitalize())
        for xx, value in zip(x + offset, vals):
            ax.text(xx, max(value, 0.8) * 1.18, f"{int(value):,}", ha="center", va="bottom", fontsize=7, color=INK)
    ax.set_yscale("log"); ax.set_xticks(x, stages); ax.set_ylabel("Triples (log scale)")
    ax.legend(frameon=False); style(ax)
    fig.savefig(OUT / "shc_diff_gwes_pipeline_attrition.png", dpi=220)
    fig.savefig(OUT / "shc_diff_gwes_pipeline_attrition.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
