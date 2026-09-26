#!/usr/bin/env python3
"""Regenerate source-backed Cinch dev2 documentation figures."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "figures"
SITE = ROOT / "site" / "assets"
BLUE = "#2457A6"
GOLD = "#D4A72C"
ORANGE = "#D66A2C"
INK = "#202631"
MUTED = "#6B7280"
GRID = "#D9DEE7"
BG = "#F7F8FA"


def save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SITE.mkdir(parents=True, exist_ok=True)
    for directory in (OUT, SITE):
        fig.savefig(directory / f"{name}.png", dpi=220, bbox_inches="tight", facecolor="white",
                    metadata={"Software": "Cinch dev2 deterministic figure generator"})
        fig.savefig(directory / f"{name}.svg", bbox_inches="tight", facecolor="white",
                    metadata={"Date": None, "Creator": "Cinch dev2 deterministic figure generator"})
    plt.close(fig)


def epidis_validation() -> None:
    frame = pd.read_csv(ROOT / "docs" / "source-results" / "epidis_correlations_with_mi.tsv", sep="\t")
    labels = ["MI v2", "MI v2 bits\n(EpiDis²)", "old weighted MI", "old weighted NMI"]
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(11, 5.8))
    fig.subplots_adjust(top=0.78)
    width = 0.34
    ax.barh(y + width / 2, frame["pearson"], height=width, color=BLUE, label="Pearson")
    ax.barh(y - width / 2, frame["spearman"], height=width, color=GOLD, label="Spearman")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Correlation across 2,609,249 frozen locus pairs")
    fig.text(0.125, 0.94, "Pairwise EpiDis is a square-root rescaling of weighted MI v2",
             fontsize=15, fontweight="bold", color=INK, va="top")
    fig.text(0.125, 0.89, "EpiDis² and MI v2 in bits agree to a maximum absolute difference of 2.87×10⁻¹⁵",
             fontsize=10.5, color=MUTED, va="top")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.legend(frameon=False, loc="lower right")
    for series, offset in (("pearson", width / 2), ("spearman", -width / 2)):
        for yi, value in zip(y, frame[series]):
            ax.text(min(value + 0.012, 1.005), yi + offset, f"{value:.3f}", va="center", fontsize=9, color=INK)
    save(fig, "FIG02_EPIDIS_VALIDATION")


def diff_gwes_validation() -> None:
    summary = json.loads((ROOT / "docs" / "source-results" / "shc_diff_gwes_validation_summary.json").read_text())
    labels = ["Planted truths", "Predeclared negatives", "Natural-data hypotheses"]
    tested = np.array([summary["toy_truth_n"], summary["toy_negative_n"], summary["true_eligible"]])
    passed = np.array([summary["toy_truth_pass"], summary["toy_negative_pass"], summary["true_significant"]])
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    fig.subplots_adjust(top=0.77)
    y = np.arange(3)
    ax.barh(y, tested, color="#E7EAF0", edgecolor="#B7BEC9", label="Tested/eligible")
    ax.barh(y, passed, color=BLUE, label="Significant/passed")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Hypotheses")
    fig.text(0.125, 0.94, "Diff-GWES benchmark separates toy truth from natural-data evidence",
             fontsize=14.5, fontweight="bold", color=INK, va="top")
    fig.text(0.125, 0.89, "Gray = tested/eligible; blue = significant/passed. Toy recovery does not establish biological power",
             fontsize=10.5, color=MUTED, va="top")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    for yi, (n, p) in enumerate(zip(tested, passed)):
        ax.text(n + max(tested) * 0.015, yi, f"{p}/{n}", va="center", fontsize=10, color=INK)
    save(fig, "FIG03_DIFF_GWES_VALIDATION")


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": INK,
                         "svg.hashsalt": "cinch-dev2-v0.1.0"})
    epidis_validation()
    diff_gwes_validation()
    print(f"Wrote figures to {OUT} and {SITE}")


if __name__ == "__main__":
    main()
