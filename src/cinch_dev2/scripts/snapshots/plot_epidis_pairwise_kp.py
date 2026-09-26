#!/usr/bin/env python3
"""Visualize Kp pairwise EpiDis, MI equivalence, and distance structure."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, LinearSegmentedColormap


INK = "#20252C"
MUTED = "#68717D"
GRID = "#E3E7EB"
BLUE = "#3976A8"
BLUE_LIGHT = "#DCEAF3"
ORANGE = "#D46A3A"
GOLD = "#B38A2E"
WHITE = "#FFFFFF"
CMAP = LinearSegmentedColormap.from_list("density", [WHITE, "#DCEAF3", "#78A9C6", "#22577A"])


def style_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#B8C0C8")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, linewidth=0.7, axis="y", zorder=0)
    ax.set_axisbelow(True)


def figure_header(fig, title: str, subtitle: str) -> None:
    fig.text(0.06, 0.975, title, ha="left", va="top", fontsize=17, weight="bold", color=INK)
    fig.text(0.06, 0.94, textwrap.fill(subtitle, 150), ha="left", va="top", fontsize=9, color=MUTED)


def distance_bins(frame: pd.DataFrame, distance: str, bins: int = 80) -> pd.DataFrame:
    positive = frame[np.isfinite(frame[distance]) & (frame[distance] >= 0)].copy()
    transformed = np.log10(positive[distance].to_numpy() + 1.0)
    edges = np.linspace(transformed.min(), transformed.max(), bins + 1)
    group = np.clip(np.digitize(transformed, edges) - 1, 0, bins - 1)
    positive["bin"] = group
    result = positive.groupby("bin", observed=True).agg(
        n_pairs=("epidis", "size"),
        distance_median=(distance, "median"),
        epidis_median=("epidis", "median"),
        epidis_q995=("epidis", lambda x: x.quantile(0.995)),
        outlier_fraction=("epidis_iqr_outlier", "mean"),
    ).reset_index()
    result["distance_field"] = distance
    return result


def plot_distance(ax, frame: pd.DataFrame, binned: pd.DataFrame, distance: str, label: str, threshold: float) -> None:
    valid = np.isfinite(frame[distance]) & np.isfinite(frame.epidis) & (frame[distance] >= 0)
    x = np.log10(frame.loc[valid, distance].to_numpy() + 1.0)
    ax.hexbin(x, frame.loc[valid, "epidis"], gridsize=(85, 55), mincnt=1, cmap=CMAP, norm=LogNorm(), linewidths=0, rasterized=True)
    bx = np.log10(binned.distance_median.to_numpy() + 1.0)
    ax.plot(bx, binned.epidis_median, color=BLUE, linewidth=1.8, label="Median")
    ax.plot(bx, binned.epidis_q995, color=ORANGE, linewidth=1.8, label="99.5th percentile")
    ax.axhline(threshold, color=INK, linestyle="--", linewidth=1.1, label=f"IQR threshold {threshold:.3f}")
    ticks = np.arange(0, int(np.ceil(x.max())) + 1)
    ax.set_xticks(ticks, [f"$10^{t}$" if t else "0" for t in ticks])
    ax.set_xlabel(f"{label} + 1 (log scale)", color=INK)
    ax.set_ylabel("EpiDis", color=INK)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    style_axis(ax)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("Leca/results/gene_gene_wgmlst/epidis_pairwise_kp/kp_pairwise_epidis_all_pairs.tsv.gz"))
    parser.add_argument("--threshold", type=Path, default=Path("Leca/results/gene_gene_wgmlst/epidis_pairwise_kp/epidis_iqr_threshold.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/results/gene_gene_wgmlst/epidis_pairwise_kp/figures"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    columns = ["gene1", "gene2", "epidis", "epidis_squared", "epidis_iqr_outlier", "mi_v2", "old_weighted_mi", "old_weighted_nmi", "mean_order_dist", "mean_bp_dist"]
    frame = pd.read_csv(args.input, sep="\t", compression="gzip", usecols=columns)
    threshold = float(pd.read_csv(args.threshold, sep="\t").threshold.iloc[0])
    order_bins = distance_bins(frame, "mean_order_dist")
    bp_bins = distance_bins(frame, "mean_bp_dist")
    pd.concat([order_bins, bp_bins], ignore_index=True).to_csv(args.outdir / "epidis_distance_binned_trends.tsv", sep="\t", index=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "bold", "axes.titlesize": 10, "axes.labelsize": 9, "savefig.facecolor": WHITE})
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor=WHITE)
    fig.subplots_adjust(left=0.06, right=0.98, bottom=0.07, top=0.88, hspace=0.34, wspace=0.23)
    figure_header(fig, "Kp pairwise EpiDis overview", "2,609,249 binary locus-presence pairs from 1,000 genomes; thesis Hamming-neighbour sample weights at tau=0.10; dashed line is the thesis default Q3 + 3 x IQR threshold.")

    ax = axes[0, 0]
    ax.hist(frame.epidis, bins=120, color=BLUE, edgecolor=WHITE, linewidth=0.2)
    ax.axvline(threshold, color=ORANGE, linestyle="--", linewidth=1.5, label=f"Threshold {threshold:.3f}")
    ax.set_yscale("log")
    ax.set_xlabel("EpiDis", color=INK)
    ax.set_ylabel("Pair count (log scale)", color=INK)
    ax.set_title("A  EpiDis distribution", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=8)
    style_axis(ax)

    ax = axes[0, 1]
    hb = ax.hexbin(frame.mi_v2, frame.epidis, gridsize=90, mincnt=1, cmap=CMAP, norm=LogNorm(), linewidths=0, rasterized=True)
    mi_grid = np.linspace(0, frame.mi_v2.max(), 300)
    ax.plot(mi_grid, np.sqrt(mi_grid / np.log(2.0)), color=ORANGE, linewidth=1.5, label=r"$\sqrt{MI_{nats}/\ln 2}$")
    ax.set_xlabel("Weighted MI v2 (nats)", color=INK)
    ax.set_ylabel("EpiDis", color=INK)
    ax.set_title("B  EpiDis versus weighted MI v2", loc="left", color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    style_axis(ax)

    plot_distance(axes[1, 0], frame, order_bins, "mean_order_dist", "Mean order distance", threshold)
    axes[1, 0].set_title("C  EpiDis versus gene-order distance", loc="left", color=INK)
    plot_distance(axes[1, 1], frame, bp_bins, "mean_bp_dist", "Mean BP distance", threshold)
    axes[1, 1].set_title("D  EpiDis versus physical distance", loc="left", color=INK)
    fig.savefig(args.outdir / "kp_epidis_overview.png", dpi=220)
    fig.savefig(args.outdir / "kp_epidis_overview.svg")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), facecolor=WHITE)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.13, top=0.80, wspace=0.24)
    figure_header(fig, "EpiDis compared with the previous weighted MI workflow", "The old workflow uses wg-content Jaccard weights and a 0.5 pseudocount; density is shown on a logarithmic count scale.")
    for ax, xcol, xlabel, panel_title, letter in [
        (axes[0], "old_weighted_mi", "Previous weighted MI (nats)", "EpiDis versus previous weighted MI", "A"),
        (axes[1], "old_weighted_nmi", "Previous weighted NMI", "EpiDis versus previous weighted NMI", "B"),
    ]:
        ax.hexbin(frame[xcol], frame.epidis, gridsize=90, mincnt=1, cmap=CMAP, norm=LogNorm(), linewidths=0, rasterized=True)
        ax.axhline(threshold, color=ORANGE, linestyle="--", linewidth=1.2)
        pearson = frame.epidis.corr(frame[xcol])
        spearman = frame.epidis.rank(method="average").corr(frame[xcol].rank(method="average"))
        ax.text(0.03, 0.96, f"Pearson = {pearson:.3f}\nSpearman = {spearman:.3f}", transform=ax.transAxes, va="top", color=INK, fontsize=8, bbox={"facecolor": WHITE, "edgecolor": "none", "alpha": 0.88})
        ax.set_xlabel(xlabel, color=INK)
        ax.set_ylabel("EpiDis", color=INK)
        ax.set_ylim(0, 1.02)
        ax.set_title(f"{letter}  {panel_title}", loc="left", color=INK)
        style_axis(ax)
    fig.savefig(args.outdir / "kp_epidis_vs_previous_mi.png", dpi=220)
    fig.savefig(args.outdir / "kp_epidis_vs_previous_mi.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
