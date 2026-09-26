#!/usr/bin/env python3
"""Color turnover-baseline plots by an external per-locus HGT frequency."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr


def plot_species(metrics: pd.DataFrame, baseline: pd.DataFrame, title: str, output: Path) -> None:
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))
    common = {
        "c": metrics["HGT_Frequency"],
        "cmap": "plasma",
        "vmin": 0.5,
        "vmax": 1.0,
        "s": 16,
        "alpha": 0.72,
        "linewidths": 0,
        "rasterized": True,
    }
    points = axes[0].scatter(
        metrics["same_allele_mpd"],
        metrics["pairwise_sdi"],
        **common,
    )
    axes[0].plot(
        baseline["model_same_allele_mpd"],
        baseline["model_sdi"],
        color="#20242C",
        linewidth=2.2,
        label="Exponential turnover baseline",
    )
    outliers = metrics[metrics["phylogeny_external_candidate"]]
    axes[0].scatter(
        outliers["same_allele_mpd"],
        outliers["pairwise_sdi"],
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.8,
        s=34,
        label=f"Phylogeny-external candidate ({len(outliers):,})",
    )
    axes[0].set_xlabel("MPD among genomes sharing the same allele")
    axes[0].set_ylabel("Pairwise SDI = 1 - P(same allele)")
    axes[0].set_title("A. HGT-colored loci and turnover baseline", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].scatter(
        metrics["baseline_arc"],
        metrics["normalized_mpd_residual"],
        **common,
    )
    axes[1].scatter(
        outliers["baseline_arc"],
        outliers["normalized_mpd_residual"],
        facecolors="none",
        edgecolors="#111111",
        linewidths=0.8,
        s=34,
    )
    center = float(metrics["residual_center"].iloc[0])
    threshold = float(metrics["positive_outlier_threshold"].iloc[0])
    axes[1].axhline(center, color="#20242C", linewidth=1.2)
    axes[1].axhline(
        threshold,
        color="#20242C",
        linestyle=":",
        linewidth=1.5,
        label=f"Positive residual threshold = {threshold:.3f}",
    )
    axes[1].set_xlabel("Position along exponential-turnover baseline (normalized arc length)")
    axes[1].set_ylabel("(Observed MPD - model MPD) / mean pairwise distance")
    axes[1].set_title("B. HGT-colored baseline-coordinate transformation", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8, loc="upper left")
    for ax in axes:
        ax.grid(color="#E7E9EF", linewidth=0.8)
        sns.despine(ax=ax)
    colorbar = fig.colorbar(points, ax=axes, pad=0.015, fraction=0.025)
    colorbar.set_label("HGT frequency")
    fig.suptitle(title, x=0.06, y=0.99, ha="left", fontsize=15, fontweight="bold")
    fig.subplots_adjust(top=0.88, bottom=0.12, left=0.06, right=0.91, wspace=0.18)
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Color rm turnover plots by HGT frequency.")
    parser.add_argument("--hgt", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()

    hgt = pd.read_csv(args.hgt, sep="\t")
    if not {"Locus", "HGT_Frequency"}.issubset(hgt.columns):
        raise ValueError("HGT table must contain Locus and HGT_Frequency columns")
    species_dirs = sorted(path for path in args.result_dir.iterdir() if path.is_dir())
    summary = []
    for species_dir in species_dirs:
        metrics_path = species_dir / "locus_turnover_baseline_metrics.tsv"
        baseline_path = species_dir / "theoretical_exponential_baseline.tsv"
        if not metrics_path.exists() or not baseline_path.exists():
            continue
        metrics = pd.read_csv(metrics_path, sep="\t").merge(
            hgt,
            left_on="gene",
            right_on="Locus",
            how="left",
            validate="one_to_one",
        )
        if metrics["HGT_Frequency"].isna().any():
            raise ValueError(
                f"{species_dir.name}: {metrics['HGT_Frequency'].isna().sum()} loci lack HGT frequency"
            )
        baseline = pd.read_csv(baseline_path, sep="\t")
        metrics.to_csv(species_dir / "locus_turnover_baseline_metrics_with_hgt.tsv", sep="\t", index=False)
        plot_species(
            metrics,
            baseline,
            f"usearch_rm | {species_dir.name} | HGT-colored",
            species_dir / "sdi_mpd_turnover_baseline_hgt_colored",
        )
        candidates = metrics["phylogeny_external_candidate"]
        correlation = spearmanr(
            metrics["HGT_Frequency"],
            metrics["normalized_mpd_residual"],
        )
        if candidates.any():
            rank_test_p = mannwhitneyu(
                metrics.loc[candidates, "HGT_Frequency"],
                metrics.loc[~candidates, "HGT_Frequency"],
                alternative="two-sided",
            ).pvalue
        else:
            rank_test_p = float("nan")
        summary.append(
            {
                "species": species_dir.name,
                "n_loci": len(metrics),
                "hgt_min": metrics["HGT_Frequency"].min(),
                "hgt_median": metrics["HGT_Frequency"].median(),
                "hgt_max": metrics["HGT_Frequency"].max(),
                "candidate_hgt_median": metrics.loc[
                    metrics["phylogeny_external_candidate"], "HGT_Frequency"
                ].median(),
                "background_hgt_median": metrics.loc[
                    ~metrics["phylogeny_external_candidate"], "HGT_Frequency"
                ].median(),
                "spearman_hgt_vs_normalized_residual": correlation.statistic,
                "spearman_p": correlation.pvalue,
                "candidate_vs_background_mannwhitney_p": rank_test_p,
            }
        )
    pd.DataFrame(summary).to_csv(args.result_dir / "hgt_coloring_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
