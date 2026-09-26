#!/usr/bin/env python3
"""Fit the exponential allelic-turnover baseline and identify positive MPD outliers."""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from allele_diversity_phylo_mixing import tree_distance_matrix


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")


def theoretical_baseline(distance: np.ndarray, points: int = 1800) -> pd.DataFrame:
    upper = distance[np.triu_indices(len(distance), k=1)].astype(float)
    positive = upper[upper > 0]
    scale = float(np.median(positive))
    scaled_lambda = np.logspace(-4, 4, points)
    lambdas = scaled_lambda / scale
    h_values = np.empty(points)
    m_values = np.empty(points)
    chunk = 100
    for start in range(0, points, chunk):
        stop = min(start + chunk, points)
        weights = np.exp(-lambdas[start:stop, None] * upper[None, :])
        denominators = weights.mean(axis=1)
        h_values[start:stop] = denominators
        m_values[start:stop] = (weights * upper[None, :]).mean(axis=1) / denominators
    sdi = 1.0 - h_values
    delta_m = np.diff(m_values)
    delta_sdi = np.diff(sdi)
    m_scale = max(np.ptp(m_values), np.finfo(float).eps)
    sdi_scale = max(np.ptp(sdi), np.finfo(float).eps)
    segment = np.sqrt((delta_m / m_scale) ** 2 + (delta_sdi / sdi_scale) ** 2)
    arc = np.r_[0.0, np.cumsum(segment)]
    arc /= arc[-1]
    return pd.DataFrame(
        {
            "lambda": lambdas,
            "scaled_lambda": scaled_lambda,
            "model_same_probability": h_values,
            "model_sdi": sdi,
            "model_same_allele_mpd": m_values,
            "baseline_arc": arc,
        }
    )


def observed_locus(values: np.ndarray, distance: np.ndarray) -> dict[str, float]:
    called_indices = np.flatnonzero(values > 0)
    called = values[called_indices]
    alleles, counts = np.unique(called, return_counts=True)
    all_pairs = len(called) * (len(called) - 1) // 2
    same_pairs = int(np.sum(counts * (counts - 1) // 2))
    if all_pairs == 0 or same_pairs == 0:
        return {}
    same_distance_sum = 0.0
    for allele, count in zip(alleles, counts):
        if count < 2:
            continue
        indices = np.flatnonzero(values == allele)
        same_distance_sum += float(np.triu(distance[np.ix_(indices, indices)], k=1).sum())
    same_probability = same_pairs / all_pairs
    frequencies = counts / counts.sum()
    return {
        "prevalence": float(len(called) / len(values)),
        "raw_n_alleles": int(len(alleles)),
        "pairwise_same_probability": same_probability,
        "pairwise_sdi": 1.0 - same_probability,
        "simpson_with_replacement": 1.0 - float(np.sum(frequencies**2)),
        "same_allele_mpd": same_distance_sum / same_pairs,
        "same_allele_pairs": same_pairs,
    }


def project_to_baseline(metrics: pd.DataFrame, baseline: pd.DataFrame, distance_mean: float) -> pd.DataFrame:
    ordered = baseline.sort_values("model_sdi")
    sdi = metrics["pairwise_sdi"].to_numpy()
    metrics = metrics.copy()
    metrics["baseline_lambda"] = np.interp(sdi, ordered["model_sdi"], ordered["lambda"])
    metrics["baseline_scaled_lambda"] = np.interp(
        sdi, ordered["model_sdi"], ordered["scaled_lambda"]
    )
    metrics["baseline_arc"] = np.interp(sdi, ordered["model_sdi"], ordered["baseline_arc"])
    metrics["model_expected_mpd"] = np.interp(
        sdi, ordered["model_sdi"], ordered["model_same_allele_mpd"]
    )
    metrics["mpd_residual"] = metrics["same_allele_mpd"] - metrics["model_expected_mpd"]
    metrics["normalized_mpd_residual"] = metrics["mpd_residual"] / distance_mean
    metrics["below_model_sdi_range"] = sdi < ordered["model_sdi"].min()
    metrics["above_model_sdi_range"] = sdi > ordered["model_sdi"].max()
    metrics["outside_model_sdi_range"] = (
        metrics["below_model_sdi_range"] | metrics["above_model_sdi_range"]
    )

    valid = metrics.loc[~metrics["outside_model_sdi_range"], "normalized_mpd_residual"]
    center = float(valid.median())
    mad = float(np.median(np.abs(valid - center)))
    robust_sigma = 1.4826 * mad
    threshold = center + 3.0 * robust_sigma
    metrics["residual_center"] = center
    metrics["residual_robust_sigma"] = robust_sigma
    metrics["positive_outlier_threshold"] = threshold
    metrics["positive_mpd_residual_outlier"] = (
        ~metrics["outside_model_sdi_range"]
        & (metrics["normalized_mpd_residual"] > threshold)
    )
    metrics["phylogeny_external_candidate"] = (
        metrics["above_model_sdi_range"] | metrics["positive_mpd_residual_outlier"]
    )
    return metrics


def plot_result(
    metrics: pd.DataFrame,
    baseline: pd.DataFrame,
    title: str,
    output: Path,
) -> None:
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))
    inlier = ~metrics["phylogeny_external_candidate"]
    axes[0].scatter(
        metrics.loc[inlier, "same_allele_mpd"],
        metrics.loc[inlier, "pairwise_sdi"],
        s=13,
        alpha=0.50,
        color="#4F7CAC",
        linewidths=0,
        rasterized=True,
        label="Background-compatible",
    )
    axes[0].scatter(
        metrics.loc[~inlier, "same_allele_mpd"],
        metrics.loc[~inlier, "pairwise_sdi"],
        s=20,
        alpha=0.80,
        color="#D66B45",
        linewidths=0,
        rasterized=True,
        label="Phylogeny-external candidate",
    )
    axes[0].plot(
        baseline["model_same_allele_mpd"],
        baseline["model_sdi"],
        color="#20242C",
        linewidth=2.2,
        label="Exponential turnover baseline",
    )
    axes[0].set_xlabel("MPD among genomes sharing the same allele")
    axes[0].set_ylabel("Pairwise SDI = 1 - P(same allele)")
    axes[0].set_title("A. Observed loci and theoretical baseline", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].scatter(
        metrics.loc[inlier, "baseline_arc"],
        metrics.loc[inlier, "normalized_mpd_residual"],
        s=13,
        alpha=0.50,
        color="#4F7CAC",
        linewidths=0,
        rasterized=True,
    )
    axes[1].scatter(
        metrics.loc[~inlier, "baseline_arc"],
        metrics.loc[~inlier, "normalized_mpd_residual"],
        s=20,
        alpha=0.80,
        color="#D66B45",
        linewidths=0,
        rasterized=True,
    )
    threshold = float(metrics["positive_outlier_threshold"].iloc[0])
    center = float(metrics["residual_center"].iloc[0])
    axes[1].axhline(center, color="#20242C", linewidth=1.2)
    axes[1].axhline(
        threshold,
        color="#D66B45",
        linestyle=":",
        linewidth=1.5,
        label=f"Median + 3 MAD-sigma = {threshold:.3f}",
    )
    axes[1].set_xlabel("Position along exponential-turnover baseline (normalized arc length)")
    axes[1].set_ylabel("(Observed MPD - model MPD) / mean pairwise distance")
    axes[1].set_title("B. Baseline-coordinate transformation", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(color="#E7E9EF", linewidth=0.8)
        sns.despine(ax=ax)
    fig.suptitle(title, x=0.06, y=0.99, ha="left", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_summary(summary: pd.DataFrame, output: Path) -> None:
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8))
    sns.barplot(
        data=summary,
        x="species",
        y="phylogeny_external_candidate_fraction",
        hue="profile",
        ax=axes[0],
    )
    axes[0].set_title("A. Fraction of core loci above the turnover baseline", loc="left", fontweight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Candidate fraction")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].legend(frameon=False, fontsize=8)

    sns.barplot(
        data=summary,
        x="species",
        y="residual_center",
        hue="profile",
        ax=axes[1],
    )
    axes[1].axhline(0, color="#20242C", linewidth=1)
    axes[1].set_title("B. Median normalized residual (model offset)", loc="left", fontweight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Median (observed MPD - model MPD) / mean D")
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        sns.despine(ax=ax)
    fig.suptitle(
        "Species-level exponential-turnover baseline assessment",
        x=0.06,
        y=0.99,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Model a one-parameter vertical-turnover baseline.")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-species", type=int, default=3)
    parser.add_argument("--core-min", type=float, default=0.95)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, sep="\t").set_index("genome")
    species = metadata["species"].value_counts().head(args.top_species).index.tolist()
    summary_rows = []
    for profile_path in args.profiles:
        profile_name = safe_name(
            profile_path.name.replace(".gene_CGAV.profile", "").replace(".profile_rm", "_rm")
        )
        profile = pd.read_csv(profile_path, sep="\t", index_col=0).fillna(0)
        profile = profile.apply(pd.to_numeric, errors="coerce").fillna(0)
        full_distance = tree_distance_matrix(args.tree, profile.index)

        for species_name in species:
            mask = metadata.loc[profile.index, "species"].eq(species_name).to_numpy()
            indices = np.flatnonzero(mask)
            subset = profile.iloc[indices]
            distance = full_distance[np.ix_(indices, indices)]
            selected = subset.loc[:, (subset > 0).mean(axis=0) >= args.core_min]
            rows = []
            for gene in selected.columns:
                result = observed_locus(selected[gene].to_numpy(dtype=np.int64), distance)
                if result:
                    rows.append({"gene": gene, **result})
            metrics = pd.DataFrame(rows)
            baseline = theoretical_baseline(distance)
            pairwise_mean = float(distance[np.triu_indices(len(distance), k=1)].mean())
            metrics = project_to_baseline(metrics, baseline, pairwise_mean)

            result_dir = args.output_dir / profile_name / safe_name(species_name)
            result_dir.mkdir(parents=True, exist_ok=True)
            metrics.to_csv(result_dir / "locus_turnover_baseline_metrics.tsv", sep="\t", index=False)
            baseline.to_csv(result_dir / "theoretical_exponential_baseline.tsv", sep="\t", index=False)
            plot_result(
                metrics,
                baseline,
                f"{profile_name} | {species_name} | n={len(indices):,}",
                result_dir / "sdi_mpd_turnover_baseline_and_outliers",
            )
            valid = ~metrics["outside_model_sdi_range"]
            summary_rows.append(
                {
                    "profile": profile_name,
                    "species": species_name,
                    "n_genomes": len(indices),
                    "n_core_loci": len(metrics),
                    "n_loci_within_model_range": int(valid.sum()),
                    "n_above_model_sdi_range": int(metrics["above_model_sdi_range"].sum()),
                    "n_positive_mpd_residual_outliers": int(
                        metrics["positive_mpd_residual_outlier"].sum()
                    ),
                    "n_phylogeny_external_candidates": int(
                        metrics["phylogeny_external_candidate"].sum()
                    ),
                    "phylogeny_external_candidate_fraction": float(
                        metrics["phylogeny_external_candidate"].mean()
                    ),
                    "positive_residual_outlier_fraction_model_range": float(
                        metrics.loc[valid, "positive_mpd_residual_outlier"].mean()
                    ),
                    "residual_center": float(metrics["residual_center"].iloc[0]),
                    "residual_robust_sigma": float(metrics["residual_robust_sigma"].iloc[0]),
                    "positive_outlier_threshold": float(
                        metrics["positive_outlier_threshold"].iloc[0]
                    ),
                    "mean_pairwise_tree_distance": pairwise_mean,
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        args.output_dir / "species_turnover_outlier_summary.tsv",
        sep="\t",
        index=False,
    )
    plot_summary(summary, args.output_dir / "species_turnover_outlier_summary")


if __name__ == "__main__":
    main()
