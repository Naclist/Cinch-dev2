#!/usr/bin/env python3
"""Run allele-diversity/phylogenetic-mixing plots for top species and complexes."""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))
from allele_diversity_phylo_mixing import locus_metrics, tree_distance_matrix


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")


def draw_panel(ax: plt.Axes, metrics: pd.DataFrame, title: str, sample_count: int) -> None:
    ax.scatter(
        metrics["same_allele_mean_tree_distance"],
        metrics["simpson_diversity"],
        c=metrics["prevalence"],
        cmap="viridis",
        vmin=0.95,
        vmax=1.0,
        s=15,
        alpha=0.68,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title(
        f"{title}\nn={sample_count:,}; loci={len(metrics):,}",
        loc="left",
        fontsize=10,
        fontweight="bold",
    )
    ax.set_xlabel("Same-allele mean tree distance")
    ax.set_ylabel("Allele diversity (Simpson)")
    ax.grid(color="#E7E9EF", linewidth=0.8)
    sns.despine(ax=ax)


def calculate_subset(
    profile: pd.DataFrame,
    distance: np.ndarray,
    indices: np.ndarray,
    core_min: float,
    minimum_allele_count: int,
) -> pd.DataFrame:
    subset = profile.iloc[indices]
    called_fraction = (subset > 0).mean(axis=0)
    selected = subset.loc[:, called_fraction >= core_min]
    subset_distance = distance[np.ix_(indices, indices)]
    rows = []
    for gene in selected.columns:
        rows.append(
            {
                "gene": gene,
                **locus_metrics(
                    selected[gene].to_numpy(dtype=np.int64),
                    subset_distance,
                    minimum_allele_count,
                ),
            }
        )
    return pd.DataFrame(rows).dropna(
        subset=["simpson_diversity", "same_allele_mean_tree_distance"]
    )


def save_individual(metrics: pd.DataFrame, title: str, sample_count: int, output: Path) -> None:
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(9.4, 7.0))
    draw_panel(ax, metrics, title, sample_count)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run top-three species and complex subsets.")
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--core-min", type=float, default=0.95)
    parser.add_argument("--min-allele-count", type=int, default=20)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, sep="\t").set_index("genome")
    top_species = metadata["species"].value_counts().head(args.top_n).index.tolist()
    valid_complex = metadata.loc[
        metadata["complex"].notna()
        & ~metadata["complex"].isin(["", "non_complex", "unknown", "Unknown", "NA"]),
        "complex",
    ]
    top_complex = valid_complex.value_counts().head(args.top_n).index.tolist()
    groups = [
        *[("species", value) for value in top_species],
        *[("complex", value) for value in top_complex],
    ]

    summary_rows = []
    for profile_path in args.profiles:
        profile_name = safe_name(profile_path.name.replace(".gene_CGAV.profile", "").replace(".profile_rm", "_rm"))
        profile = pd.read_csv(profile_path, sep="\t", index_col=0).fillna(0)
        profile = profile.apply(pd.to_numeric, errors="coerce").fillna(0)
        missing_metadata = profile.index.difference(metadata.index)
        if len(missing_metadata):
            raise ValueError(f"{profile_path.name}: {len(missing_metadata)} genomes missing metadata")
        distance = tree_distance_matrix(args.tree, profile.index)
        profile_dir = args.output_dir / profile_name
        profile_dir.mkdir(parents=True, exist_ok=True)
        panel_results = []

        for group_type, group_value in groups:
            mask = metadata.loc[profile.index, group_type].eq(group_value).to_numpy()
            indices = np.flatnonzero(mask)
            metrics = calculate_subset(
                profile,
                distance,
                indices,
                args.core_min,
                args.min_allele_count,
            )
            group_dir = profile_dir / group_type / safe_name(group_value)
            group_dir.mkdir(parents=True, exist_ok=True)
            metrics.to_csv(group_dir / "gene_allele_diversity_phylo_mixing.tsv", sep="\t", index=False)
            save_individual(
                metrics,
                f"{group_type}: {group_value}",
                len(indices),
                group_dir / "allele_diversity_vs_same_allele_tree_distance",
            )
            panel_results.append((group_type, group_value, len(indices), metrics))
            summary_rows.append(
                {
                    "profile": profile_name,
                    "group_type": group_type,
                    "group": group_value,
                    "n_genomes": len(indices),
                    "n_loci_plotted": len(metrics),
                    "core_min": args.core_min,
                    "min_allele_count": args.min_allele_count,
                }
            )

        sns.set_theme(style="whitegrid")
        fig, axes = plt.subplots(2, 3, figsize=(17, 10))
        for ax, (group_type, group_value, sample_count, metrics) in zip(axes.flat, panel_results):
            draw_panel(ax, metrics, f"{group_type}: {group_value}", sample_count)
        fig.suptitle(
            f"{profile_name}: top species and species complexes",
            x=0.06,
            y=0.99,
            ha="left",
            fontsize=16,
            fontweight="bold",
        )
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(profile_dir / "top3_species_top3_complex_overview.png", dpi=260, bbox_inches="tight")
        fig.savefig(profile_dir / "top3_species_top3_complex_overview.svg", bbox_inches="tight")
        plt.close(fig)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.output_dir / "top_group_run_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"group_type": group_type, "rank": rank, "group": value}
            for group_type, values in [("species", top_species), ("complex", top_complex)]
            for rank, value in enumerate(values, start=1)
        ]
    ).to_csv(args.output_dir / "selected_top_groups.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
