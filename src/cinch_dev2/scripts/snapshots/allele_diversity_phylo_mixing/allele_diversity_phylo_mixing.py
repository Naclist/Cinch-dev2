#!/usr/bin/env python3
"""Plot per-locus allele diversity against same-allele phylogenetic mixing."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def canonical_genome(name: str) -> str:
    match = re.search(r"(GC[AF]_\d+)", str(name))
    if match:
        return match.group(1)
    return str(name).removesuffix(".fna").removesuffix(".fa").removesuffix(".fasta")


class Node:
    def __init__(self) -> None:
        self.name = ""
        self.length = 0.0
        self.children: list["Node"] = []


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")

    def parse_label(node: Node, position: int) -> int:
        start = position
        while position < len(text) and text[position] not in ",():":
            position += 1
        node.name = text[start:position].strip().strip("'\"")
        if position < len(text) and text[position] == ":":
            position += 1
            start = position
            while position < len(text) and text[position] not in ",()":
                position += 1
            try:
                node.length = float(text[start:position])
            except ValueError:
                node.length = 0.0
        return position

    def parse_node(position: int) -> tuple[Node, int]:
        node = Node()
        if text[position] == "(":
            position += 1
            while True:
                child, position = parse_node(position)
                node.children.append(child)
                if position >= len(text):
                    break
                if text[position] == ",":
                    position += 1
                    continue
                if text[position] == ")":
                    position += 1
                    break
        return node, parse_label(node, position)

    root, _ = parse_node(0)
    return root


def tree_distance_matrix(tree_path: Path, genomes: pd.Index) -> np.ndarray:
    root = parse_newick(tree_path.read_text())
    tips: dict[str, Node] = {}
    depths: dict[Node, float] = {}
    paths: dict[Node, list[Node]] = {}

    def walk(node: Node, depth: float, path: list[Node]) -> None:
        depths[node] = depth
        paths[node] = [*path, node]
        if not node.children and node.name:
            tips[canonical_genome(node.name)] = node
        for child in node.children:
            walk(child, depth + child.length, paths[node])

    walk(root, 0.0, [])
    nodes = []
    missing = []
    for genome in genomes:
        key = canonical_genome(genome)
        if key not in tips:
            missing.append(str(genome))
        else:
            nodes.append(tips[key])
    if missing:
        raise ValueError(
            f"Tree is missing {len(missing)} profile genomes; first missing: {missing[:5]}"
        )

    distance = np.zeros((len(nodes), len(nodes)), dtype=np.float32)
    ancestor_sets = [set(paths[node]) for node in nodes]
    for left in range(len(nodes) - 1):
        for right in range(left + 1, len(nodes)):
            lca = root
            for ancestor in paths[nodes[left]]:
                if ancestor in ancestor_sets[right]:
                    lca = ancestor
            value = depths[nodes[left]] + depths[nodes[right]] - 2.0 * depths[lca]
            distance[left, right] = distance[right, left] = value
    return distance


def locus_metrics(
    values: np.ndarray,
    tree_distance: np.ndarray,
    minimum_allele_count: int,
) -> dict[str, float]:
    called = values[values > 0]
    alleles, counts = np.unique(called, return_counts=True)
    frequencies = counts / counts.sum()
    simpson = 1.0 - float(np.sum(frequencies**2))

    distance_sum = 0.0
    pair_count = 0
    alleles_used = 0
    for allele, count in zip(alleles, counts):
        if count < minimum_allele_count:
            continue
        indices = np.flatnonzero(values == allele)
        pairs = len(indices) * (len(indices) - 1) // 2
        if pairs == 0:
            continue
        distance_sum += float(np.triu(tree_distance[np.ix_(indices, indices)], k=1).sum())
        pair_count += pairs
        alleles_used += 1

    return {
        "prevalence": float(len(called) / len(values)),
        "raw_n_alleles": int(len(alleles)),
        "simpson_diversity": simpson,
        "same_allele_mean_tree_distance": (
            distance_sum / pair_count if pair_count else np.nan
        ),
        "same_allele_pairs_used": int(pair_count),
        "same_allele_alleles_used": int(alleles_used),
    }


def plot_metrics(metrics: pd.DataFrame, output: Path) -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": "#FCFCFD",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#D8DCE5",
        },
    )
    fig, ax = plt.subplots(figsize=(10.5, 7.5))
    points = ax.scatter(
        metrics["same_allele_mean_tree_distance"],
        metrics["simpson_diversity"],
        c=metrics["prevalence"],
        cmap="viridis",
        s=18,
        alpha=0.68,
        linewidths=0,
        rasterized=True,
    )
    colorbar = fig.colorbar(points, ax=ax, pad=0.015)
    colorbar.set_label("Locus prevalence")
    ax.set_xlabel("Mean phylogenetic distance among genomes sharing the same allele")
    ax.set_ylabel("Allele diversity (Simpson)")
    ax.set_title(
        "Per-locus allele diversity vs phylogenetic mixing",
        loc="left",
        fontsize=15,
        fontweight="bold",
    )
    ax.grid(axis="both", color="#E7E9EF", linewidth=0.8)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_two_panels(metrics: pd.DataFrame, output: Path) -> None:
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": "#FCFCFD",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#D8DCE5",
        },
    )
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.8), sharex=True, sharey=True)
    panels = [
        ("cgMLST_like", "A. cgMLST-like loci (prevalence >=95%)"),
        ("wgMLST_accessory", "B. wgMLST accessory loci (1-95%)"),
    ]
    points = None
    for ax, (locus_set, title) in zip(axes, panels):
        block = metrics[metrics["locus_set"] == locus_set]
        points = ax.scatter(
            block["same_allele_mean_tree_distance"],
            block["simpson_diversity"],
            c=block["prevalence"],
            cmap="viridis",
            vmin=0,
            vmax=1,
            s=16,
            alpha=0.67,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{title}\n{len(block):,} loci plotted", loc="left", fontsize=12, fontweight="bold")
        ax.set_xlabel("Mean phylogenetic distance among genomes sharing the same allele")
        ax.grid(axis="both", color="#E7E9EF", linewidth=0.8)
        sns.despine(ax=ax)
    axes[0].set_ylabel("Allele diversity (Simpson)")
    colorbar = fig.colorbar(points, ax=axes, pad=0.015, fraction=0.025)
    colorbar.set_label("Locus prevalence")
    fig.suptitle(
        "Allele diversity vs phylogenetic mixing across cgMLST and wgMLST loci",
        x=0.08,
        y=0.98,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.86, bottom=0.12, left=0.08, right=0.93, wspace=0.08)
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate Simpson allele diversity and the pair-weighted mean tree "
            "distance among genomes sharing the same allele for every locus."
        )
    )
    parser.add_argument("--profile", type=Path, required=True, help="Tab-delimited allelic profile.")
    parser.add_argument("--tree", type=Path, required=True, help="Newick tree with branch lengths.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--core-min",
        type=float,
        default=0.95,
        help="Minimum called fraction required for a locus (default: 0.95).",
    )
    parser.add_argument(
        "--wg-min",
        type=float,
        default=0.01,
        help="Minimum called fraction for a wgMLST accessory locus (default: 0.01).",
    )
    parser.add_argument(
        "--min-allele-count",
        type=int,
        default=20,
        help="Minimum carriers required for an allele to contribute to the X axis (default: 20).",
    )
    parser.add_argument(
        "--include-wg-panel",
        action="store_true",
        help="Additionally draw the non-overlapping cgMLST/wgMLST two-panel figure.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(args.profile, sep="\t", index_col=0).fillna(0)
    profile = profile.apply(pd.to_numeric, errors="coerce").fillna(0)
    called_fraction = (profile > 0).mean(axis=0)
    selected = profile.loc[:, called_fraction >= args.wg_min]
    if selected.empty:
        raise ValueError("No loci pass --wg-min")

    distance = tree_distance_matrix(args.tree, selected.index)
    all_pair_mean = float(distance[np.triu_indices(len(distance), k=1)].mean())
    rows = []
    for gene in selected.columns:
        values = selected[gene].to_numpy(dtype=np.int64)
        rows.append(
            {
                "gene": gene,
                **locus_metrics(values, distance, args.min_allele_count),
            }
        )
    metrics = pd.DataFrame(rows).dropna(
        subset=["simpson_diversity", "same_allele_mean_tree_distance"]
    )
    metrics["locus_set"] = np.where(
        metrics["prevalence"] >= args.core_min,
        "cgMLST_like",
        "wgMLST_accessory",
    )
    metrics.to_csv(args.output_dir / "gene_allele_diversity_phylo_mixing.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "profile_genomes": len(profile),
                "profile_loci": profile.shape[1],
                "loci_passing_wg_min": selected.shape[1],
                "loci_plotted": len(metrics),
                "core_min": args.core_min,
                "wg_min": args.wg_min,
                "cgmlst_like_loci_plotted": int((metrics["locus_set"] == "cgMLST_like").sum()),
                "wgmlst_accessory_loci_plotted": int((metrics["locus_set"] == "wgMLST_accessory").sum()),
                "min_allele_count": args.min_allele_count,
                "all_genome_pair_mean_tree_distance": all_pair_mean,
            }
        ]
    ).to_csv(args.output_dir / "run_summary.tsv", sep="\t", index=False)
    plot_metrics(
        metrics[metrics["locus_set"] == "cgMLST_like"],
        args.output_dir / "allele_diversity_vs_same_allele_tree_distance",
    )
    if args.include_wg_panel:
        plot_two_panels(
            metrics,
            args.output_dir / "allele_diversity_vs_same_allele_tree_distance_cg_wg_panels",
        )


if __name__ == "__main__":
    main()
