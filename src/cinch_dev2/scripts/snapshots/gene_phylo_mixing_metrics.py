#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import gaussian_filter1d
from sklearn.cluster import KMeans


def canonical_genome(x: str) -> str:
    match = re.search(r"(GC[AF]_\d+)", str(x))
    if match:
        return match.group(1)
    return str(x).removesuffix(".fna").removesuffix(".fa").removesuffix(".fasta")


class Node:
    def __init__(self) -> None:
        self.name = ""
        self.length = 0.0
        self.children: list["Node"] = []
        self.parent: "Node | None" = None


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")

    def parse_name_length(node: Node, pos: int) -> int:
        start = pos
        while pos < len(text) and text[pos] not in ",():":
            pos += 1
        node.name = text[start:pos].strip()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            try:
                node.length = float(text[start:pos])
            except ValueError:
                node.length = 0.0
        return pos

    def parse_node(pos: int) -> tuple[Node, int]:
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                child, pos = parse_node(pos)
                child.parent = node
                node.children.append(child)
                if pos >= len(text):
                    break
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        pos = parse_name_length(node, pos)
        return node, pos

    root, _ = parse_node(0)
    return root


def tree_distance_matrix(tree_path: Path, genome_order: pd.Index) -> np.ndarray:
    root = parse_newick(tree_path.read_text())
    tips: dict[str, Node] = {}
    depths: dict[Node, float] = {}
    ancestors: dict[Node, list[Node]] = {}

    def walk(node: Node, depth: float, path: list[Node]) -> None:
        depths[node] = depth
        ancestors[node] = path + [node]
        if not node.children and node.name:
            tips[canonical_genome(node.name)] = node
        for child in node.children:
            walk(child, depth + child.length, ancestors[node])

    walk(root, 0.0, [])
    nodes = []
    missing = []
    for genome in genome_order:
        key = canonical_genome(genome)
        if key in tips:
            nodes.append(tips[key])
        else:
            missing.append(str(genome))
    if missing:
        raise ValueError(f"Tree is missing {len(missing)} profile genomes, first missing: {missing[:5]}")

    n = len(nodes)
    dist = np.zeros((n, n), dtype=np.float32)
    ancestor_sets = [set(ancestors[node]) for node in nodes]
    for i in range(n - 1):
        ai = ancestors[nodes[i]]
        for j in range(i + 1, n):
            common = ancestor_sets[j]
            lca = root
            for a in ai:
                if a in common:
                    lca = a
            d = depths[nodes[i]] + depths[nodes[j]] - 2.0 * depths[lca]
            dist[i, j] = dist[j, i] = d
    return dist


def mean_upper_triangle(dist: np.ndarray) -> float:
    idx = np.triu_indices(dist.shape[0], k=1)
    return float(dist[idx].mean())


def mean_same_allele_distance(values: np.ndarray, dist: np.ndarray, min_allele_count: int) -> tuple[float, int, int]:
    values = values.astype(np.int64, copy=False)
    alleles, counts = np.unique(values[values > 0], return_counts=True)
    total_sum = 0.0
    total_pairs = 0
    used_alleles = 0
    for allele, count in zip(alleles, counts):
        if count < min_allele_count:
            continue
        idx = np.flatnonzero(values == allele)
        sub = dist[np.ix_(idx, idx)]
        pair_count = len(idx) * (len(idx) - 1) // 2
        if pair_count <= 0:
            continue
        total_sum += float(np.triu(sub, k=1).sum())
        total_pairs += pair_count
        used_alleles += 1
    if total_pairs == 0:
        return np.nan, 0, 0
    return total_sum / total_pairs, total_pairs, used_alleles


def allele_diversity(values: np.ndarray) -> dict[str, float]:
    pos = values[values > 0]
    if len(pos) == 0:
        return {
            "prevalence": 0.0,
            "raw_n_alleles": 0,
            "simpson_diversity": np.nan,
            "shannon_entropy": np.nan,
            "effective_alleles_simpson": np.nan,
            "effective_alleles_shannon": np.nan,
        }
    _, counts = np.unique(pos, return_counts=True)
    p = counts / counts.sum()
    simpson = 1.0 - float(np.sum(p * p))
    shannon = -float(np.sum(p * np.log(p)))
    return {
        "prevalence": float(len(pos) / len(values)),
        "raw_n_alleles": int(len(counts)),
        "simpson_diversity": simpson,
        "shannon_entropy": shannon,
        "effective_alleles_simpson": float(1.0 / np.sum(p * p)),
        "effective_alleles_shannon": float(np.exp(shannon)),
    }


def add_vertical_kmeans_layers(metrics: pd.DataFrame, n_layers: int, n_baseline_bins: int) -> pd.DataFrame:
    metrics = metrics.copy()
    x_col = "same_allele_mean_tree_distance"
    y_col = "simpson_diversity"
    work = metrics[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    n_bins = min(n_baseline_bins, max(8, work[x_col].nunique()))
    work["_bin"] = pd.qcut(work[x_col].rank(method="first"), q=n_bins, duplicates="drop")
    baseline = (
        work.groupby("_bin", observed=True)
        .agg(x_mid=(x_col, "median"), y_med=(y_col, "median"))
        .sort_values("x_mid")
        .reset_index(drop=True)
    )
    sigma = max(2.0, len(baseline) / 28.0)
    baseline["y_smooth"] = gaussian_filter1d(baseline["y_med"].to_numpy(dtype=float), sigma=sigma, mode="nearest")
    metrics["trend_baseline"] = np.interp(metrics[x_col], baseline["x_mid"], baseline["y_smooth"])
    metrics["vertical_residual"] = metrics[y_col] - metrics["trend_baseline"]

    km = KMeans(n_clusters=n_layers, random_state=42, n_init=50)
    labels = km.fit_predict(metrics[["vertical_residual"]])
    centers = pd.Series(km.cluster_centers_.ravel(), index=range(n_layers)).sort_values()
    name_by_label = {int(label): f"layer_{rank + 1}_{name}" for rank, (label, name) in enumerate(zip(centers.index, ["low", "middle", "high"]))}
    metrics["vertical_layer"] = [name_by_label[int(label)] for label in labels]
    metrics["vertical_layer_rank"] = metrics["vertical_layer"].str.extract(r"layer_(\d+)")[0].astype(int)
    metrics.attrs["baseline_curve"] = baseline
    metrics.attrs["layer_centers"] = centers.rename("vertical_residual_center").reset_index().rename(columns={"index": "kmeans_label"})
    return metrics


def _least_dense_boundary(y: np.ndarray, low_center: float, high_center: float) -> float:
    if not np.isfinite(low_center) or not np.isfinite(high_center) or high_center <= low_center:
        return float((low_center + high_center) / 2.0)
    candidates = np.linspace(low_center, high_center, 151)
    iqr = np.subtract(*np.quantile(y, [0.75, 0.25]))
    bandwidth = max(0.008, float(iqr) / 18.0)
    density = np.exp(-0.5 * ((y[:, None] - candidates[None, :]) / bandwidth) ** 2).sum(axis=0)
    return float(candidates[int(np.argmin(density))])


def add_local_valley_layers(metrics: pd.DataFrame, n_bins: int, smooth_sigma: float) -> pd.DataFrame:
    metrics = metrics.copy()
    x_col = "same_allele_mean_tree_distance"
    y_col = "simpson_diversity"
    work = metrics[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    n_bins = min(n_bins, max(8, work[x_col].nunique()))
    work["_bin"] = pd.qcut(work[x_col].rank(method="first"), q=n_bins, duplicates="drop")

    rows = []
    for _, block in work.groupby("_bin", observed=True):
        if len(block) < 18 or block[y_col].nunique() < 3:
            continue
        y = block[y_col].to_numpy(dtype=float)
        km = KMeans(n_clusters=3, random_state=42, n_init=30)
        km.fit(y.reshape(-1, 1))
        centers = np.sort(km.cluster_centers_.ravel())
        lower_boundary = _least_dense_boundary(y, centers[0], centers[1])
        upper_boundary = _least_dense_boundary(y, centers[1], centers[2])
        rows.append(
            {
                "x_mid": float(block[x_col].median()),
                "lower_boundary_raw": lower_boundary,
                "upper_boundary_raw": upper_boundary,
                "n": int(len(block)),
                "center_low": float(centers[0]),
                "center_middle": float(centers[1]),
                "center_high": float(centers[2]),
            }
        )

    boundaries = pd.DataFrame(rows).sort_values("x_mid").reset_index(drop=True)
    if boundaries.empty:
        raise ValueError("Could not estimate local valley layer boundaries")
    sigma = smooth_sigma if smooth_sigma > 0 else max(1.5, len(boundaries) / 22.0)
    boundaries["lower_boundary"] = gaussian_filter1d(boundaries["lower_boundary_raw"].to_numpy(dtype=float), sigma=sigma, mode="nearest")
    boundaries["upper_boundary"] = gaussian_filter1d(boundaries["upper_boundary_raw"].to_numpy(dtype=float), sigma=sigma, mode="nearest")
    swapped = boundaries["lower_boundary"] >= boundaries["upper_boundary"]
    if swapped.any():
        mid = (boundaries.loc[swapped, "lower_boundary"] + boundaries.loc[swapped, "upper_boundary"]) / 2.0
        boundaries.loc[swapped, "lower_boundary"] = mid - 0.01
        boundaries.loc[swapped, "upper_boundary"] = mid + 0.01

    lower = np.interp(metrics[x_col], boundaries["x_mid"], boundaries["lower_boundary"])
    upper = np.interp(metrics[x_col], boundaries["x_mid"], boundaries["upper_boundary"])
    metrics["valley_lower_boundary"] = lower
    metrics["valley_upper_boundary"] = upper
    metrics["vertical_layer"] = "layer_2_middle"
    metrics.loc[metrics[y_col] <= lower, "vertical_layer"] = "layer_1_low"
    metrics.loc[metrics[y_col] > upper, "vertical_layer"] = "layer_3_high"
    metrics["vertical_layer_rank"] = metrics["vertical_layer"].str.extract(r"layer_(\d+)")[0].astype(int)
    metrics["layer_method"] = "local_valley_boundary"
    metrics.attrs["valley_boundaries"] = boundaries
    return metrics


def plot_metrics(metrics: pd.DataFrame, out_base: Path, lineage_threshold: float, mixed_threshold: float) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9.2, 7.0))
    palette = {
        "lineage_strong": "#c23b22",
        "mixed_high": "#1f77b4",
        "intermediate": "#555555",
    }
    sns.scatterplot(
        data=metrics,
        x="same_allele_mean_tree_distance",
        y="simpson_diversity",
        hue="phylo_mixing_class",
        style="phylo_mixing_class",
        palette=palette,
        markers={"intermediate": "o", "mixed_high": "^", "lineage_strong": "X"},
        s=26,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.15,
        ax=ax,
    )
    ax.axvline(lineage_threshold, color="#c23b22", ls="--", lw=1.3)
    ax.axvline(mixed_threshold, color="#1f77b4", ls="--", lw=1.3)
    ax.set_xlabel("Mean patristic distance among genome pairs\nsharing the same allele")
    ax.set_ylabel("Allele diversity (Simpson index)")
    fig.suptitle("Per-gene allele diversity vs phylogenetic mixing", y=0.985, fontsize=19)
    fig.text(
        0.5,
        0.94,
        f"{len(metrics):,} loci; dashed lines mark empirical 10th and 90th percentiles",
        ha="center",
        va="top",
        fontsize=9.5,
        color="#555555",
    )
    ax.legend(
        title="Empirical class",
        frameon=False,
        markerscale=1.2,
        loc="upper right",
        fontsize=10,
        title_fontsize=10,
    )
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_vertical_layers(metrics: pd.DataFrame, out_base: Path, genome_pair_mean_distance: float) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9.2, 7.0))
    palette = {
        "layer_1_low": "#2a9d8f",
        "layer_2_middle": "#555555",
        "layer_3_high": "#e76f51",
    }
    sns.scatterplot(
        data=metrics,
        x="same_allele_mean_tree_distance",
        y="simpson_diversity",
        hue="vertical_layer",
        palette=palette,
        s=26,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.15,
        ax=ax,
    )

    boundaries = metrics.attrs.get("valley_boundaries")
    if boundaries is not None:
        ax.plot(boundaries["x_mid"], boundaries["lower_boundary"], color="#202020", lw=2.5, alpha=0.94)
        ax.plot(boundaries["x_mid"], boundaries["upper_boundary"], color="#202020", lw=2.5, alpha=0.94)
    baseline = metrics.attrs.get("baseline_curve")
    centers = metrics.attrs.get("layer_centers")
    if boundaries is None and baseline is not None and centers is not None:
        centers = centers.sort_values("vertical_residual_center")["vertical_residual_center"].to_numpy()
        boundaries = (centers[:-1] + centers[1:]) / 2.0
        for boundary in boundaries:
            ax.plot(
                baseline["x_mid"],
                baseline["y_smooth"] + boundary,
                color="#16a6df",
                lw=2.0,
                alpha=0.95,
            )
    ax.axvline(genome_pair_mean_distance, color="#202020", ls=":", lw=1.5, label="all genome-pair mean distance")
    ax.set_xlabel("mean phylogenetic distance among genomes sharing the same allele")
    ax.set_ylabel("allele diversity (Simpson)")
    ax.set_title("Per-gene vertical-layer clustering after phylogenetic trend correction")
    ax.legend(title="", frameon=False, markerscale=1.2)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_same_allele_distance_stacked_hist(metrics: pd.DataFrame, out_base: Path, bins: int, genome_pair_mean_distance: float) -> pd.DataFrame:
    x_col = "same_allele_mean_tree_distance"
    layers = ["layer_1_low", "layer_2_middle", "layer_3_high"]
    counts, edges = [], None
    for layer in layers:
        vals = metrics.loc[metrics["vertical_layer"] == layer, x_col].to_numpy(dtype=float)
        hist, edges = np.histogram(vals, bins=bins, range=(float(metrics[x_col].min()), float(metrics[x_col].max())))
        counts.append(hist)
    count_arr = np.vstack(counts)
    centers = (edges[:-1] + edges[1:]) / 2.0
    hist_table = pd.DataFrame({"bin_left": edges[:-1], "bin_right": edges[1:], "bin_mid": centers})
    for layer, hist in zip(layers, count_arr):
        hist_table[layer] = hist

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(11, 5.8))
    bottom = np.zeros(len(centers))
    colors = ["#2a9d8f", "#555555", "#e76f51"]
    width = np.diff(edges)
    for layer, hist, color in zip(layers, count_arr, colors):
        ax.bar(centers, hist, width=width, bottom=bottom, color=color, alpha=0.82, linewidth=0, label=layer)
        bottom += hist
    ax.axvline(genome_pair_mean_distance, color="#202020", ls=":", lw=1.6, label="all genome-pair mean distance")
    ax.set_xlabel("same-allele mean phylogenetic distance")
    ax.set_ylabel("number of genes")
    ax.set_title(f"Stacked histogram of same-allele mean distance ({bins} bins)")
    ax.legend(frameon=False, fontsize=9)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)
    return hist_table


def boundary_crossing_summary(metrics: pd.DataFrame, margins: list[float]) -> pd.DataFrame:
    rows = []
    for margin in margins:
        near_lower = (metrics["simpson_diversity"] - metrics["valley_lower_boundary"]).abs() <= margin
        near_upper = (metrics["simpson_diversity"] - metrics["valley_upper_boundary"]).abs() <= margin
        rows.append(
            {
                "margin": margin,
                "n_near_lower_boundary": int(near_lower.sum()),
                "n_near_upper_boundary": int(near_upper.sum()),
                "n_near_either_boundary": int((near_lower | near_upper).sum()),
                "fraction_near_either_boundary": float((near_lower | near_upper).mean()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Per-gene phylogenetic mixing metric from cgMLST alleles.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument("--tree", type=Path, default=Path("Leca/source/tree/Kp.labelled.nwk"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/results/gene_gene_cgmlst/phylo_mixing"))
    parser.add_argument("--core-min", type=float, default=0.95)
    parser.add_argument("--min-allele-count", type=int, default=20)
    parser.add_argument("--lineage-quantile", type=float, default=0.10)
    parser.add_argument("--mixed-quantile", type=float, default=0.90)
    parser.add_argument("--vertical-layers", type=int, default=3)
    parser.add_argument("--baseline-bins", type=int, default=160)
    parser.add_argument("--layer-bins", type=int, default=72)
    parser.add_argument("--layer-smooth-sigma", type=float, default=0.0)
    parser.add_argument("--hist-bins", type=int, default=1000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    profile = pd.read_csv(args.profile, sep="\t", index_col=0).fillna(0)
    profile = profile.loc[:, (profile > 0).mean(axis=0) >= args.core_min].copy()
    dist = tree_distance_matrix(args.tree, profile.index)
    genome_pair_mean_distance = mean_upper_triangle(dist)

    rows = []
    for gene in profile.columns:
        values = pd.to_numeric(profile[gene], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
        row = {"gene": gene, **allele_diversity(values)}
        mean_dist, n_pairs, n_alleles_used = mean_same_allele_distance(values, dist, args.min_allele_count)
        row.update(
            {
                "same_allele_mean_tree_distance": mean_dist,
                "same_allele_pairs_used": n_pairs,
                "same_allele_alleles_used": n_alleles_used,
            }
        )
        rows.append(row)

    metrics = pd.DataFrame(rows).dropna(subset=["same_allele_mean_tree_distance", "simpson_diversity"])
    lineage_threshold = float(metrics["same_allele_mean_tree_distance"].quantile(args.lineage_quantile))
    mixed_threshold = float(metrics["same_allele_mean_tree_distance"].quantile(args.mixed_quantile))
    metrics["phylo_mixing_class"] = "intermediate"
    metrics.loc[metrics["same_allele_mean_tree_distance"] <= lineage_threshold, "phylo_mixing_class"] = "lineage_strong"
    metrics.loc[metrics["same_allele_mean_tree_distance"] >= mixed_threshold, "phylo_mixing_class"] = "mixed_high"
    metrics = add_local_valley_layers(metrics, args.layer_bins, args.layer_smooth_sigma)
    metrics.to_csv(args.outdir / "gene_phylo_mixing_metrics.tsv", sep="\t", index=False)
    metrics.loc[metrics["vertical_layer"] == "layer_1_low", ["gene"]].to_csv(
        args.outdir / "vertical_layer_1_low_genes.tsv", sep="\t", index=False
    )
    metrics.loc[metrics["vertical_layer"] == "layer_2_middle", ["gene"]].to_csv(
        args.outdir / "vertical_layer_2_middle_genes.tsv", sep="\t", index=False
    )
    metrics.loc[metrics["vertical_layer"] == "layer_3_high", ["gene"]].to_csv(
        args.outdir / "vertical_layer_3_high_genes.tsv", sep="\t", index=False
    )
    if "baseline_curve" in metrics.attrs:
        metrics.attrs["baseline_curve"].to_csv(args.outdir / "vertical_layer_baseline_curve.tsv", sep="\t", index=False)
    if "layer_centers" in metrics.attrs:
        metrics.attrs["layer_centers"].to_csv(args.outdir / "vertical_layer_kmeans_centers.tsv", sep="\t", index=False)
    if "valley_boundaries" in metrics.attrs:
        metrics.attrs["valley_boundaries"].to_csv(args.outdir / "vertical_layer_valley_boundaries.tsv", sep="\t", index=False)
        boundary_crossing_summary(metrics, [0.005, 0.01, 0.02]).to_csv(
            args.outdir / "vertical_layer_boundary_crossing_summary.tsv", sep="\t", index=False
        )
    metrics.loc[metrics["phylo_mixing_class"] == "lineage_strong", ["gene"]].to_csv(
        args.outdir / "lineage_strong_genes.tsv", sep="\t", index=False
    )
    metrics.loc[metrics["phylo_mixing_class"] == "mixed_high", ["gene"]].to_csv(
        args.outdir / "mixed_high_genes.tsv", sep="\t", index=False
    )
    plot_metrics(metrics, args.outdir / "allele_diversity_vs_same_allele_tree_distance", lineage_threshold, mixed_threshold)
    plot_vertical_layers(metrics, args.outdir / "allele_diversity_vs_same_allele_tree_distance_vertical_layers", genome_pair_mean_distance)
    hist_table = plot_same_allele_distance_stacked_hist(
        metrics,
        args.outdir / "same_allele_mean_tree_distance_1000bin_stacked_hist",
        args.hist_bins,
        genome_pair_mean_distance,
    )
    hist_table.to_csv(args.outdir / "same_allele_mean_tree_distance_1000bin_stacked_hist.tsv", sep="\t", index=False)
    layer_counts = metrics.groupby("vertical_layer", observed=True).size().rename("n_genes").reset_index()
    layer_counts.to_csv(args.outdir / "vertical_layer_summary.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_genomes": int(profile.shape[0]),
                "n_core_genes_evaluated": int(len(metrics)),
                "all_genome_pair_mean_tree_distance": genome_pair_mean_distance,
                "min_allele_count": args.min_allele_count,
                "lineage_quantile": args.lineage_quantile,
                "lineage_threshold_same_allele_tree_distance": lineage_threshold,
                "n_lineage_strong_genes": int((metrics["phylo_mixing_class"] == "lineage_strong").sum()),
                "mixed_quantile": args.mixed_quantile,
                "mixed_threshold_same_allele_tree_distance": mixed_threshold,
                "n_mixed_high_genes": int((metrics["phylo_mixing_class"] == "mixed_high").sum()),
                "vertical_layers": args.vertical_layers,
                "layer_method": "local_valley_boundary",
                "layer_bins": args.layer_bins,
                "hist_bins": args.hist_bins,
            }
        ]
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
