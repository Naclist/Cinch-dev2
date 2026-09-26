#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import normalized_mutual_info_score, silhouette_score

from spydrpick_locus_wg_nmi import aracne_prune


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def weighted_mi(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    total = weights.sum()
    if total <= 0:
        return 0.0
    joint = np.array(
        [
            weights[(x == 0) & (y == 0)].sum(),
            weights[(x == 0) & (y == 1)].sum(),
            weights[(x == 1) & (y == 0)].sum(),
            weights[(x == 1) & (y == 1)].sum(),
        ],
        dtype=float,
    ).reshape(2, 2)
    joint /= total
    px = joint.sum(axis=1)
    py = joint.sum(axis=0)
    expected = px[:, None] * py[None, :]
    mask = joint > 0
    return float(np.sum(joint[mask] * np.log(joint[mask] / expected[mask])))


def conditional_mi(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    clusters: np.ndarray,
    eligible_clusters: list[int] | None = None,
) -> float:
    total = weights.sum()
    value = 0.0
    cluster_values = np.unique(clusters) if eligible_clusters is None else eligible_clusters
    for cluster in cluster_values:
        mask = clusters == cluster
        cluster_weight = weights[mask].sum()
        if cluster_weight > 0:
            value += cluster_weight / total * weighted_mi(x[mask], y[mask], weights[mask])
    return value


def informative_clusters(
    x: np.ndarray,
    y: np.ndarray,
    clusters: np.ndarray,
    minimum_size: int,
) -> tuple[list[int], int]:
    informative = []
    positive = 0
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        size = int(mask.sum())
        if size < minimum_size:
            continue
        sx = int(x[mask].sum())
        sy = int(y[mask].sum())
        if not (0 < sx < size and 0 < sy < size):
            continue
        informative.append(int(cluster))
        n11 = int(np.sum(x[mask] & y[mask]))
        n10 = int(np.sum(x[mask] & ~y[mask]))
        n01 = int(np.sum(~x[mask] & y[mask]))
        n00 = int(np.sum(~x[mask] & ~y[mask]))
        if n11 * n00 > n10 * n01:
            positive += 1
    return informative, positive


def phylo_permutation(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    clusters: np.ndarray,
    eligible_clusters: list[int],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float]:
    observed = conditional_mi(x, y, weights, clusters, eligible_clusters)
    null = np.empty(permutations, dtype=float)
    eligible_mask = np.isin(clusters, eligible_clusters)
    for permutation in range(permutations):
        shuffled = y.copy()
        for cluster in eligible_clusters:
            indices = np.flatnonzero(clusters == cluster)
            shuffled[indices] = rng.permutation(shuffled[indices])
        null[permutation] = conditional_mi(
            x[eligible_mask],
            shuffled[eligible_mask],
            weights[eligible_mask],
            clusters[eligible_mask],
            eligible_clusters,
        )
    pvalue = (1.0 + np.sum(null >= observed)) / (permutations + 1.0)
    standard_deviation = float(null.std(ddof=1))
    zscore = (observed - float(null.mean())) / standard_deviation if standard_deviation > 0 else np.inf
    return observed, float(null.mean()), zscore, float(pvalue)


def build_hiercc(
    similarity: np.ndarray,
    scheme_loci: int,
    singleton_limit: float,
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray, float, tuple[float, float]]:
    distance = scheme_loci * (1.0 - similarity.astype(float))
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method="single")
    thresholds = np.unique(
        np.r_[0, 1, 2, 5, 10, 20, 50, 100, 200, 400, np.rint(tree[:, 2])]
    )
    thresholds = thresholds[(thresholds >= 0) & (thresholds < scheme_loci)]
    labels = [fcluster(tree, threshold, criterion="distance") for threshold in thresholds]

    n_levels = len(thresholds)
    partition_nmi = np.eye(n_levels)
    for left in range(n_levels):
        for right in range(left + 1, n_levels):
            value = normalized_mutual_info_score(labels[left], labels[right])
            partition_nmi[left, right] = partition_nmi[right, left] = value

    blocks = []
    start = 0
    while start < n_levels:
        end = start
        while end + 1 < n_levels and partition_nmi[start : end + 2, start : end + 2].min() >= 0.9:
            end += 1
        if end - start + 1 >= 3:
            blocks.append((start, end))
        start = max(start + 1, end + 1)

    metrics = []
    for threshold, cluster_labels in zip(thresholds, labels):
        sizes = pd.Series(cluster_labels).value_counts()
        metrics.append(
            {
                "hc_threshold": float(threshold),
                "n_clusters": len(sizes),
                "singleton_fraction": float((sizes == 1).sum() / len(cluster_labels)),
                "largest_cluster": int(sizes.max()),
                "shannon": float(
                    -np.sum((sizes / sizes.sum()) * np.log(sizes / sizes.sum())) / np.log(len(cluster_labels))
                ),
            }
        )
    metrics = pd.DataFrame(metrics)

    selected_block = None
    for start, end in blocks:
        if metrics.iloc[start]["singleton_fraction"] <= singleton_limit:
            selected_block = (start, end)
            break
    if selected_block is None:
        raise RuntimeError("No stable HierCC block passed the singleton criterion")

    start, end = selected_block
    best = None
    for index in range(start, end + 1):
        cluster_labels = labels[index]
        if len(np.unique(cluster_labels)) in (1, len(cluster_labels)):
            continue
        score = silhouette_score(distance, cluster_labels, metric="precomputed")
        metrics.loc[index, "silhouette"] = score
        if best is None or score > best[0]:
            best = (score, index)
    selected_score, selected_index = best
    metrics["selected_shc"] = False
    metrics.loc[selected_index, "selected_shc"] = True
    metrics["stable_block"] = ""
    for block_number, (left, right) in enumerate(blocks, start=1):
        metrics.loc[left:right, "stable_block"] = f"block_{block_number}"
    return (
        labels[selected_index],
        metrics,
        partition_nmi,
        float(thresholds[selected_index]),
        (float(thresholds[start]), float(thresholds[end])),
    )


def analyse_edges(
    edge_path: Path,
    selection: str,
    profile: pd.DataFrame,
    weights: np.ndarray,
    clusters: np.ndarray,
    permutations: int,
    minimum_cluster_size: int,
    minimum_informative_clusters: int,
    outdir: Path,
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    edges = pd.read_csv(edge_path, sep="\t")
    genes = sorted(set(edges["gene1"]) | set(edges["gene2"]))
    binary = (profile[genes].to_numpy() > 0)
    index = {gene: position for position, gene in enumerate(genes)}
    records = []
    rng = np.random.default_rng(seed)

    for row in edges.itertuples(index=False):
        x = binary[:, index[row.gene1]]
        y = binary[:, index[row.gene2]]
        informative, positive = informative_clusters(x, y, clusters, minimum_cluster_size)
        record = {
            "gene1": row.gene1,
            "gene2": row.gene2,
            "weighted_mi": row.weighted_mi,
            "mean_order_dist": row.mean_order_dist,
            "mean_bp_dist": row.mean_bp_dist,
            "informative_shc": len(informative),
            "positive_shc": positive,
            "eligible_for_permutation": (
                len(informative) >= minimum_informative_clusters
                and positive >= minimum_informative_clusters
            ),
        }
        if record["eligible_for_permutation"]:
            observed, null_mean, zscore, pvalue = phylo_permutation(
                x,
                y,
                weights,
                clusters,
                informative,
                permutations,
                rng,
            )
            record.update(
                {
                    "conditional_mi": observed,
                    "conditional_mi_null_mean": null_mean,
                    "conditional_mi_zscore": zscore,
                    "phylo_perm_p": pvalue,
                }
            )
        records.append(record)

    results = pd.DataFrame(records)
    results["phylo_perm_q"] = np.nan
    eligible = results["eligible_for_permutation"]
    if eligible.any():
        results.loc[eligible, "phylo_perm_q"] = bh_adjust(
            results.loc[eligible, "phylo_perm_p"].to_numpy()
        )
    results["passes_phylo"] = (
        results["eligible_for_permutation"]
        & (results["phylo_perm_q"] < 0.05)
        & (results["conditional_mi_zscore"] > 0)
    )
    results.to_csv(outdir / f"{selection}_shc_phylo_test_all_candidates.tsv.gz", sep="\t", index=False)

    passed = results[results["passes_phylo"]].copy()
    if len(passed):
        passed = aracne_prune(passed, "conditional_mi")
    else:
        passed["aracne_direct"] = pd.Series(dtype=bool)
    passed.to_csv(outdir / f"{selection}_shc_strict_edges_aracne.tsv", sep="\t", index=False)
    direct = passed[passed["aracne_direct"]] if len(passed) else passed
    direct.to_csv(outdir / f"{selection}_shc_strict_direct_edges.tsv", sep="\t", index=False)

    summary = {
        "selection": selection,
        "distance_candidates": len(results),
        "eligible_two_positive_informative_shc": int(results["eligible_for_permutation"].sum()),
        "phylo_q_lt_0.05": int(results["passes_phylo"].sum()),
        "strict_direct_after_aracne": len(direct),
        "strict_indirect_after_aracne": int(len(passed) - len(direct)),
    }
    return results, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="HierCC-like SHC correction for weighted MI candidate edges.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument(
        "--similarity",
        type=Path,
        default=Path(
            "Leca/results/profile_similarity/wg_cg_profile_0427_1/wg_cg_similarity_matrices.npz"
        ),
    )
    parser.add_argument(
        "--weights",
        type=Path,
        default=Path(
            "Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv"
        ),
    )
    parser.add_argument(
        "--candidate-dir",
        type=Path,
        default=Path(
            "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_nmi_distance_comparison"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(
            "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction"
        ),
    )
    parser.add_argument("--scheme-loci", type=int, default=3632)
    parser.add_argument("--singleton-limit", type=float, default=0.20)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--minimum-cluster-size", type=int, default=5)
    parser.add_argument("--minimum-informative-clusters", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260709)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    matrices = np.load(args.similarity)
    genomes = matrices["genomes"].astype(str)
    similarity = matrices["cg_allele_concordance_callable"]
    clusters, metrics, partition_nmi, shc_threshold, stable_block = build_hiercc(
        similarity, args.scheme_loci, args.singleton_limit
    )
    metrics.to_csv(args.outdir / "hiercc_level_metrics.tsv", sep="\t", index=False)
    np.save(args.outdir / "hiercc_partition_nmi.npy", partition_nmi)
    pd.DataFrame({"genome": genomes, "shc": clusters}).to_csv(
        args.outdir / "genome_shc_assignments.tsv", sep="\t", index=False
    )

    profile = pd.read_csv(args.profile, sep="\t", index_col=0).loc[genomes]
    weight_table = pd.read_csv(args.weights, sep="\t").set_index("genome")
    weights = weight_table.loc[genomes, "wg_content_weight"].to_numpy(dtype=float)

    summaries = []
    results = {}
    for offset, selection in enumerate(["order_dist", "bp_dist"]):
        result, summary = analyse_edges(
            args.candidate_dir / f"{selection}_mi_upper_tail_outliers_aracne.tsv",
            selection,
            profile,
            weights,
            clusters,
            args.permutations,
            args.minimum_cluster_size,
            args.minimum_informative_clusters,
            args.outdir,
            args.seed + offset,
        )
        results[selection] = result
        summaries.append(summary)
    summary_table = pd.DataFrame(summaries)
    summary_table.insert(0, "shc_threshold", shc_threshold)
    summary_table.insert(1, "stable_block_start", stable_block[0])
    summary_table.insert(2, "stable_block_end", stable_block[1])
    summary_table.to_csv(args.outdir / "shc_filter_summary.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    sns.heatmap(partition_nmi, cmap="viridis", vmin=0.8, vmax=1.0, ax=axes[0, 0], cbar_kws={"label": "Partition NMI"})
    axes[0, 0].set_title("A. HierCC partition stability")
    axes[0, 0].set_xlabel("HC level index")
    axes[0, 0].set_ylabel("HC level index")

    axes[0, 1].plot(metrics["hc_threshold"], metrics["singleton_fraction"], label="Singleton fraction")
    silhouette = metrics.dropna(subset=["silhouette"])
    axes[0, 1].plot(silhouette["hc_threshold"], silhouette["silhouette"], label="Silhouette")
    axes[0, 1].axvline(shc_threshold, color="#C45D3A", linestyle=":", label=f"SHC=HC{shc_threshold:.0f}")
    axes[0, 1].set_title("B. Selection of the first population-scale stable block")
    axes[0, 1].set_xlabel("HC allelic threshold")
    axes[0, 1].legend(frameon=False)

    x = np.arange(len(summaries))
    stages = [
        ("distance_candidates", "Distance candidates"),
        ("eligible_two_positive_informative_shc", "Replicated within SHCs"),
        ("phylo_q_lt_0.05", "Permutation q<0.05"),
        ("strict_direct_after_aracne", "Strict direct"),
    ]
    width = 0.18
    for stage_index, (column, label) in enumerate(stages):
        axes[1, 0].bar(
            x + (stage_index - 1.5) * width,
            [summary[column] for summary in summaries],
            width,
            label=label,
        )
    axes[1, 0].set_xticks(x, [summary["selection"] for summary in summaries])
    axes[1, 0].set_title("C. Edge-filtering waterfall")
    axes[1, 0].set_ylabel("Pairs")
    axes[1, 0].legend(frameon=False, fontsize=8)

    colors = {"order_dist": "#4267B2", "bp_dist": "#267A6B"}
    for selection, result in results.items():
        eligible = result[result["eligible_for_permutation"]]
        axes[1, 1].scatter(
            eligible["weighted_mi"],
            eligible["conditional_mi"],
            s=12,
            alpha=0.5,
            label=selection,
            color=colors[selection],
        )
    axes[1, 1].set_title("D. Raw MI versus SHC-conditional MI")
    axes[1, 1].set_xlabel("Weighted MI")
    axes[1, 1].set_ylabel("Conditional MI")
    axes[1, 1].legend(frameon=False)
    sns.despine(fig=fig)
    fig.suptitle(
        "HierCC-like smallest stable cluster correction of wgMLST MI edges",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(args.outdir / "hiercc_shc_correction_diagnostics.png", dpi=260)
    fig.savefig(args.outdir / "hiercc_shc_correction_diagnostics.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
