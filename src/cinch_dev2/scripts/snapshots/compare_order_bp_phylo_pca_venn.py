#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
from sklearn.preprocessing import StandardScaler


def edge_id(g1: str, g2: str) -> str:
    return "||".join(sorted((str(g1), str(g2))))


def short_gene_name(name: str) -> str:
    return name.split(":")[-1]


def read_thresholds(path: Path) -> dict[str, float]:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["metric"], df["threshold"]))


def build_matrix(edges: pd.DataFrame, r2_col: str) -> pd.DataFrame:
    genes = pd.Index(sorted(set(edges["gene1"]).union(edges["gene2"])))
    mat = pd.DataFrame(0.0, index=genes, columns=genes)
    for row in edges.itertuples(index=False):
        g1 = getattr(row, "gene1")
        g2 = getattr(row, "gene2")
        val = float(getattr(row, r2_col))
        mat.loc[g1, g2] = val
        mat.loc[g2, g1] = val
    np.fill_diagonal(mat.values, 0.0)
    return mat


def make_pca(mat: pd.DataFrame, random_state: int) -> tuple[pd.DataFrame, np.ndarray]:
    x = StandardScaler().fit_transform(mat.values)
    pca = PCA(n_components=2, random_state=random_state)
    xy = pca.fit_transform(x)
    degree = (mat.values > 0).sum(axis=1)
    strength = mat.values.sum(axis=1)
    coords = pd.DataFrame(
        {
            "gene": mat.index,
            "gene_short": [short_gene_name(g) for g in mat.index],
            "PCA1": xy[:, 0],
            "PCA2": xy[:, 1],
            "degree": degree,
            "r2_strength": strength,
        }
    )
    return coords, pca.explained_variance_ratio_


def plot_pca(coords: pd.DataFrame, var_ratio: np.ndarray, title: str, out_base: Path) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    points = ax.scatter(
        coords["PCA1"],
        coords["PCA2"],
        c=coords["r2_strength"],
        s=np.clip(18 + coords["degree"] * 1.8, 26, 160),
        cmap="viridis",
        alpha=0.84,
        edgecolors="white",
        linewidths=0.25,
    )
    cbar = fig.colorbar(points, ax=ax)
    cbar.set_label("r2_strength")
    for row in coords.sort_values(["degree", "r2_strength"], ascending=False).head(12).itertuples(index=False):
        ax.text(row.PCA1, row.PCA2, row.gene_short, fontsize=8, ha="left", va="bottom", alpha=0.9)
    ax.set_title(title)
    ax.set_xlabel(f"PCA1 ({var_ratio[0] * 100:.1f}%)")
    ax.set_ylabel(f"PCA2 ({var_ratio[1] * 100:.1f}%)")
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def evaluate_kmeans(coords: pd.DataFrame, x_col: str, y_col: str, k_values: list[int], random_state: int) -> pd.DataFrame:
    xy = coords[[x_col, y_col]].to_numpy(dtype=float)
    xy = StandardScaler().fit_transform(xy)
    rows = []
    for k in k_values:
        if k < 2 or k >= len(coords):
            continue
        model = KMeans(n_clusters=k, random_state=random_state, n_init=50)
        labels = model.fit_predict(xy)
        counts = pd.Series(labels).value_counts().sort_index()
        rows.append(
            {
                "k": k,
                "silhouette": silhouette_score(xy, labels),
                "calinski_harabasz": calinski_harabasz_score(xy, labels),
                "davies_bouldin": davies_bouldin_score(xy, labels),
                "inertia": model.inertia_,
                "min_cluster_size": int(counts.min()),
                "max_cluster_size": int(counts.max()),
                "cluster_sizes": ",".join(str(int(counts.get(i, 0))) for i in range(k)),
            }
        )
    if not rows:
        raise ValueError("No valid k values for KMeans evaluation")
    eval_df = pd.DataFrame(rows).sort_values(["silhouette", "calinski_harabasz"], ascending=[False, False])
    eval_df["selected"] = False
    eval_df.loc[eval_df.index[0], "selected"] = True
    return eval_df.sort_values("k")


def add_kmeans_clusters(
    coords: pd.DataFrame,
    x_col: str,
    y_col: str,
    k_values: list[int],
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    eval_df = evaluate_kmeans(coords, x_col, y_col, k_values, random_state)
    best_k = int(eval_df.loc[eval_df["selected"], "k"].iloc[0])
    xy = StandardScaler().fit_transform(coords[[x_col, y_col]].to_numpy(dtype=float))
    model = KMeans(n_clusters=best_k, random_state=random_state, n_init=50)
    clustered = coords.copy()
    clustered["kmeans_k"] = best_k
    clustered["kmeans_cluster"] = model.fit_predict(xy) + 1
    return clustered, eval_df


def plot_pca_clusters(coords: pd.DataFrame, var_ratio: np.ndarray, title: str, out_base: Path) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    palette = sns.color_palette("tab10", n_colors=int(coords["kmeans_cluster"].nunique()))
    sns.scatterplot(
        data=coords,
        x="PCA1",
        y="PCA2",
        hue="kmeans_cluster",
        palette=palette,
        s=62,
        alpha=0.86,
        edgecolor="white",
        linewidth=0.25,
        ax=ax,
    )
    for row in coords.sort_values(["degree", "r2_strength"], ascending=False).head(12).itertuples(index=False):
        ax.text(row.PCA1, row.PCA2, row.gene_short, fontsize=8, ha="left", va="bottom", alpha=0.9)
    best_k = int(coords["kmeans_k"].iloc[0])
    ax.set_title(f"{title} KMeans clusters (k={best_k})")
    ax.set_xlabel(f"PCA1 ({var_ratio[0] * 100:.1f}%)")
    ax.set_ylabel(f"PCA2 ({var_ratio[1] * 100:.1f}%)")
    ax.legend(title="cluster", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_venn(order_set: set[str], bp_set: set[str], out_base: Path) -> dict[str, int]:
    only_order = len(order_set - bp_set)
    only_bp = len(bp_set - order_set)
    both = len(order_set & bp_set)
    total = len(order_set | bp_set)

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(8, 6.5))
    circle1 = plt.Circle((0.43, 0.5), 0.32, color="#4c78a8", alpha=0.48)
    circle2 = plt.Circle((0.62, 0.5), 0.32, color="#f58518", alpha=0.48)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.text(0.28, 0.5, f"{only_order:,}", ha="center", va="center", fontsize=22, weight="bold")
    ax.text(0.525, 0.5, f"{both:,}", ha="center", va="center", fontsize=22, weight="bold")
    ax.text(0.77, 0.5, f"{only_bp:,}", ha="center", va="center", fontsize=22, weight="bold")
    ax.text(0.31, 0.86, "order + phylo", ha="center", va="center", fontsize=15)
    ax.text(0.74, 0.86, "bp + phylo", ha="center", va="center", fontsize=15)
    ax.text(0.525, 0.09, f"union = {total:,}", ha="center", va="center", fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)
    return {"only_order_phylo": only_order, "only_bp_phylo": only_bp, "overlap": both, "union": total}


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare order+phylo and bp+phylo filtered LD-like edge sets.")
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("--thresholds", type=Path, default=Path("Leca/secondary_filter_0427_1/elbow_thresholds.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1"))
    parser.add_argument("--r2-col", default="r2_adj")
    parser.add_argument("--p-col", default="p_adj")
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--kmeans-min-k", type=int, default=2)
    parser.add_argument("--kmeans-max-k", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    thresholds = read_thresholds(args.thresholds)
    df = pd.read_csv(args.input_tsv, sep="\t")
    for col in [args.r2_col, args.p_col, "mean_order_dist", "mean_bp_dist", "mean_phylo_dist"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    sig = df[(df[args.p_col] <= args.p_cutoff) & np.isfinite(df[args.r2_col])].copy()
    order_edges = sig[
        (sig["mean_order_dist"] >= thresholds["mean_order_dist"])
        & (sig["mean_phylo_dist"] >= thresholds["mean_phylo_dist"])
    ].copy()
    bp_edges = sig[
        (sig["mean_bp_dist"] >= thresholds["mean_bp_dist"])
        & (sig["mean_phylo_dist"] >= thresholds["mean_phylo_dist"])
    ].copy()

    order_edges["edge_id"] = [edge_id(a, b) for a, b in zip(order_edges["gene1"], order_edges["gene2"])]
    bp_edges["edge_id"] = [edge_id(a, b) for a, b in zip(bp_edges["gene1"], bp_edges["gene2"])]
    order_set = set(order_edges["edge_id"])
    bp_set = set(bp_edges["edge_id"])

    order_edges.to_csv(args.outdir / "order_phylo_filtered_edges.tsv.gz", sep="\t", index=False)
    bp_edges.to_csv(args.outdir / "bp_phylo_filtered_edges.tsv.gz", sep="\t", index=False)

    order_mat = build_matrix(order_edges, args.r2_col)
    bp_mat = build_matrix(bp_edges, args.r2_col)
    order_mat.to_csv(args.outdir / "order_phylo_gene_by_gene_r2_matrix.tsv", sep="\t")
    bp_mat.to_csv(args.outdir / "bp_phylo_gene_by_gene_r2_matrix.tsv", sep="\t")

    order_coords, order_var = make_pca(order_mat, args.random_state)
    bp_coords, bp_var = make_pca(bp_mat, args.random_state)
    k_values = list(range(args.kmeans_min_k, args.kmeans_max_k + 1))
    order_coords, order_cluster_eval = add_kmeans_clusters(order_coords, "PCA1", "PCA2", k_values, args.random_state)
    bp_coords, bp_cluster_eval = add_kmeans_clusters(bp_coords, "PCA1", "PCA2", k_values, args.random_state)
    order_coords.to_csv(args.outdir / "order_phylo_pca_coords.tsv", sep="\t", index=False)
    bp_coords.to_csv(args.outdir / "bp_phylo_pca_coords.tsv", sep="\t", index=False)
    order_cluster_eval.to_csv(args.outdir / "order_phylo_pca_kmeans_eval.tsv", sep="\t", index=False)
    bp_cluster_eval.to_csv(args.outdir / "bp_phylo_pca_kmeans_eval.tsv", sep="\t", index=False)
    plot_pca(order_coords, order_var, "Order + phylo filtered r2 profile PCA", args.outdir / "pca_order_phylo")
    plot_pca(bp_coords, bp_var, "bp + phylo filtered r2 profile PCA", args.outdir / "pca_bp_phylo")
    plot_pca_clusters(
        order_coords,
        order_var,
        "Order + phylo filtered r2 profile PCA",
        args.outdir / "pca_order_phylo_kmeans",
    )
    plot_pca_clusters(
        bp_coords,
        bp_var,
        "bp + phylo filtered r2 profile PCA",
        args.outdir / "pca_bp_phylo_kmeans",
    )

    venn_counts = plot_venn(order_set, bp_set, args.outdir / "venn_order_phylo_vs_bp_phylo")
    summary = {
        "input_pairs": len(df),
        "significant_pairs": len(sig),
        "order_phylo_edges": len(order_edges),
        "bp_phylo_edges": len(bp_edges),
        "order_phylo_genes": len(order_mat),
        "bp_phylo_genes": len(bp_mat),
        "order_phylo_pca1_var": order_var[0],
        "order_phylo_pca2_var": order_var[1],
        "bp_phylo_pca1_var": bp_var[0],
        "bp_phylo_pca2_var": bp_var[1],
        "order_phylo_pca_kmeans_k": int(order_cluster_eval.loc[order_cluster_eval["selected"], "k"].iloc[0]),
        "order_phylo_pca_kmeans_silhouette": float(
            order_cluster_eval.loc[order_cluster_eval["selected"], "silhouette"].iloc[0]
        ),
        "bp_phylo_pca_kmeans_k": int(bp_cluster_eval.loc[bp_cluster_eval["selected"], "k"].iloc[0]),
        "bp_phylo_pca_kmeans_silhouette": float(
            bp_cluster_eval.loc[bp_cluster_eval["selected"], "silhouette"].iloc[0]
        ),
        **venn_counts,
        "order_threshold": thresholds["mean_order_dist"],
        "bp_threshold": thresholds["mean_bp_dist"],
        "phylo_threshold": thresholds["mean_phylo_dist"],
    }
    pd.Series(summary).to_csv(args.outdir / "summary.tsv", sep="\t", header=False)
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
