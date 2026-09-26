#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.manifold import MDS
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, pairwise_distances, silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD


def short_gene_name(name: str) -> str:
    return str(name).split(":")[-1]


def build_plot_meta(mat: pd.DataFrame, coords: np.ndarray, x_col: str, y_col: str) -> pd.DataFrame:
    values = mat.values
    return pd.DataFrame(
        {
            "gene": mat.index,
            "gene_short": [short_gene_name(g) for g in mat.index],
            x_col: coords[:, 0],
            y_col: coords[:, 1],
            "degree": (values > 0).sum(axis=1),
            "r2_strength": values.sum(axis=1),
        }
    )


def classical_pcoa(distance_matrix: np.ndarray, n_components: int = 2) -> np.ndarray:
    n = distance_matrix.shape[0]
    d2 = np.square(distance_matrix)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    return eigvecs[:, :n_components] * np.sqrt(np.maximum(eigvals[:n_components], 0))


def fit_vectors(coords: pd.DataFrame, assoc: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    merged = coords[["gene", x_col, y_col]].merge(
        assoc.reset_index().rename(columns={"index": "gene"}),
        on="gene",
        how="inner",
    )
    x = merged[x_col].to_numpy(dtype=float)
    y = merged[y_col].to_numpy(dtype=float)
    vectors = []
    for col in assoc.columns:
        z = merged[col].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        if valid.sum() < 5 or np.nanstd(z[valid]) == 0:
            continue
        rx = np.corrcoef(x[valid], z[valid])[0, 1]
        ry = np.corrcoef(y[valid], z[valid])[0, 1]
        vectors.append(
            {
                "metadata": col,
                "axis1_corr": rx,
                "axis2_corr": ry,
                "r2_on_axes": rx * rx + ry * ry,
                "n_genes": int(valid.sum()),
                "mean_score": float(np.nanmean(z[valid])),
                "max_score": float(np.nanmax(z[valid])),
            }
        )
    return pd.DataFrame(vectors).sort_values("r2_on_axes", ascending=False)


def plot_ordination(
    coords: pd.DataFrame,
    vectors: pd.DataFrame | None,
    x_col: str,
    y_col: str,
    title: str,
    out_base: Path,
    top_n: int = 16,
    arrow_scale: float = 0.78,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9.2, 7.6))
    pts = ax.scatter(
        coords[x_col],
        coords[y_col],
        c=coords["r2_strength"] if "r2_strength" in coords.columns else np.zeros(len(coords)),
        s=np.clip(20 + coords.get("degree", pd.Series(1, index=coords.index)) * 1.5, 24, 150),
        cmap="viridis",
        alpha=0.62,
        edgecolors="white",
        linewidths=0.18,
    )
    cbar = fig.colorbar(pts, ax=ax)
    cbar.set_label("r2_strength")

    if vectors is not None and not vectors.empty:
        x_span = coords[x_col].max() - coords[x_col].min()
        y_span = coords[y_col].max() - coords[y_col].min()
        scale = arrow_scale * min(x_span, y_span)
        for row in vectors.head(top_n).itertuples(index=False):
            dx = row.axis1_corr * scale
            dy = row.axis2_corr * scale
            ax.arrow(
                0,
                0,
                dx,
                dy,
                color="#c23b22",
                width=0.0025 * min(x_span, y_span),
                head_width=0.026 * min(x_span, y_span),
                length_includes_head=True,
                alpha=0.9,
                zorder=5,
            )
            ax.text(dx * 1.08, dy * 1.08, row.metadata, color="#8f1d14", fontsize=9, ha="center", va="center")

    ax.axhline(0, color="#777777", linewidth=0.8, alpha=0.35)
    ax.axvline(0, color="#777777", linewidth=0.8, alpha=0.35)
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
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


def plot_ordination_clusters(coords: pd.DataFrame, x_col: str, y_col: str, title: str, out_base: Path) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9.2, 7.6))
    palette = sns.color_palette("tab10", n_colors=int(coords["kmeans_cluster"].nunique()))
    sns.scatterplot(
        data=coords,
        x=x_col,
        y=y_col,
        hue="kmeans_cluster",
        palette=palette,
        s=62,
        alpha=0.86,
        edgecolor="white",
        linewidth=0.25,
        ax=ax,
    )
    for row in coords.sort_values(["degree", "r2_strength"], ascending=False).head(12).itertuples(index=False):
        ax.text(getattr(row, x_col), getattr(row, y_col), row.gene_short, fontsize=8, ha="left", va="bottom", alpha=0.9)
    best_k = int(coords["kmeans_k"].iloc[0])
    ax.axhline(0, color="#777777", linewidth=0.8, alpha=0.35)
    ax.axvline(0, color="#777777", linewidth=0.8, alpha=0.35)
    ax.set_title(f"{title} KMeans clusters (k={best_k})")
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.legend(title="cluster", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def run_nmds(mat: pd.DataFrame, random_state: int, max_iter: int, n_init: int) -> tuple[pd.DataFrame, float]:
    x = StandardScaler().fit_transform(mat.values)
    dist = pairwise_distances(x, metric="cosine")
    dist = np.nan_to_num(dist, nan=0.0, posinf=1.0, neginf=0.0)
    model = MDS(
        n_components=2,
        metric_mds=False,
        metric="precomputed",
        random_state=random_state,
        max_iter=max_iter,
        n_init=n_init,
        normalized_stress="auto",
        n_jobs=1,
        eps=1e-3,
    )
    init = classical_pcoa(dist, 2)
    xy = model.fit_transform(dist, init=init)
    xy = xy - xy.mean(axis=0)
    return build_plot_meta(mat, xy, "NMDS1", "NMDS2"), float(model.stress_)


def constrained_rda(mat: pd.DataFrame, assoc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    common = mat.index.intersection(assoc.index)
    y = mat.loc[common].values.astype(np.float64)
    x = assoc.loc[common].values.astype(np.float64)

    y_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(y)
    x_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(x)

    q, _ = np.linalg.qr(x_scaled)
    y_hat = q @ (q.T @ y_scaled)

    svd = TruncatedSVD(n_components=2, random_state=42)
    site_scores = svd.fit_transform(y_hat)
    coords = build_plot_meta(mat.loc[common], site_scores, "RDA1", "RDA2")

    total_ss = float(np.sum(y_scaled * y_scaled))
    constrained_ss = float(np.sum(y_hat * y_hat))
    axis_var = svd.explained_variance_ratio_
    stats = {
        "total_variance_explained_by_constraints": constrained_ss / total_ss if total_ss else np.nan,
        "rda1_fraction_of_total_variance": float(axis_var[0]),
        "rda2_fraction_of_total_variance": float(axis_var[1]),
        "genes_used": int(len(common)),
        "constraints_used": int(assoc.shape[1]),
    }

    vectors = fit_vectors(coords, assoc.loc[common], "RDA1", "RDA2")
    return coords, vectors, stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Run NMDS and constrained RDA for r2-profile gene networks.")
    parser.add_argument("--order-matrix", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/order_phylo_gene_by_gene_r2_matrix.tsv"))
    parser.add_argument("--bp-matrix", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/bp_phylo_gene_by_gene_r2_matrix.tsv"))
    parser.add_argument("--assoc", type=Path, default=Path("Leca/pcoa_metadata_0427_1/gene_by_metadata_expanded_assoc_matrix.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/nmds_rda_0427_1"))
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--nmds-max-iter", type=int, default=300)
    parser.add_argument("--nmds-n-init", type=int, default=2)
    parser.add_argument("--kmeans-min-k", type=int, default=2)
    parser.add_argument("--kmeans-max-k", type=int, default=10)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    assoc = pd.read_csv(args.assoc, sep="\t", index_col=0)
    k_values = list(range(args.kmeans_min_k, args.kmeans_max_k + 1))

    summaries = []
    for label, matrix_path in [("order_phylo", args.order_matrix), ("bp_phylo", args.bp_matrix)]:
        mat = pd.read_csv(matrix_path, sep="\t", index_col=0)

        rda_coords, rda_vectors, stats = constrained_rda(mat, assoc)
        rda_coords, rda_cluster_eval = add_kmeans_clusters(rda_coords, "RDA1", "RDA2", k_values, args.random_state)
        rda_coords.to_csv(args.outdir / f"{label}_constrained_rda_coords.tsv", sep="\t", index=False)
        rda_vectors.to_csv(args.outdir / f"{label}_constrained_rda_vectors.tsv", sep="\t", index=False)
        rda_cluster_eval.to_csv(args.outdir / f"{label}_constrained_rda_kmeans_eval.tsv", sep="\t", index=False)
        plot_ordination(
            rda_coords,
            rda_vectors,
            "RDA1",
            "RDA2",
            f"{label} constrained RDA ordination",
            args.outdir / f"{label}_constrained_rda",
        )
        plot_ordination_clusters(
            rda_coords,
            "RDA1",
            "RDA2",
            f"{label} constrained RDA ordination",
            args.outdir / f"{label}_constrained_rda_kmeans",
        )

        nmds_coords, stress = run_nmds(mat, args.random_state, args.nmds_max_iter, args.nmds_n_init)
        nmds_vectors = fit_vectors(nmds_coords, assoc, "NMDS1", "NMDS2")
        nmds_coords.to_csv(args.outdir / f"{label}_nmds_coords.tsv", sep="\t", index=False)
        nmds_vectors.to_csv(args.outdir / f"{label}_nmds_vectors.tsv", sep="\t", index=False)
        plot_ordination(
            nmds_coords,
            nmds_vectors,
            "NMDS1",
            "NMDS2",
            f"{label} NMDS with metadata arrows",
            args.outdir / f"{label}_nmds_metadata_arrows",
        )

        summaries.append(
            {
                "dataset": label,
                "nmds_stress": stress,
                "nmds_top_metadata": nmds_vectors.iloc[0]["metadata"] if len(nmds_vectors) else "",
                "nmds_top_r2_on_axes": nmds_vectors.iloc[0]["r2_on_axes"] if len(nmds_vectors) else np.nan,
                "rda_top_metadata": rda_vectors.iloc[0]["metadata"] if len(rda_vectors) else "",
                "rda_top_r2_on_axes": rda_vectors.iloc[0]["r2_on_axes"] if len(rda_vectors) else np.nan,
                "rda_kmeans_k": int(rda_cluster_eval.loc[rda_cluster_eval["selected"], "k"].iloc[0]),
                "rda_kmeans_silhouette": float(
                    rda_cluster_eval.loc[rda_cluster_eval["selected"], "silhouette"].iloc[0]
                ),
                "rda_kmeans_calinski_harabasz": float(
                    rda_cluster_eval.loc[rda_cluster_eval["selected"], "calinski_harabasz"].iloc[0]
                ),
                "rda_kmeans_davies_bouldin": float(
                    rda_cluster_eval.loc[rda_cluster_eval["selected"], "davies_bouldin"].iloc[0]
                ),
                **stats,
            }
        )

    summary = pd.DataFrame(summaries)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
