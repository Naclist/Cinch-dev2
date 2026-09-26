#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap


def short_gene_name(name: str) -> str:
    return name.split(":")[-1]


def build_r2_matrix(df: pd.DataFrame, p_cutoff: float, r2_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"gene1", "gene2", "p_adj", r2_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")

    work = df.copy()
    work["p_adj"] = pd.to_numeric(work["p_adj"], errors="coerce")
    work[r2_col] = pd.to_numeric(work[r2_col], errors="coerce")
    sig = work[(work["p_adj"] <= p_cutoff) & np.isfinite(work[r2_col])].copy()
    if sig.empty:
        raise ValueError(f"No significant pairs found with p_adj <= {p_cutoff}")

    genes = pd.Index(sorted(set(sig["gene1"]).union(sig["gene2"])))
    mat = pd.DataFrame(0.0, index=genes, columns=genes)

    for row in sig.itertuples(index=False):
        g1 = getattr(row, "gene1")
        g2 = getattr(row, "gene2")
        val = float(getattr(row, r2_col))
        mat.loc[g1, g2] = val
        mat.loc[g2, g1] = val

    np.fill_diagonal(mat.values, 0.0)
    return mat, sig


def make_embeddings(mat: pd.DataFrame, random_state: int) -> pd.DataFrame:
    x = mat.values
    x_scaled = StandardScaler().fit_transform(x)

    pca = PCA(n_components=2, random_state=random_state)
    pca_xy = pca.fit_transform(x_scaled)

    n_neighbors = min(15, max(2, len(mat) - 1))
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.15,
        metric="euclidean",
        random_state=random_state,
    )
    umap_xy = reducer.fit_transform(x_scaled)

    degree = (mat.values > 0).sum(axis=1)
    strength = mat.values.sum(axis=1)
    max_r2 = mat.values.max(axis=1)

    out = pd.DataFrame(
        {
            "gene": mat.index,
            "gene_short": [short_gene_name(g) for g in mat.index],
            "PCA1": pca_xy[:, 0],
            "PCA2": pca_xy[:, 1],
            "UMAP1": umap_xy[:, 0],
            "UMAP2": umap_xy[:, 1],
            "degree": degree,
            "r2_strength": strength,
            "max_r2_adj": max_r2,
            "pca1_var_ratio": pca.explained_variance_ratio_[0],
            "pca2_var_ratio": pca.explained_variance_ratio_[1],
        }
    )
    return out


def plot_embedding(
    emb: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    title: str,
    out_base: Path,
    xlabel: str | None = None,
    ylabel: str | None = None,
    label_top: int = 15,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    pts = ax.scatter(
        emb[x_col],
        emb[y_col],
        c=emb[color_col],
        s=np.clip(20 + emb["degree"] * 3, 28, 180),
        cmap="viridis",
        alpha=0.82,
        linewidths=0.25,
        edgecolors="white",
    )
    cbar = fig.colorbar(pts, ax=ax)
    cbar.set_label(color_col)

    if label_top > 0:
        label_df = emb.sort_values(["degree", "r2_strength"], ascending=False).head(label_top)
        for row in label_df.itertuples(index=False):
            ax.text(
                getattr(row, x_col),
                getattr(row, y_col),
                row.gene_short,
                fontsize=8,
                ha="left",
                va="bottom",
                alpha=0.9,
            )

    ax.set_title(title)
    ax.set_xlabel(xlabel or x_col)
    ax.set_ylabel(ylabel or y_col)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build an r2_adj profile matrix from significant ld-like gene pairs and plot PCA/UMAP."
    )
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("Leca/r2_significant_ordination"))
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--r2-col", default="r2_adj", choices=["r2_adj", "r2_raw"])
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--label-top", type=int, default=15)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input_tsv, sep="\t")
    mat, sig = build_r2_matrix(df, args.p_cutoff, args.r2_col)
    emb = make_embeddings(mat, args.random_state)

    mat.to_csv(args.outdir / f"significant_gene_by_gene_{args.r2_col}_matrix.tsv", sep="\t")
    sig.to_csv(args.outdir / "significant_pairs.tsv", sep="\t", index=False)
    emb.to_csv(args.outdir / "gene_r2_profile_pca_umap.tsv", sep="\t", index=False)

    pca1 = emb["pca1_var_ratio"].iloc[0] * 100
    pca2 = emb["pca2_var_ratio"].iloc[0] * 100
    plot_embedding(
        emb,
        "PCA1",
        "PCA2",
        "r2_strength",
        f"Significant {args.r2_col} profile PCA",
        args.outdir / "pca_r2_profile_by_strength",
        xlabel=f"PCA1 ({pca1:.1f}%)",
        ylabel=f"PCA2 ({pca2:.1f}%)",
        label_top=args.label_top,
    )
    plot_embedding(
        emb,
        "PCA1",
        "PCA2",
        "degree",
        f"Significant {args.r2_col} profile PCA",
        args.outdir / "pca_r2_profile_by_degree",
        xlabel=f"PCA1 ({pca1:.1f}%)",
        ylabel=f"PCA2 ({pca2:.1f}%)",
        label_top=args.label_top,
    )
    plot_embedding(
        emb,
        "UMAP1",
        "UMAP2",
        "r2_strength",
        f"Significant {args.r2_col} profile UMAP",
        args.outdir / "umap_r2_profile_by_strength",
        label_top=args.label_top,
    )
    plot_embedding(
        emb,
        "UMAP1",
        "UMAP2",
        "degree",
        f"Significant {args.r2_col} profile UMAP",
        args.outdir / "umap_r2_profile_by_degree",
        label_top=args.label_top,
    )

    summary = {
        "input_rows": len(df),
        "significant_pairs": len(sig),
        "genes_in_significant_pairs": len(mat),
        "p_cutoff": args.p_cutoff,
        "r2_column": args.r2_col,
        "pca1_variance_ratio": emb["pca1_var_ratio"].iloc[0],
        "pca2_variance_ratio": emb["pca2_var_ratio"].iloc[0],
        "mean_degree": float(emb["degree"].mean()),
        "median_degree": float(emb["degree"].median()),
        "max_degree": int(emb["degree"].max()),
    }
    pd.Series(summary).to_csv(args.outdir / "summary.tsv", sep="\t", header=False)
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
