#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_coords(path: Path, x_col: str, y_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"gene", x_col, y_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    return df


def get_axis_labels(df: pd.DataFrame, x_col: str, y_col: str) -> tuple[str, str]:
    if x_col.startswith("PCA") and {"PCA1", "PCA2"}.issubset(df.columns):
        return x_col, y_col
    return x_col, y_col


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
        r2 = rx * rx + ry * ry
        vectors.append(
            {
                "metadata": col,
                "axis1_corr": rx,
                "axis2_corr": ry,
                "r2_on_axes": r2,
                "n_genes": int(valid.sum()),
                "mean_nmi": float(np.nanmean(z[valid])),
                "max_nmi": float(np.nanmax(z[valid])),
            }
        )

    vec = pd.DataFrame(vectors).sort_values("r2_on_axes", ascending=False)
    return vec


def plot_with_vectors(
    coords: pd.DataFrame,
    vectors: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    out_base: Path,
    color_col: str = "r2_strength",
    top_n: int = 18,
    arrow_scale: float = 0.8,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9, 7.5))

    plot_color = coords[color_col] if color_col in coords.columns else np.zeros(len(coords))
    size_col = "degree" if "degree" in coords.columns else color_col
    plot_size = coords[size_col] if size_col in coords.columns else np.ones(len(coords))
    pts = ax.scatter(
        coords[x_col],
        coords[y_col],
        c=plot_color,
        s=np.clip(20 + plot_size * 1.6, 24, 150),
        cmap="viridis",
        alpha=0.62,
        edgecolors="white",
        linewidths=0.18,
    )
    cbar = fig.colorbar(pts, ax=ax)
    cbar.set_label(color_col if color_col in coords.columns else "value")

    x_span = coords[x_col].max() - coords[x_col].min()
    y_span = coords[y_col].max() - coords[y_col].min()
    # Arrow length is a visual scaling of axis correlations; the real fitted
    # strength is reported in the companion *_vectors.tsv files.
    scale = arrow_scale * min(x_span, y_span)
    vec = vectors.head(top_n).copy()

    for row in vec.itertuples(index=False):
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
        ax.text(
            dx * 1.08,
            dy * 1.08,
            row.metadata,
            color="#8f1d14",
            fontsize=10,
            ha="center",
            va="center",
            zorder=6,
        )

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


def process_plot(
    coords_path: Path,
    assoc: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    out_base: Path,
) -> pd.DataFrame:
    coords = load_coords(coords_path, x_col, y_col)
    vectors = fit_vectors(coords, assoc, x_col, y_col)
    plot_with_vectors(coords, vectors, x_col, y_col, title, out_base)
    vectors.to_csv(out_base.with_name(out_base.name + "_vectors.tsv"), sep="\t", index=False)
    return vectors


def main() -> None:
    parser = argparse.ArgumentParser(description="Overlay gene-metadata RDA/envfit-style arrows on PCA/PCoA plots.")
    parser.add_argument("--assoc", type=Path, default=Path("Leca/pcoa_metadata_0427_1/gene_by_metadata_nmi_matrix.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/rda_arrows_0427_1"))
    parser.add_argument("--order-pca", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/order_phylo_pca_coords.tsv"))
    parser.add_argument("--bp-pca", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/bp_phylo_pca_coords.tsv"))
    parser.add_argument("--order-pcoa", type=Path, default=Path("Leca/pcoa_metadata_0427_1/order_phylo_r2_profile_pcoa_coords.tsv"))
    parser.add_argument("--bp-pcoa", type=Path, default=Path("Leca/pcoa_metadata_0427_1/bp_phylo_r2_profile_pcoa_coords.tsv"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    assoc = pd.read_csv(args.assoc, sep="\t", index_col=0)
    assoc.index.name = "gene"

    results = []
    for label, path, x_col, y_col, title, out_name in [
        ("order_phylo_pca", args.order_pca, "PCA1", "PCA2", "Order + phylo PCA with metadata RDA arrows", "rda_order_phylo_pca"),
        ("bp_phylo_pca", args.bp_pca, "PCA1", "PCA2", "bp + phylo PCA with metadata RDA arrows", "rda_bp_phylo_pca"),
        ("order_phylo_pcoa", args.order_pcoa, "PCoA1", "PCoA2", "Order + phylo PCoA with metadata RDA arrows", "rda_order_phylo_pcoa"),
        ("bp_phylo_pcoa", args.bp_pcoa, "PCoA1", "PCoA2", "bp + phylo PCoA with metadata RDA arrows", "rda_bp_phylo_pcoa"),
    ]:
        vectors = process_plot(path, assoc, x_col, y_col, title, args.outdir / out_name)
        top = vectors.head(1).iloc[0] if not vectors.empty else None
        results.append(
            {
                "plot": label,
                "n_vectors": len(vectors),
                "top_metadata": top["metadata"] if top is not None else "",
                "top_r2_on_axes": top["r2_on_axes"] if top is not None else np.nan,
            }
        )

    summary = pd.DataFrame(results)
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
