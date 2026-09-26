#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler


UNKNOWN_VALUES = {
    "",
    "unknown",
    "not_available",
    "not available",
    "missing",
    "na",
    "nan",
    "none",
    "null",
    "unspecified",
}

DEFAULT_DROP_METADATA = {
    "Query",
    "Sample_Name",
    "Other_Notes",
    "SRA_Accession",
    "Assembly",
    "EBI_Accession",
    "Biosample",
    "Description",
    "Coordinates",
}

CONTINUOUS_METADATA = {"Isolate_Date"}

RELEVANT_METADATA = {
    "Isolation_Source",
    "Host",
    "Geography",
    "Region",
    "Isolate_Date",
}


def canonical_genome_name(name: str) -> str:
    match = re.search(r"(GC[AF]_\d+)", str(name))
    return match.group(1) if match else str(name)


def short_gene_name(name: str) -> str:
    return str(name).split(":")[-1]


def clean_value(value: object) -> object:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if text.lower() in UNKNOWN_VALUES:
        return np.nan
    return text


def classical_pcoa(distance_matrix: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    n = distance_matrix.shape[0]
    d2 = np.square(distance_matrix)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 0
    coords = eigvecs[:, :n_components] * np.sqrt(np.maximum(eigvals[:n_components], 0))
    denom = eigvals[positive].sum()
    var_ratio = eigvals[:n_components] / denom if denom > 0 else np.full(n_components, np.nan)
    return coords, var_ratio


def make_pcoa_from_matrix(mat: pd.DataFrame, metric: str = "cosine") -> tuple[pd.DataFrame, np.ndarray]:
    distances = pairwise_distances(mat.values, metric=metric)
    distances = np.nan_to_num(distances, nan=0.0, posinf=1.0, neginf=0.0)
    coords, var_ratio = classical_pcoa(distances, 2)
    degree = (mat.values > 0).sum(axis=1)
    strength = mat.values.sum(axis=1)
    out = pd.DataFrame(
        {
            "gene": mat.index,
            "gene_short": [short_gene_name(g) for g in mat.index],
            "PCoA1": coords[:, 0],
            "PCoA2": coords[:, 1],
            "degree": degree,
            "r2_strength": strength,
        }
    )
    return out, var_ratio


def plot_ordination(
    coords: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    size_col: str,
    title: str,
    out_base: Path,
    var_ratio: np.ndarray | None = None,
    label_top: int = 12,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.6, 7.2))
    points = ax.scatter(
        coords[x_col],
        coords[y_col],
        c=coords[color_col],
        s=np.clip(22 + coords[size_col] * 1.8, 26, 160),
        cmap="viridis",
        alpha=0.84,
        edgecolors="white",
        linewidths=0.25,
    )
    cbar = fig.colorbar(points, ax=ax)
    cbar.set_label(color_col)
    for row in coords.sort_values([size_col, color_col], ascending=False).head(label_top).itertuples(index=False):
        ax.text(getattr(row, x_col), getattr(row, y_col), row.gene_short, fontsize=8, ha="left", va="bottom", alpha=0.9)
    ax.set_title(title)
    if var_ratio is not None and np.all(np.isfinite(var_ratio)):
        ax.set_xlabel(f"{x_col} ({var_ratio[0] * 100:.1f}%)")
        ax.set_ylabel(f"{y_col} ({var_ratio[1] * 100:.1f}%)")
    else:
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip", index_col=False)
    meta["genome_key"] = meta["Query"].map(canonical_genome_name)
    meta = meta[meta["genome_key"].str.match(r"GC[AF]_\d+", na=False)].copy()
    meta = meta.drop_duplicates("genome_key", keep="first")
    for col in meta.columns:
        if col != "genome_key":
            meta[col] = meta[col].map(clean_value)
    return meta


def select_metadata_columns(
    meta: pd.DataFrame,
    candidate_columns: list[str],
    min_nonmissing: int,
    min_level_count: int,
    max_categories: int,
) -> tuple[list[str], pd.DataFrame]:
    rows = []
    selected = []
    for col in candidate_columns:
        if col in DEFAULT_DROP_METADATA or col == "genome_key" or col not in RELEVANT_METADATA:
            continue
        s = meta[col].dropna().astype(str)
        counts = s.value_counts()
        counts = counts[counts >= min_level_count]
        n_nonmissing = int(s.size)
        n_categories = int(counts.size)
        keep = n_nonmissing >= min_nonmissing and 2 <= n_categories <= max_categories
        rows.append(
            {
                "metadata": col,
                "nonmissing": n_nonmissing,
                "categories_after_min_count": n_categories,
                "top_value": counts.index[0] if len(counts) else "",
                "top_count": int(counts.iloc[0]) if len(counts) else 0,
                "selected": keep,
            }
        )
        if keep:
            selected.append(col)
    return selected, pd.DataFrame(rows)


def binary_by_categorical_nmi_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    if n == 0:
        return np.zeros(x.shape[1], dtype=np.float32)

    levels = np.unique(y)
    n1 = x.sum(axis=0).astype(np.float64)
    n0 = n - n1
    px1 = n1 / n
    px0 = n0 / n
    hx = -(px1 * np.log(px1 + 1e-12) + px0 * np.log(px0 + 1e-12))

    level_counts = np.array([(y == level).sum() for level in levels], dtype=np.float64)
    py = level_counts / n
    hy = -np.sum(py * np.log(py + 1e-12))
    mi = np.zeros(x.shape[1], dtype=np.float64)

    for level, cy in zip(levels, level_counts):
        mask = y == level
        c1y = x[mask].sum(axis=0).astype(np.float64)
        c0y = cy - c1y
        p1y = c1y / n
        p0y = c0y / n
        py_scalar = cy / n

        term1 = p1y * np.log((p1y + 1e-12) / (px1 * py_scalar + 1e-12))
        term0 = p0y * np.log((p0y + 1e-12) / (px0 * py_scalar + 1e-12))
        term1[p1y <= 0] = 0.0
        term0[p0y <= 0] = 0.0
        mi += term1 + term0

    denom = (hx + hy) / 2.0
    nmi = np.where(denom > 0, mi / denom, 0.0)
    return np.nan_to_num(nmi, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def binary_by_binary_nmi_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return binary_by_categorical_nmi_matrix(x, y.astype(np.int8))


def binary_by_continuous_corr_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    y = y.astype(float)
    y_centered = y - np.nanmean(y)
    y_sd = np.sqrt(np.nanmean(y_centered * y_centered))
    if not np.isfinite(y_sd) or y_sd == 0:
        return np.zeros(x.shape[1], dtype=np.float32)

    x_mean = x.mean(axis=0)
    x_centered = x - x_mean
    x_sd = np.sqrt(np.mean(x_centered * x_centered, axis=0))
    cov = np.mean(x_centered * y_centered[:, None], axis=0)
    corr = np.divide(cov, x_sd * y_sd, out=np.zeros_like(cov, dtype=float), where=x_sd > 0)
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def extract_year(series: pd.Series) -> pd.Series:
    years = series.astype(str).str.extract(r"((?:19|20)\d{2})")[0]
    return pd.to_numeric(years, errors="coerce")


def build_gene_metadata_nmi(
    profile_path: Path,
    metadata_path: Path,
    genes: list[str],
    min_nonmissing: int,
    min_level_count: int,
    max_categories: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    profile = pd.read_csv(profile_path, sep="\t")
    profile["genome_key"] = profile["genome"].map(canonical_genome_name)
    genes = [g for g in genes if g in profile.columns]
    profile = profile[["genome_key", *genes]].copy()
    profile[genes] = (profile[genes] > 0).astype(np.int8)

    meta = load_metadata(metadata_path)
    candidate_metadata_cols = [c for c in meta.columns if c != "genome_key"]
    merged = profile.merge(meta, on="genome_key", how="inner")
    selected_cols, meta_summary = select_metadata_columns(
        merged,
        candidate_metadata_cols,
        min_nonmissing,
        min_level_count,
        max_categories,
    )
    if not selected_cols:
        raise ValueError("No metadata columns passed filtering. Relax min counts or max categories.")

    x = merged[genes].to_numpy(dtype=np.int8)
    assoc_columns = []
    assoc_arrays = []
    assoc_meta = []

    for col in selected_cols:
        if col in CONTINUOUS_METADATA:
            y_num = extract_year(merged[col])
            mask = y_num.notna().to_numpy()
            if mask.sum() >= min_nonmissing:
                vals = binary_by_continuous_corr_matrix(x[mask, :], y_num[mask].to_numpy(dtype=float))
                assoc_columns.append(f"{col}:year")
                assoc_arrays.append(vals)
                assoc_meta.append(
                    {
                        "association_variable": f"{col}:year",
                        "source_metadata": col,
                        "encoding": "continuous_year",
                        "level": "",
                        "n_genomes": int(mask.sum()),
                        "score": "pearson_r_gene_presence_vs_year",
                    }
                )
            continue

        counts = merged[col].dropna().astype(str).value_counts()
        levels = counts[counts >= min_level_count].index.tolist()
        for level in levels:
            mask = merged[col].notna().to_numpy()
            y = (merged.loc[mask, col].astype(str).to_numpy() == level).astype(np.int8)
            if y.min() == y.max():
                vals = np.zeros(len(genes), dtype=np.float32)
            else:
                vals = binary_by_binary_nmi_matrix(x[mask, :], y)
            safe_level = str(level).replace("\t", " ").replace("\n", " ")
            assoc_columns.append(f"{col}={safe_level}")
            assoc_arrays.append(vals)
            assoc_meta.append(
                {
                    "association_variable": f"{col}={safe_level}",
                    "source_metadata": col,
                    "encoding": "one_vs_rest_binary",
                    "level": safe_level,
                    "n_genomes": int(mask.sum()),
                    "level_count": int((y == 1).sum()),
                    "score": "nmi_gene_presence_vs_level_presence",
                }
            )

    if not assoc_arrays:
        raise ValueError("No metadata association variables could be constructed.")

    assoc = np.vstack(assoc_arrays).T.astype(np.float32)
    assoc_df = pd.DataFrame(assoc, index=genes, columns=assoc_columns)
    assoc_df.index.name = "gene"
    assoc_var_summary = pd.DataFrame(assoc_meta)

    sample_summary = pd.DataFrame(
        {
            "aligned_genomes": [len(merged)],
            "genes_used": [len(genes)],
            "metadata_variables_used": [len(selected_cols)],
            "association_variables_used": [len(assoc_columns)],
        }
    )
    return assoc_df, meta_summary, sample_summary, assoc_var_summary


def make_pca_from_assoc(assoc: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    x = StandardScaler().fit_transform(assoc.values)
    pca = PCA(n_components=2, random_state=42)
    xy = pca.fit_transform(x)
    max_nmi = assoc.max(axis=1).to_numpy()
    strength = assoc.sum(axis=1).to_numpy()
    top_var = assoc.idxmax(axis=1).to_numpy()
    coords = pd.DataFrame(
        {
            "gene": assoc.index,
            "gene_short": [short_gene_name(g) for g in assoc.index],
            "PCA1": xy[:, 0],
            "PCA2": xy[:, 1],
            "max_metadata_nmi": max_nmi,
            "metadata_assoc_strength": strength,
            "top_metadata": top_var,
        }
    )
    return coords, pca.explained_variance_ratio_


def make_pcoa_from_assoc(assoc: pd.DataFrame, metric: str) -> tuple[pd.DataFrame, np.ndarray]:
    distances = pairwise_distances(assoc.values, metric=metric)
    distances = np.nan_to_num(distances, nan=0.0, posinf=1.0, neginf=0.0)
    xy, var_ratio = classical_pcoa(distances, 2)
    coords = pd.DataFrame(
        {
            "gene": assoc.index,
            "gene_short": [short_gene_name(g) for g in assoc.index],
            "PCoA1": xy[:, 0],
            "PCoA2": xy[:, 1],
            "max_metadata_nmi": assoc.max(axis=1).to_numpy(),
            "metadata_assoc_strength": assoc.sum(axis=1).to_numpy(),
            "top_metadata": assoc.idxmax(axis=1).to_numpy(),
        }
    )
    return coords, var_ratio


def main() -> None:
    parser = argparse.ArgumentParser(description="Add PCoA and gene-metadata association ordinations.")
    parser.add_argument("--order-matrix", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/order_phylo_gene_by_gene_r2_matrix.tsv"))
    parser.add_argument("--bp-matrix", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/bp_phylo_gene_by_gene_r2_matrix.tsv"))
    parser.add_argument("--profile", type=Path, default=Path("Leca/0427_1.profile"))
    parser.add_argument("--metadata", type=Path, default=Path("Leca/Kp.mp.merge"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/pcoa_metadata_0427_1"))
    parser.add_argument("--distance-metric", default="cosine")
    parser.add_argument("--min-nonmissing", type=int, default=40)
    parser.add_argument("--min-level-count", type=int, default=8)
    parser.add_argument("--max-categories", type=int, default=40)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    order_mat = pd.read_csv(args.order_matrix, sep="\t", index_col=0)
    bp_mat = pd.read_csv(args.bp_matrix, sep="\t", index_col=0)
    genes = sorted(set(order_mat.index).union(bp_mat.index))

    order_pcoa, order_var = make_pcoa_from_matrix(order_mat, args.distance_metric)
    bp_pcoa, bp_var = make_pcoa_from_matrix(bp_mat, args.distance_metric)
    order_pcoa.to_csv(args.outdir / "order_phylo_r2_profile_pcoa_coords.tsv", sep="\t", index=False)
    bp_pcoa.to_csv(args.outdir / "bp_phylo_r2_profile_pcoa_coords.tsv", sep="\t", index=False)
    plot_ordination(order_pcoa, "PCoA1", "PCoA2", "r2_strength", "degree", "Order + phylo filtered r2 profile PCoA", args.outdir / "pcoa_order_phylo_r2_profile", order_var)
    plot_ordination(bp_pcoa, "PCoA1", "PCoA2", "r2_strength", "degree", "bp + phylo filtered r2 profile PCoA", args.outdir / "pcoa_bp_phylo_r2_profile", bp_var)

    assoc, meta_summary, sample_summary, assoc_var_summary = build_gene_metadata_nmi(
        args.profile,
        args.metadata,
        genes,
        args.min_nonmissing,
        args.min_level_count,
        args.max_categories,
    )
    assoc.to_csv(args.outdir / "gene_by_metadata_nmi_matrix.tsv", sep="\t")
    assoc.to_csv(args.outdir / "gene_by_metadata_expanded_assoc_matrix.tsv", sep="\t")
    meta_summary.to_csv(args.outdir / "metadata_variable_filter_summary.tsv", sep="\t", index=False)
    assoc_var_summary.to_csv(args.outdir / "metadata_association_variable_summary.tsv", sep="\t", index=False)
    sample_summary.to_csv(args.outdir / "gene_metadata_sample_summary.tsv", sep="\t", index=False)

    pca_coords, pca_var = make_pca_from_assoc(assoc)
    pcoa_coords, pcoa_var = make_pcoa_from_assoc(assoc, args.distance_metric)
    pca_coords.to_csv(args.outdir / "gene_metadata_nmi_pca_coords.tsv", sep="\t", index=False)
    pcoa_coords.to_csv(args.outdir / "gene_metadata_nmi_pcoa_coords.tsv", sep="\t", index=False)

    plot_ordination(
        pca_coords,
        "PCA1",
        "PCA2",
        "max_metadata_nmi",
        "metadata_assoc_strength",
        "Gene x metadata NMI PCA",
        args.outdir / "pca_gene_metadata_nmi",
        pca_var,
    )
    plot_ordination(
        pcoa_coords,
        "PCoA1",
        "PCoA2",
        "max_metadata_nmi",
        "metadata_assoc_strength",
        "Gene x metadata NMI PCoA",
        args.outdir / "pcoa_gene_metadata_nmi",
        pcoa_var,
    )

    top_assoc = (
        assoc.stack()
        .rename("nmi")
        .reset_index()
        .rename(columns={"level_1": "metadata"})
        .sort_values("nmi", ascending=False)
        .head(5000)
    )
    top_assoc.to_csv(args.outdir / "top_gene_metadata_nmi_pairs.tsv", sep="\t", index=False)

    summary = pd.Series(
        {
            "order_pcoa1_var": order_var[0],
            "order_pcoa2_var": order_var[1],
            "bp_pcoa1_var": bp_var[0],
            "bp_pcoa2_var": bp_var[1],
            "gene_metadata_pca1_var": pca_var[0],
            "gene_metadata_pca2_var": pca_var[1],
            "gene_metadata_pcoa1_var": pcoa_var[0],
            "gene_metadata_pcoa2_var": pcoa_var[1],
            "genes_used": assoc.shape[0],
            "association_variables_used": assoc.shape[1],
            "metadata_variables_used": int(sample_summary["metadata_variables_used"].iloc[0]),
            "aligned_genomes": int(sample_summary["aligned_genomes"].iloc[0]),
            "distance_metric": args.distance_metric,
        }
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", header=False)
    print(summary.to_string())
    print("\nConstructed association variables:")
    print(", ".join(assoc.columns[:80]) + (" ..." if len(assoc.columns) > 80 else ""))


if __name__ == "__main__":
    main()
