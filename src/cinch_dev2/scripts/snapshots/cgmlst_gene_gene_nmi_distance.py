#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numba import njit, prange
from scipy.stats import chi2


def ref_num(gene: str) -> int:
    match = re.search(r"ref_gene_(\d+)$", str(gene))
    return int(match.group(1)) if match else -1


def normalize_gene(gene: str, prefix: str) -> str:
    n = ref_num(gene)
    return f"{prefix}:ref_gene_{n}" if n >= 0 else str(gene)


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


@njit
def entropy_from_counts(counts, n):
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            h -= p * np.log(p)
    return h


@njit
def pair_nmi_gtest(x, y, kx, ky):
    n = x.shape[0]
    cx = np.zeros(kx, dtype=np.int32)
    cy = np.zeros(ky, dtype=np.int32)
    cxy = np.zeros((kx, ky), dtype=np.int32)
    for i in range(n):
        a = x[i]
        b = y[i]
        cx[a] += 1
        cy[b] += 1
        cxy[a, b] += 1

    hx = entropy_from_counts(cx, n)
    hy = entropy_from_counts(cy, n)
    hxy = 0.0
    for a in range(kx):
        for b in range(ky):
            c = cxy[a, b]
            if c > 0:
                p = c / n
                hxy -= p * np.log(p)
    mi = max(0.0, hx + hy - hxy)
    g_stat = 2.0 * n * mi

    used_x = 0
    used_y = 0
    for a in range(kx):
        if cx[a] > 0:
            used_x += 1
    for b in range(ky):
        if cy[b] > 0:
            used_y += 1
    df = max(1, (used_x - 1) * (used_y - 1))
    if hx + hy <= 0:
        return 0.0, g_stat, df
    return mi / ((hx + hy) / 2.0), g_stat, df


@njit(parallel=True)
def all_pair_nmi_gtest(encoded, k_eff):
    n_genes = encoded.shape[1]
    n_pairs = n_genes * (n_genes - 1) // 2
    g1 = np.empty(n_pairs, dtype=np.int32)
    g2 = np.empty(n_pairs, dtype=np.int32)
    nmi = np.empty(n_pairs, dtype=np.float32)
    g_stat = np.empty(n_pairs, dtype=np.float32)
    dfs = np.empty(n_pairs, dtype=np.int32)

    idx = 0
    for i in range(n_genes - 1):
        for j in range(i + 1, n_genes):
            g1[idx] = i
            g2[idx] = j
            idx += 1

    for p in prange(n_pairs):
        i = g1[p]
        j = g2[p]
        score, g, df = pair_nmi_gtest(encoded[:, i], encoded[:, j], int(k_eff[i]), int(k_eff[j]))
        nmi[p] = score
        g_stat[p] = g
        dfs[p] = df
    return g1, g2, nmi, g_stat, dfs


def encode_core_alleles(profile: pd.DataFrame, core_min: float, allele_min_count: int) -> tuple[np.ndarray, list[str], np.ndarray, pd.DataFrame]:
    x = profile.to_numpy()
    prev = (x > 0).mean(axis=0)
    keep = prev >= core_min
    genes = profile.columns[keep].tolist()
    x = x[:, keep]
    encoded = np.zeros_like(x, dtype=np.int16)
    k_eff = np.zeros(x.shape[1], dtype=np.int16)
    rows = []

    for j, gene in enumerate(genes):
        vals = x[:, j]
        pos = vals[vals > 0]
        alleles, counts = np.unique(pos, return_counts=True)
        kept = alleles[counts >= allele_min_count]
        mapping = {int(a): i + 1 for i, a in enumerate(kept)}
        out = np.zeros(vals.shape[0], dtype=np.int16)
        for allele, code in mapping.items():
            out[vals == allele] = code
        encoded[:, j] = out
        k_eff[j] = len(kept) + 1
        rows.append(
            {
                "gene": gene,
                "prevalence": float((vals > 0).mean()),
                "raw_n_alleles": int(len(alleles)),
                "kept_alleles": int(len(kept)),
                "encoded_categories_including_rare_or_missing": int(k_eff[j]),
                "singletons": int((counts == 1).sum()),
            }
        )
    return encoded, genes, k_eff, pd.DataFrame(rows)


def read_exclude_genes(path: Path | None) -> set[str]:
    if path is None:
        return set()
    table = pd.read_csv(path, sep="\t")
    if "gene" in table.columns:
        return set(table["gene"].dropna().astype(str))
    if table.shape[1] == 0:
        return set()
    return set(table.iloc[:, 0].dropna().astype(str))


def read_positions(path: Path, gene_prefix: str) -> pd.DataFrame:
    pos = pd.read_csv(path, sep="\t")
    required = {"gene", "contig", "start", "end", "contig_length"}
    missing = required - set(pos.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
    pos = pos.copy()
    pos["gene"] = pos["gene"].map(lambda g: normalize_gene(g, gene_prefix))
    pos["start"] = pd.to_numeric(pos["start"], errors="coerce")
    pos["end"] = pd.to_numeric(pos["end"], errors="coerce")
    pos["contig_length"] = pd.to_numeric(pos["contig_length"], errors="coerce")
    pos = pos.dropna(subset=["start", "end", "contig_length"]).drop_duplicates("gene", keep="first")

    contigs = pos[["contig", "contig_length"]].drop_duplicates("contig", keep="first").reset_index(drop=True)
    contigs["offset"] = contigs["contig_length"].shift(fill_value=0).cumsum()
    pos = pos.merge(contigs[["contig", "offset"]], on="contig", how="left")
    pos["midpoint"] = (pos["start"] + pos["end"]) / 2.0
    pos["flat_midpoint"] = pos["offset"] + pos["midpoint"]
    pos["order_index"] = pos["flat_midpoint"].rank(method="first").astype(int) - 1
    return pos[["gene", "contig", "start", "end", "contig_length", "offset", "midpoint", "flat_midpoint", "order_index"]]


def add_distances(pairs: pd.DataFrame, pos: pd.DataFrame) -> pd.DataFrame:
    p1 = pos.add_prefix("gene1_").rename(columns={"gene1_gene": "gene1"})
    p2 = pos.add_prefix("gene2_").rename(columns={"gene2_gene": "gene2"})
    out = pairs.merge(p1, on="gene1", how="inner").merge(p2, on="gene2", how="inner")
    out["order_dist"] = (out["gene1_order_index"] - out["gene2_order_index"]).abs()
    out["bp_dist"] = (out["gene1_flat_midpoint"] - out["gene2_flat_midpoint"]).abs()
    out["same_contig"] = out["gene1_contig"].eq(out["gene2_contig"])
    return out


def smooth_series(y: pd.Series) -> pd.Series:
    window = min(9, max(3, len(y) // 6 * 2 + 1))
    return y.rolling(window=window, min_periods=1, center=True).mean()


def elbow_threshold(df: pd.DataFrame, x_col: str, y_col: str, out_base: Path, log_x: bool, bins: int, q: float, sample_n: int) -> tuple[float, pd.DataFrame]:
    work = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    work = work[work[x_col] >= 0].copy()
    if work[x_col].nunique() < 5:
        threshold = float(work[x_col].median()) if len(work) else 0.0
        return threshold, pd.DataFrame()
    n_bins = min(bins, max(8, work[x_col].nunique()))
    work["_bin"] = pd.qcut(work[x_col].rank(method="first"), q=n_bins, duplicates="drop")
    grouped = work.groupby("_bin", observed=True)
    binned = grouped.agg(x_mid=(x_col, "median"), y_median=(y_col, "median"), n=(y_col, "size")).reset_index(drop=True)
    binned["y_q"] = grouped[y_col].quantile(q).to_numpy()
    binned["y_q_smooth"] = smooth_series(binned["y_q"])
    binned["y_median_smooth"] = smooth_series(binned["y_median"])

    x = binned["x_mid"].to_numpy(dtype=float)
    y = binned["y_q_smooth"].to_numpy(dtype=float)
    tx = np.log10(x + 1.0) if log_x else x
    valid = np.isfinite(tx) & np.isfinite(y)
    if valid.sum() < 5 or tx[valid].max() == tx[valid].min() or y[valid].max() == y[valid].min():
        threshold = float(np.nanmedian(x))
    else:
        txv = tx[valid]
        yv = y[valid]
        x_norm = (txv - txv.min()) / (txv.max() - txv.min())
        y_norm = (yv - yv.min()) / (yv.max() - yv.min())
        distance = (1.0 - x_norm) - y_norm
        distance[:2] = -np.inf
        distance[-2:] = -np.inf
        idx = np.where(valid)[0][int(np.nanargmax(distance))]
        threshold = float(binned.loc[idx, "x_mid"])

    plot_distance_decay(work, binned, x_col, y_col, threshold, out_base, log_x, sample_n)
    return threshold, binned


def plot_distance_decay(work: pd.DataFrame, binned: pd.DataFrame, x_col: str, y_col: str, threshold: float, out_base: Path, log_x: bool, sample_n: int) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    sample = work.sample(min(len(work), sample_n), random_state=42) if len(work) > sample_n else work
    fig, ax = plt.subplots(figsize=(9.0, 6.4))
    ax.scatter(sample[x_col], sample[y_col], s=6, alpha=0.08, linewidths=0, color="#355c7d", rasterized=True)
    if not binned.empty:
        ax.plot(binned["x_mid"], binned["y_q_smooth"], color="#c23b22", lw=2.3, label="smoothed high-quantile NMI")
        ax.plot(binned["x_mid"], binned["y_median_smooth"], color="#2a9d8f", lw=1.8, alpha=0.8, label="smoothed median NMI")
    ax.axvline(threshold, color="#202020", ls="--", lw=1.5, label=f"threshold={threshold:.4g}")
    if log_x:
        ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.legend(frameon=False, fontsize=9)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_top_links(top: pd.DataFrame, pos: pd.DataFrame, out_base: Path, top_n: int) -> None:
    show = top.head(top_n).copy()
    if show.empty:
        return
    contigs = pos[["contig", "contig_length", "offset"]].drop_duplicates("contig", keep="first").sort_values("offset")
    genome_len = float((contigs["offset"] + contigs["contig_length"]).max())
    y_base = 0.0

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(14, 5.2))
    for i, row in enumerate(contigs.itertuples(index=False)):
        color = "#e8e8e8" if i % 2 == 0 else "#d7d7d7"
        ax.add_patch(plt.Rectangle((row.offset, y_base - 0.04), row.contig_length, 0.08, color=color, ec="#777777", lw=0.3))
        ax.text(row.offset + row.contig_length / 2, y_base - 0.11, row.contig, ha="center", va="top", fontsize=8, rotation=30)

    norm = plt.Normalize(show["nmi"].min(), show["nmi"].max())
    cmap = plt.get_cmap("magma")
    for rank, row in enumerate(show.itertuples(index=False), start=1):
        x1 = float(row.gene1_flat_midpoint)
        x2 = float(row.gene2_flat_midpoint)
        left, right = sorted((x1, x2))
        height = 0.18 + 0.72 * (top_n + 1 - rank) / top_n
        xs = np.linspace(left, right, 120)
        if right == left:
            ys = np.full_like(xs, y_base)
        else:
            ys = y_base + height * np.sin(np.pi * (xs - left) / (right - left))
        color = cmap(norm(float(row.nmi)))
        ax.plot(xs, ys, color=color, lw=1.4 + 1.0 * (top_n + 1 - rank) / top_n, alpha=0.82)
        ax.scatter([x1, x2], [y_base, y_base], s=26, color=color, zorder=5, edgecolor="white", linewidth=0.4)

    ax.set_xlim(-genome_len * 0.01, genome_len * 1.01)
    ax.set_ylim(-0.28, 1.15)
    ax.set_yticks([])
    ax.set_xlabel("flattened genome coordinate (bp)")
    ax.set_title(f"Top {len(show)} distance-filtered cgMLST gene-gene NMI links")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label("NMI")
    sns.despine(fig=fig, ax=ax, left=True)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="cgMLST gene-gene NMI, approximate chi-square, distance filtering, and genome-link plotting.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument("--positions", type=Path, default=Path("Leca/source/positions/test0427_1_GCA_905232925.fna.positions.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/results/gene_gene_cgmlst"))
    parser.add_argument("--gene-prefix", default="0427_1")
    parser.add_argument("--core-min", type=float, default=0.95)
    parser.add_argument("--allele-min-count", type=int, default=20)
    parser.add_argument("--p-adj-cutoff", type=float, default=0.05)
    parser.add_argument("--order-threshold", type=float, default=None)
    parser.add_argument("--bp-threshold", type=float, default=None)
    parser.add_argument("--exclude-genes", type=Path, default=None)
    parser.add_argument("--elbow-bins", type=int, default=120)
    parser.add_argument("--elbow-quantile", type=float, default=0.90)
    parser.add_argument("--plot-sample-n", type=int, default=200000)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--write-all-pairs", action="store_true")
    parser.add_argument("--write-significant-pairs", action="store_true")
    parser.add_argument("--write-filtered-pairs", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    profile = pd.read_csv(args.profile, sep="\t", index_col=0).fillna(0)
    encoded, genes, k_eff, allele_summary = encode_core_alleles(profile, args.core_min, args.allele_min_count)
    exclude_genes = read_exclude_genes(args.exclude_genes)
    n_core_genes_before_exclusion = len(genes)
    excluded_core_genes = sorted(set(genes) & exclude_genes, key=ref_num)
    if excluded_core_genes:
        keep = np.array([gene not in exclude_genes for gene in genes], dtype=bool)
        encoded = encoded[:, keep]
        k_eff = k_eff[keep]
        genes = [gene for gene in genes if gene not in exclude_genes]
        allele_summary = allele_summary[~allele_summary["gene"].isin(exclude_genes)].copy()
        pd.DataFrame({"gene": excluded_core_genes}).to_csv(args.outdir / "excluded_genes_used.tsv", sep="\t", index=False)
    allele_summary.to_csv(args.outdir / "core_gene_allele_summary.tsv", sep="\t", index=False)

    positions = read_positions(args.positions, args.gene_prefix)
    positions.to_csv(args.outdir / "flattened_gene_positions.tsv", sep="\t", index=False)

    g1, g2, scores, g_stats, dfs = all_pair_nmi_gtest(encoded, k_eff)
    p_values = chi2.sf(g_stats.astype(np.float64), dfs)
    q_values = bh_qvalues(p_values)
    gene_array = np.asarray(genes, dtype=object)
    pairs = pd.DataFrame(
        {
            "gene1_idx": g1,
            "gene2_idx": g2,
            "gene1": gene_array[g1],
            "gene2": gene_array[g2],
            "nmi": scores,
            "approx_chisq": g_stats,
            "approx_df": dfs,
            "approx_p": p_values,
            "approx_q": q_values,
            "gene1_kept_alleles": k_eff[g1] - 1,
            "gene2_kept_alleles": k_eff[g2] - 1,
        }
    )
    pairs = add_distances(pairs, positions)
    if args.write_all_pairs:
        pairs.to_csv(args.outdir / "all_gene_gene_nmi_chisq_distances.tsv.gz", sep="\t", index=False)

    sig = pairs[pairs["approx_q"] <= args.p_adj_cutoff].copy()
    if args.write_significant_pairs:
        sig.to_csv(args.outdir / "significant_gene_gene_pairs.tsv.gz", sep="\t", index=False)
    if sig.empty:
        raise ValueError("No significant pairs after approximate q-value filtering")

    order_threshold = args.order_threshold
    bp_threshold = args.bp_threshold
    binned_tables = []
    if order_threshold is None:
        order_threshold, binned = elbow_threshold(
            sig,
            "order_dist",
            "nmi",
            args.outdir / "decay_nmi_vs_order_dist",
            False,
            args.elbow_bins,
            args.elbow_quantile,
            args.plot_sample_n,
        )
        if not binned.empty:
            binned.insert(0, "metric", "order_dist")
            binned_tables.append(binned)
    if bp_threshold is None:
        bp_threshold, binned = elbow_threshold(
            sig,
            "bp_dist",
            "nmi",
            args.outdir / "decay_nmi_vs_bp_dist",
            True,
            args.elbow_bins,
            args.elbow_quantile,
            args.plot_sample_n,
        )
        if not binned.empty:
            binned.insert(0, "metric", "bp_dist")
            binned_tables.append(binned)
    if binned_tables:
        pd.concat(binned_tables, ignore_index=True).to_csv(args.outdir / "distance_decay_binned_trends.tsv", sep="\t", index=False)

    thresholds = pd.DataFrame(
        [
            {"metric": "order_dist", "threshold": order_threshold, "direction": "keep >= threshold"},
            {"metric": "bp_dist", "threshold": bp_threshold, "direction": "keep >= threshold"},
        ]
    )
    thresholds.to_csv(args.outdir / "distance_filter_thresholds.tsv", sep="\t", index=False)

    filtered = sig[(sig["order_dist"] >= order_threshold) & (sig["bp_dist"] >= bp_threshold)].copy()
    filtered = filtered.sort_values(["approx_q", "approx_p", "nmi"], ascending=[True, True, False])
    if args.write_filtered_pairs:
        filtered.to_csv(args.outdir / "distance_filtered_gene_gene_pairs.tsv.gz", sep="\t", index=False)
    top = filtered.head(args.top_n).copy()
    top.to_csv(args.outdir / f"top{args.top_n}_distance_filtered_pairs.tsv", sep="\t", index=False)
    plot_top_links(top, positions, args.outdir / f"top{args.top_n}_flattened_genome_links", args.top_n)

    summary = pd.DataFrame(
        [
            {
                "n_genomes": int(encoded.shape[0]),
                "n_core_genes": int(encoded.shape[1]),
                "n_core_genes_before_exclusion": int(n_core_genes_before_exclusion),
                "n_excluded_core_genes": int(len(excluded_core_genes)),
                "gene_pairs_total": int(len(pairs)),
                "gene_pairs_with_positions": int(len(pairs)),
                "significant_pairs_q_le_cutoff": int(len(sig)),
                "p_adj_cutoff": args.p_adj_cutoff,
                "order_dist_threshold": float(order_threshold),
                "bp_dist_threshold": float(bp_threshold),
                "distance_filtered_pairs": int(len(filtered)),
                "top_links_plotted": int(len(top)),
                "min_filtered_q": float(filtered["approx_q"].min()) if len(filtered) else np.nan,
                "max_top20_nmi": float(top["nmi"].max()) if len(top) else np.nan,
            }
        ]
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("\nDistance thresholds:")
    print(thresholds.to_string(index=False))


if __name__ == "__main__":
    main()
