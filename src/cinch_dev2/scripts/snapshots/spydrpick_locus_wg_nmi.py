#!/usr/bin/env python3
import argparse
import re
import textwrap
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

from secondary_filter_ld_decay import find_elbow, make_quantile_bins, smooth_series


TOKENS = {
    "surface": "#FCFCFD",
    "panel": "#FFFFFF",
    "ink": "#1F2430",
    "muted": "#6F768A",
    "grid": "#E6E8F0",
    "axis": "#D7DBE7",
    "blue": "#5477C4",
    "blue_light": "#A3BEFA",
    "grey": "#C5CAD3",
    "grey_dark": "#464C55",
    "orange": "#CC6F47",
}


def binary_weighted_information(
    c11: np.ndarray,
    m1: np.ndarray,
    m2: np.ndarray,
    total_weight: float,
    pseudocount: float,
) -> tuple[np.ndarray, np.ndarray]:
    counts = np.column_stack(
        [
            total_weight - m1 - m2 + c11,
            m2 - c11,
            m1 - c11,
            c11,
        ]
    )
    counts = np.maximum(counts, 0.0) + pseudocount
    probabilities = counts / (total_weight + 4.0 * pseudocount)
    px0 = probabilities[:, 0] + probabilities[:, 1]
    px1 = probabilities[:, 2] + probabilities[:, 3]
    py0 = probabilities[:, 0] + probabilities[:, 2]
    py1 = probabilities[:, 1] + probabilities[:, 3]
    expected = np.column_stack([px0 * py0, px0 * py1, px1 * py0, px1 * py1])
    mi = np.sum(probabilities * np.log(probabilities / expected), axis=1)
    hx = -(px0 * np.log(px0) + px1 * np.log(px1))
    hy = -(py0 * np.log(py0) + py1 * np.log(py1))
    nmi = np.divide(2.0 * mi, hx + hy, out=np.zeros_like(mi), where=(hx + hy) > 0)
    return mi, nmi


def aracne_prune(edges: pd.DataFrame, weight_col: str) -> pd.DataFrame:
    edge_weight = {}
    neighbours = {}
    for row in edges[["gene1", "gene2", weight_col]].itertuples(index=False):
        key = tuple(sorted((row.gene1, row.gene2)))
        edge_weight[key] = float(getattr(row, weight_col))
        neighbours.setdefault(row.gene1, set()).add(row.gene2)
        neighbours.setdefault(row.gene2, set()).add(row.gene1)

    indirect = set()
    for center, adjacent in neighbours.items():
        adjacent = sorted(adjacent)
        for left, right in combinations(adjacent, 2):
            closing = tuple(sorted((left, right)))
            if closing not in edge_weight:
                continue
            triangle = [
                (tuple(sorted((center, left))), edge_weight[tuple(sorted((center, left)))]),
                (tuple(sorted((center, right))), edge_weight[tuple(sorted((center, right)))]),
                (closing, edge_weight[closing]),
            ]
            values = np.array([value for _, value in triangle])
            minimum = values.min()
            if np.sum(np.isclose(values, minimum, rtol=1e-12, atol=1e-12)) == 1:
                indirect.add(triangle[int(np.argmin(values))][0])
    out = edges.copy()
    out["aracne_direct"] = [
        tuple(sorted((gene1, gene2))) not in indirect
        for gene1, gene2 in zip(out["gene1"], out["gene2"])
    ]
    return out


def read_annotations(path: Path, genes: set[str]) -> dict[str, str]:
    annotations = {}
    pattern = re.compile(r"ref_gene_(\d+)")
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            match = pattern.search(fields[0])
            if not match:
                continue
            gene = f"0427_1:ref_gene_{match.group(1)}"
            if gene not in genes:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            annotations[gene] = attrs.get("gene", attrs.get("Name", gene.split(":")[-1]))
    return annotations


def add_header(fig, ax, title: str, subtitle: str) -> None:
    fig.subplots_adjust(top=0.80, bottom=0.12, left=0.07, right=0.98, wspace=0.17)
    left = ax.get_position().x0
    fig.text(left, 0.975, textwrap.fill(title, 90), ha="left", va="top", fontsize=15, fontweight="bold", color=TOKENS["ink"])
    fig.text(left, 0.92, textwrap.fill(subtitle, 150), ha="left", va="top", fontsize=9, color=TOKENS["muted"])


def plot_panel(
    ax,
    background: pd.DataFrame,
    outliers: pd.DataFrame,
    ld_threshold: float,
    outlier_threshold: float,
    extreme_threshold: float,
    x_limit: float | None,
    title: str,
) -> None:
    sns.scatterplot(
        data=background,
        x="mean_order_dist",
        y="weighted_nmi",
        color=TOKENS["grey"],
        edgecolor="none",
        s=7,
        alpha=0.28,
        rasterized=True,
        ax=ax,
    )
    indirect = outliers[~outliers["aracne_direct"]]
    direct = outliers[outliers["aracne_direct"]]
    if len(indirect):
        sns.scatterplot(
            data=indirect,
            x="mean_order_dist",
            y="weighted_nmi",
            color=TOKENS["grey_dark"],
            edgecolor="none",
            s=13,
            alpha=0.65,
            label="Indirect outlier",
            ax=ax,
        )
    if len(direct):
        sns.scatterplot(
            data=direct,
            x="mean_order_dist",
            y="weighted_nmi",
            color=TOKENS["blue"],
            edgecolor=TOKENS["blue"],
            s=16,
            alpha=0.85,
            label="Direct outlier",
            ax=ax,
        )
    ax.axvline(ld_threshold, color=TOKENS["orange"], linestyle=":", linewidth=1.2, label=f"LD elbow={ld_threshold:.1f}")
    ax.axhline(outlier_threshold, color=TOKENS["orange"], linestyle=":", linewidth=1.2, label=f"Outlier={outlier_threshold:.3f}")
    ax.axhline(extreme_threshold, color=TOKENS["orange"], linestyle="--", linewidth=1.2, label=f"Extreme={extreme_threshold:.3f}")
    if x_limit is not None:
        ax.set_xlim(0, x_limit)
    ax.set_xlabel("Mean gene order distance")
    ax.set_ylabel("Weighted binary NMI")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.grid(axis="y", color=TOKENS["grid"], linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    sns.despine(ax=ax)


def main() -> None:
    parser = argparse.ArgumentParser(description="SpydrPick-like weighted wgMLST locus NMI and ARACNE analysis.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument("--pairs", type=Path, default=Path("Leca/source/profile/0427_1.ldlike_pairs.tsv.gz"))
    parser.add_argument(
        "--weights",
        type=Path,
        default=Path("Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv"),
    )
    parser.add_argument("--gff", type=Path, default=Path("Leca/source/annotation/klebsiella_reference.gff"))
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi"),
    )
    parser.add_argument("--pseudocount", type=float, default=0.5)
    parser.add_argument("--bins", type=int, default=120)
    parser.add_argument("--trend-quantile", type=float, default=0.90)
    parser.add_argument("--plot-sample", type=int, default=250_000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    usecols = [
        "gene1",
        "gene2",
        "r2_adj",
        "p_adj",
        "n_joint",
        "mean_order_dist",
        "mean_bp_dist",
        "same_replicon_rate",
    ]
    pairs = pd.read_csv(args.pairs, sep="\t", compression="gzip", usecols=usecols)
    genes = sorted(set(pairs["gene1"]) | set(pairs["gene2"]))
    profile_header = pd.read_csv(args.profile, sep="\t", nrows=0).columns
    profile = pd.read_csv(args.profile, sep="\t", usecols=[profile_header[0], *genes], index_col=0).fillna(0)
    binary = (profile[genes].to_numpy() > 0).astype(np.float64)

    weight_table = pd.read_csv(args.weights, sep="\t").set_index("genome")
    weights = weight_table.loc[profile.index, "wg_content_weight"].to_numpy(dtype=np.float64)
    total_weight = float(weights.sum())
    weighted_marginal = weights @ binary
    weighted_joint = (binary * weights[:, None]).T @ binary
    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    idx1 = pairs["gene1"].map(gene_index).to_numpy()
    idx2 = pairs["gene2"].map(gene_index).to_numpy()
    mi, nmi = binary_weighted_information(
        weighted_joint[idx1, idx2],
        weighted_marginal[idx1],
        weighted_marginal[idx2],
        total_weight,
        args.pseudocount,
    )
    pairs["weighted_mi"] = mi
    pairs["weighted_nmi"] = nmi

    binned = make_quantile_bins(pairs, "mean_order_dist", "weighted_nmi", args.bins, args.trend_quantile)
    binned["y_median_smooth"] = smooth_series(binned["y_median"])
    binned["y_q_smooth"] = smooth_series(binned["y_q"])
    ld_threshold, binned = find_elbow(binned, log_x=False)
    binned.to_csv(args.outdir / "weighted_nmi_order_distance_binned.tsv", sep="\t", index=False)

    non_ld = pairs[pairs["mean_order_dist"] >= ld_threshold]
    maxima = pd.concat(
        [
            non_ld[["gene1", "weighted_nmi"]].rename(columns={"gene1": "gene"}),
            non_ld[["gene2", "weighted_nmi"]].rename(columns={"gene2": "gene"}),
        ],
        ignore_index=True,
    ).groupby("gene", as_index=False)["weighted_nmi"].max()
    q1, q3 = maxima["weighted_nmi"].quantile([0.25, 0.75])
    iqr = q3 - q1
    outlier_threshold = float(q3 + 1.5 * iqr)
    extreme_threshold = float(q3 + 3.0 * iqr)
    maxima.to_csv(args.outdir / "non_ld_per_locus_maximum_nmi.tsv", sep="\t", index=False)

    outliers = non_ld[non_ld["weighted_nmi"] > outlier_threshold].copy()
    outliers = aracne_prune(outliers, "weighted_nmi")
    outliers["extreme_outlier"] = outliers["weighted_nmi"] > extreme_threshold
    outliers = outliers.sort_values(["aracne_direct", "weighted_nmi"], ascending=[False, False])
    annotations = read_annotations(args.gff, set(genes))
    outliers["gene1_annotation"] = outliers["gene1"].map(annotations).fillna("")
    outliers["gene2_annotation"] = outliers["gene2"].map(annotations).fillna("")
    outliers.to_csv(args.outdir / "wgmlst_weighted_nmi_outliers_aracne.tsv", sep="\t", index=False)
    pairs.to_csv(args.outdir / "wgmlst_all_weighted_mi_nmi_pairs.tsv.gz", sep="\t", index=False)

    r2_binned = make_quantile_bins(pairs, "mean_order_dist", "r2_adj", args.bins, args.trend_quantile)
    r2_binned["y_median_smooth"] = smooth_series(r2_binned["y_median"])
    r2_binned["y_q_smooth"] = smooth_series(r2_binned["y_q"])
    r2_ld_threshold, r2_binned = find_elbow(r2_binned, log_x=False)
    r2_binned.to_csv(args.outdir / "r2_adj_order_distance_binned.tsv", sep="\t", index=False)

    direct_keys = set(
        zip(
            outliers.loc[outliers["aracne_direct"], "gene1"],
            outliers.loc[outliers["aracne_direct"], "gene2"],
        )
    )
    pair_keys = list(zip(pairs["gene1"], pairs["gene2"]))
    pairs["new_nmi_direct_outlier"] = [key in direct_keys for key in pair_keys]
    pairs["old_r2_candidate"] = (pairs["p_adj"] <= 0.05) & (pairs["mean_order_dist"] >= 50.09)
    overlap_table = pd.DataFrame(
        [
            {
                "set": "old_r2_candidate_p_adj_le_0.05_order_ge_50.09",
                "edges": int(pairs["old_r2_candidate"].sum()),
            },
            {
                "set": "new_weighted_nmi_direct_outlier",
                "edges": int(pairs["new_nmi_direct_outlier"].sum()),
            },
            {
                "set": "intersection",
                "edges": int((pairs["old_r2_candidate"] & pairs["new_nmi_direct_outlier"]).sum()),
            },
            {
                "set": "new_only",
                "edges": int((~pairs["old_r2_candidate"] & pairs["new_nmi_direct_outlier"]).sum()),
            },
        ]
    )
    overlap_table.to_csv(args.outdir / "r2_adj_vs_weighted_nmi_overlap.tsv", sep="\t", index=False)

    valid_compare = pairs[["weighted_nmi", "r2_adj"]].replace([np.inf, -np.inf], np.nan).dropna()
    comparison_summary = pd.DataFrame(
        [
            {
                "n_pairs": len(valid_compare),
                "pearson_weighted_nmi_vs_r2_adj": pearsonr(
                    valid_compare["weighted_nmi"], valid_compare["r2_adj"]
                ).statistic,
                "spearman_weighted_nmi_vs_r2_adj": spearmanr(
                    valid_compare["weighted_nmi"], valid_compare["r2_adj"]
                ).statistic,
                "weighted_nmi_order_elbow": ld_threshold,
                "r2_adj_order_elbow": r2_ld_threshold,
                "new_direct_edges": int(pairs["new_nmi_direct_outlier"].sum()),
                "new_direct_in_old_r2_candidates": int(
                    (pairs["old_r2_candidate"] & pairs["new_nmi_direct_outlier"]).sum()
                ),
            }
        ]
    )
    comparison_summary.to_csv(args.outdir / "r2_adj_vs_weighted_nmi_summary.tsv", sep="\t", index=False)

    direct = outliers[outliers["aracne_direct"]].nlargest(8, "weighted_nmi")
    plot_source = pairs.sample(min(args.plot_sample, len(pairs)), random_state=42)
    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": TOKENS["surface"],
            "axes.facecolor": TOKENS["panel"],
            "axes.edgecolor": TOKENS["axis"],
            "font.family": "sans-serif",
        },
    )
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.6))
    plot_panel(
        axes[0],
        plot_source,
        outliers,
        ld_threshold,
        outlier_threshold,
        extreme_threshold,
        None,
        "A. Full order-distance range",
    )
    plot_panel(
        axes[1],
        plot_source[plot_source["mean_order_dist"] <= 2.0 * ld_threshold],
        outliers[outliers["mean_order_dist"] <= 2.0 * ld_threshold],
        ld_threshold,
        outlier_threshold,
        extreme_threshold,
        2.0 * ld_threshold,
        "B. Linkage-decay region",
    )
    top_labels = []
    for rank, row in enumerate(direct.itertuples(index=False), start=1):
        label1 = row.gene1_annotation or row.gene1.split(":")[-1]
        label2 = row.gene2_annotation or row.gene2.split(":")[-1]
        top_labels.append(f"{rank}. {label1} - {label2}")
    axes[0].text(
        0.98,
        0.43,
        "Top direct outliers\n" + "\n".join(top_labels),
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=6.5,
        color=TOKENS["ink"],
        family="monospace",
        linespacing=1.25,
    )
    add_header(
        fig,
        axes[0],
        "SpydrPick-style wgMLST locus coupling by gene order distance",
        (
            f"{len(pairs):,} locus pairs; wg-content Jaccard reweighting (effective n={total_weight:.1f}); "
            f"Jeffreys pseudocount={args.pseudocount}; LD elbow={ld_threshold:.1f}. "
            "Horizontal thresholds use Tukey fences on per-locus maximum non-LD NMI; ARACNE marks direct and indirect outliers."
        ),
    )
    out_base = args.outdir / "wgmlst_weighted_nmi_vs_order_distance"
    fig.savefig(out_base.with_suffix(".png"), dpi=260, bbox_inches="tight", facecolor=TOKENS["surface"])
    fig.savefig(out_base.with_suffix(".svg"), bbox_inches="tight", facecolor=TOKENS["surface"])
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.4))
    compare_sample = pairs.sample(min(args.plot_sample, len(pairs)), random_state=42)
    axes[0].hexbin(
        compare_sample["weighted_nmi"],
        compare_sample["r2_adj"],
        gridsize=45,
        mincnt=1,
        bins="log",
        cmap=sns.blend_palette([TOKENS["panel"], TOKENS["blue_light"], TOKENS["blue"]], as_cmap=True),
    )
    direct_compare = pairs[pairs["new_nmi_direct_outlier"]]
    axes[0].scatter(
        direct_compare["weighted_nmi"],
        direct_compare["r2_adj"],
        s=10,
        facecolors="none",
        edgecolors=TOKENS["orange"],
        linewidths=0.55,
        alpha=0.75,
        label="New direct NMI outliers",
    )
    axes[0].set_xlabel("Weighted binary NMI")
    axes[0].set_ylabel("Previous R2_adj")
    axes[0].set_title("A. Pair-level agreement", loc="left", fontsize=11, fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    nmi_x = binned["x_mid"].to_numpy()
    nmi_y = binned["y_q_smooth"].to_numpy()
    r2_x = r2_binned["x_mid"].to_numpy()
    r2_y = r2_binned["y_q_smooth"].to_numpy()
    nmi_y = (nmi_y - nmi_y.min()) / max(nmi_y.max() - nmi_y.min(), 1e-12)
    r2_y = (r2_y - r2_y.min()) / max(r2_y.max() - r2_y.min(), 1e-12)
    axes[1].plot(nmi_x, nmi_y, color=TOKENS["blue"], linewidth=1.7, label="Weighted NMI q90")
    axes[1].plot(r2_x, r2_y, color=TOKENS["orange"], linewidth=1.7, label="R2_adj q90")
    axes[1].axvline(ld_threshold, color=TOKENS["blue"], linestyle=":", linewidth=1.1, label=f"NMI elbow={ld_threshold:.1f}")
    axes[1].axvline(r2_ld_threshold, color=TOKENS["orange"], linestyle="--", linewidth=1.1, label=f"R2 elbow={r2_ld_threshold:.1f}")
    axes[1].set_xlabel("Mean gene order distance")
    axes[1].set_ylabel("Normalized smoothed 90th percentile")
    axes[1].set_title("B. Distance-decay comparison", loc="left", fontsize=11, fontweight="bold")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(axis="y", color=TOKENS["grid"], linewidth=0.8)
        ax.grid(axis="x", visible=False)
        sns.despine(ax=ax)
    add_header(
        fig,
        axes[0],
        "Weighted wgMLST NMI compared with the previous R2_adj scan",
        (
            f"{len(valid_compare):,} matched gene pairs. "
            f"Pearson={comparison_summary.iloc[0]['pearson_weighted_nmi_vs_r2_adj']:.3f}; "
            f"Spearman={comparison_summary.iloc[0]['spearman_weighted_nmi_vs_r2_adj']:.3f}. "
            "Right panel normalizes each q90 trend to [0,1] so decay shapes, rather than absolute scales, are compared."
        ),
    )
    compare_base = args.outdir / "weighted_nmi_vs_previous_r2_adj"
    fig.savefig(compare_base.with_suffix(".png"), dpi=260, bbox_inches="tight", facecolor=TOKENS["surface"])
    fig.savefig(compare_base.with_suffix(".svg"), bbox_inches="tight", facecolor=TOKENS["surface"])
    plt.close(fig)

    thresholds = pd.DataFrame(
        [
            {"threshold": "ld_order_distance_elbow", "value": ld_threshold},
            {"threshold": "nmi_outlier_q3_plus_1.5_iqr", "value": outlier_threshold},
            {"threshold": "nmi_extreme_q3_plus_3_iqr", "value": extreme_threshold},
        ]
    )
    thresholds.to_csv(args.outdir / "thresholds.tsv", sep="\t", index=False)
    sensitivity_rows = []
    for candidate_ld in [50.0, 75.0, 100.0, 150.0, 200.0]:
        candidate_non_ld = pairs[pairs["mean_order_dist"] >= candidate_ld]
        candidate_maxima = pd.concat(
            [
                candidate_non_ld[["gene1", "weighted_nmi"]].rename(columns={"gene1": "gene"}),
                candidate_non_ld[["gene2", "weighted_nmi"]].rename(columns={"gene2": "gene"}),
            ],
            ignore_index=True,
        ).groupby("gene", as_index=False)["weighted_nmi"].max()
        candidate_q1, candidate_q3 = candidate_maxima["weighted_nmi"].quantile([0.25, 0.75])
        candidate_iqr = candidate_q3 - candidate_q1
        candidate_outlier_threshold = float(candidate_q3 + 1.5 * candidate_iqr)
        candidate_extreme_threshold = float(candidate_q3 + 3.0 * candidate_iqr)
        candidate_outliers = candidate_non_ld[
            candidate_non_ld["weighted_nmi"] > candidate_outlier_threshold
        ].copy()
        candidate_outliers = aracne_prune(candidate_outliers, "weighted_nmi")
        sensitivity_rows.append(
            {
                "ld_order_threshold": candidate_ld,
                "nmi_outlier_threshold": candidate_outlier_threshold,
                "nmi_extreme_threshold": candidate_extreme_threshold,
                "outlier_edges": len(candidate_outliers),
                "direct_outlier_edges": int(candidate_outliers["aracne_direct"].sum()),
                "indirect_outlier_edges": int((~candidate_outliers["aracne_direct"]).sum()),
                "extreme_outlier_edges": int(
                    (candidate_outliers["weighted_nmi"] > candidate_extreme_threshold).sum()
                ),
            }
        )
    pd.DataFrame(sensitivity_rows).to_csv(
        args.outdir / "ld_threshold_sensitivity.tsv",
        sep="\t",
        index=False,
    )
    summary = pd.DataFrame(
        [
            {
                "n_genomes": len(profile),
                "n_loci": len(genes),
                "n_pairs": len(pairs),
                "effective_sample_size": total_weight,
                "ld_order_threshold": ld_threshold,
                "nmi_outlier_threshold": outlier_threshold,
                "nmi_extreme_threshold": extreme_threshold,
                "r2_adj_order_threshold": r2_ld_threshold,
                "outlier_edges": len(outliers),
                "direct_outlier_edges": int(outliers["aracne_direct"].sum()),
                "indirect_outlier_edges": int((~outliers["aracne_direct"]).sum()),
                "extreme_outlier_edges": int(outliers["extreme_outlier"].sum()),
            }
        ]
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
