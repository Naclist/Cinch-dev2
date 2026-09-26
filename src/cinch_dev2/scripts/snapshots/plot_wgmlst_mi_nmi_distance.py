#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from secondary_filter_ld_decay import find_elbow, make_quantile_bins, smooth_series
from spydrpick_locus_wg_nmi import aracne_prune, read_annotations


COLORS = {
    "surface": "#FCFCFD",
    "panel": "#FFFFFF",
    "ink": "#20242C",
    "muted": "#697080",
    "grid": "#E5E8EE",
    "mi": "#267A6B",
    "nmi": "#4267B2",
    "threshold": "#C45D3A",
}


def metric_threshold(
    pairs: pd.DataFrame,
    metric: str,
    distance: str,
    distance_threshold: float,
) -> tuple[float, float, pd.DataFrame]:
    non_ld = pairs[pairs[distance] >= distance_threshold]
    maxima = pd.concat(
        [
            non_ld[["gene1", metric]].rename(columns={"gene1": "gene"}),
            non_ld[["gene2", metric]].rename(columns={"gene2": "gene"}),
        ],
        ignore_index=True,
    ).groupby("gene", as_index=False)[metric].max()
    q1, q3 = maxima[metric].quantile([0.25, 0.75])
    iqr = q3 - q1
    return float(q3 + 1.5 * iqr), float(q3 + 3.0 * iqr), maxima


def distance_elbow(pairs: pd.DataFrame, metric: str, distance: str) -> tuple[float, pd.DataFrame]:
    table = make_quantile_bins(pairs, distance, metric, 120, 0.90)
    table["y_median_smooth"] = smooth_series(table["y_median"])
    table["y_q_smooth"] = smooth_series(table["y_q"])
    elbow, table = find_elbow(table, log_x=False)
    return float(elbow), table


def plot_panel(
    ax: plt.Axes,
    background: pd.DataFrame,
    outliers: pd.DataFrame,
    metric: str,
    distance: str,
    ld_threshold: float,
    outlier_threshold: float,
    title: str,
) -> None:
    color = COLORS["mi"] if metric == "weighted_mi" else COLORS["nmi"]
    ax.scatter(
        background[distance],
        background[metric],
        s=3,
        alpha=0.14,
        color=color,
        edgecolors="none",
        rasterized=True,
    )
    if len(outliers):
        ax.scatter(
            outliers[distance],
            outliers[metric],
            s=9,
            alpha=0.72,
            color=COLORS["threshold"],
            edgecolors="none",
            rasterized=True,
            label="Tukey outlier",
        )
    ax.axvline(ld_threshold, color=COLORS["threshold"], linestyle=":", linewidth=1.2)
    ax.axhline(outlier_threshold, color=COLORS["threshold"], linestyle=":", linewidth=1.2)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel("Mean gene order distance" if distance == "mean_order_dist" else "Mean physical distance (bp)")
    ax.set_ylabel("Weighted MI (nats)" if metric == "weighted_mi" else "Weighted NMI")
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.8)
    ax.grid(axis="x", visible=False)
    sns.despine(ax=ax)


def quantile_trend(
    pairs: pd.DataFrame,
    distance: str,
    quantile: float,
    bins: int = 120,
) -> pd.DataFrame:
    ordered = pairs.sort_values(distance).reset_index(drop=True)
    chunks = np.array_split(np.arange(len(ordered)), bins)
    return pd.DataFrame(
        {
            "x": [ordered.iloc[index][distance].median() for index in chunks],
            "y": [ordered.iloc[index]["weighted_mi"].quantile(quantile) for index in chunks],
            "quantile": quantile,
        }
    )


def segmented_fit(trend: pd.DataFrame) -> tuple[float, np.ndarray, np.ndarray, int]:
    x = trend["x"].to_numpy(dtype=float)
    y = trend["y"].to_numpy(dtype=float)
    best = None
    for split in range(8, len(x) - 8):
        left = np.polyfit(x[: split + 1], y[: split + 1], 1)
        right = np.polyfit(x[split:], y[split:], 1)
        error = np.square(y[: split + 1] - np.polyval(left, x[: split + 1])).sum()
        error += np.square(y[split:] - np.polyval(right, x[split:])).sum()
        if best is None or error < best[0]:
            best = (error, float(x[split]), left, right, split)
    return best[1], best[2], best[3], best[4]


def make_mi_direct_figure(
    pairs: pd.DataFrame,
    background: pd.DataFrame,
    outdir: Path,
    gff: Path,
    distance: str,
    distance_label: str,
    prefix: str,
) -> None:
    valid = pairs.dropna(subset=[distance, "weighted_mi"])
    central_trend = quantile_trend(valid, distance, 0.90)
    upper_trend = quantile_trend(valid, distance, 0.995)
    central_elbow, central_left, central_right, central_split = segmented_fit(central_trend)
    upper_elbow, upper_left, upper_right, upper_split = segmented_fit(upper_trend)

    outlier_threshold, extreme_threshold, maxima = metric_threshold(
        valid, "weighted_mi", distance, upper_elbow
    )
    outliers = valid[
        (valid[distance] >= upper_elbow)
        & (valid["weighted_mi"] > outlier_threshold)
    ].copy()
    outliers = aracne_prune(outliers, "weighted_mi")
    outliers["extreme_outlier"] = outliers["weighted_mi"] > extreme_threshold
    outliers = outliers.sort_values(
        ["aracne_direct", "weighted_mi"], ascending=[False, False]
    )

    annotations = read_annotations(gff, set(valid["gene1"]) | set(valid["gene2"]))
    outliers["gene1_annotation"] = outliers["gene1"].map(annotations).fillna("")
    outliers["gene2_annotation"] = outliers["gene2"].map(annotations).fillna("")
    outliers.to_csv(outdir / f"{prefix}_mi_upper_tail_outliers_aracne.tsv", sep="\t", index=False)

    direct = outliers[outliers["aracne_direct"]]
    top20 = direct.nlargest(20, "weighted_mi").copy()
    top20.insert(0, "rank", np.arange(1, len(top20) + 1))
    top20.to_csv(outdir / f"{prefix}_mi_top20_direct_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "central_q90_segmented_elbow": central_elbow,
                "upper_q995_segmented_elbow": upper_elbow,
                "mi_outlier_threshold": outlier_threshold,
                "mi_extreme_threshold": extreme_threshold,
                "outlier_pairs": len(outliers),
                "direct_pairs": int(outliers["aracne_direct"].sum()),
                "indirect_pairs": int((~outliers["aracne_direct"]).sum()),
            }
        ]
    ).to_csv(outdir / f"{prefix}_mi_aracne_elbow_summary.tsv", sep="\t", index=False)

    fig = plt.figure(figsize=(18, 8.2))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.55, 1.05, 0.92], wspace=0.22)
    ax = fig.add_subplot(grid[0, 0])
    ax_trend = fig.add_subplot(grid[0, 1])
    ax_table = fig.add_subplot(grid[0, 2])

    bg = background.dropna(subset=[distance, "weighted_mi"])
    ax.scatter(
        bg[distance],
        bg["weighted_mi"],
        s=3,
        alpha=0.13,
        color=COLORS["mi"],
        edgecolors="none",
        rasterized=True,
    )
    indirect = outliers[~outliers["aracne_direct"]]
    ax.scatter(
        indirect[distance],
        indirect["weighted_mi"],
        s=13,
        alpha=0.65,
        color="#8A909A",
        edgecolors="none",
        label=f"Indirect ({len(indirect):,})",
    )
    ax.scatter(
        direct[distance],
        direct["weighted_mi"],
        s=15,
        alpha=0.80,
        color=COLORS["nmi"],
        edgecolors="none",
        label=f"Direct ({len(direct):,})",
    )
    ax.axvline(upper_elbow, color=COLORS["threshold"], linestyle=":", linewidth=1.4)
    ax.axhline(outlier_threshold, color=COLORS["threshold"], linestyle=":", linewidth=1.4)
    for row in top20.itertuples(index=False):
        ax.annotate(
            str(row.rank),
            (getattr(row, distance), row.weighted_mi),
            xytext=(2, 2),
            textcoords="offset points",
            fontsize=6.5,
            color=COLORS["ink"],
        )
    ax.set_title("A. Weighted MI and ARACNE classification", loc="left", fontsize=11, fontweight="bold")
    ax.set_xlabel(distance_label)
    ax.set_ylabel("Weighted MI (nats)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.8)
    ax.grid(axis="x", visible=False)
    sns.despine(ax=ax)

    for trend, color, label, elbow, left, right, split in [
        (central_trend, COLORS["mi"], "90th percentile", central_elbow, central_left, central_right, central_split),
        (upper_trend, COLORS["nmi"], "99.5th percentile", upper_elbow, upper_left, upper_right, upper_split),
    ]:
        x = trend["x"].to_numpy()
        ax_trend.scatter(x, trend["y"], s=13, alpha=0.55, color=color)
        ax_trend.plot(x[: split + 1], np.polyval(left, x[: split + 1]), color=color, linewidth=2)
        ax_trend.plot(x[split:], np.polyval(right, x[split:]), color=color, linewidth=2)
        ax_trend.axvline(elbow, color=color, linestyle=":", linewidth=1.3, label=f"{label}: {elbow:.1f}")
    ax_trend.set_title("B. Segmented fits to MI-distance decay", loc="left", fontsize=11, fontweight="bold")
    ax_trend.set_xlabel(distance_label)
    ax_trend.set_ylabel("Binned weighted MI")
    ax_trend.legend(frameon=False, fontsize=8)
    ax_trend.grid(axis="y", color=COLORS["grid"], linewidth=0.8)
    ax_trend.grid(axis="x", visible=False)
    sns.despine(ax=ax_trend)

    ax_table.axis("off")
    ax_table.set_title("C. Top 20 direct pairs", loc="left", fontsize=11, fontweight="bold")
    lines = []
    for row in top20.itertuples(index=False):
        left = row.gene1_annotation or row.gene1.split(":")[-1]
        right = row.gene2_annotation or row.gene2.split(":")[-1]
        lines.append(
            f"{row.rank:>2}. {left[:16]:<16}  {right[:16]:<16}  "
            f"MI={row.weighted_mi:.3f}  d={getattr(row, distance):.1f}"
        )
    ax_table.text(
        0,
        0.98,
        "\n".join(lines),
        transform=ax_table.transAxes,
        va="top",
        ha="left",
        fontsize=7.2,
        family="monospace",
        linespacing=1.35,
        color=COLORS["ink"],
    )

    fig.subplots_adjust(top=0.84, bottom=0.10, left=0.055, right=0.985)
    fig.suptitle(
        f"wgMLST weighted MI by {distance_label.lower()}: upper-tail linkage boundary and direct associations",
        x=0.055,
        y=0.965,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.text(
        0.055,
        0.91,
        f"Jaccard threshold=0.90. Central decay elbow={central_elbow:.1f}; "
        f"upper-tail elbow={upper_elbow:.1f} is used for LD exclusion. "
        f"MI Tukey fence={outlier_threshold:.3f}; ARACNE separates direct from triangle-mediated edges.",
        ha="left",
        fontsize=9.5,
        color=COLORS["muted"],
    )
    output = outdir / f"{prefix}_weighted_mi_direct_indirect_top20_and_elbow"
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight", facecolor=COLORS["surface"])
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor=COLORS["surface"])
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare weighted MI and NMI over order and physical distance.")
    parser.add_argument(
        "--pairs",
        type=Path,
        default=Path(
            "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/"
            "wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_nmi_distance_comparison"),
    )
    parser.add_argument("--plot-sample", type=int, default=350_000)
    parser.add_argument(
        "--gff",
        type=Path,
        default=Path("Leca/source/annotation/klebsiella_reference.gff"),
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    usecols = [
        "gene1",
        "gene2",
        "weighted_mi",
        "weighted_nmi",
        "mean_order_dist",
        "mean_bp_dist",
    ]
    pairs = pd.read_csv(args.pairs, sep="\t", usecols=usecols)
    background = pairs.sample(min(args.plot_sample, len(pairs)), random_state=42)

    rows = []
    trends = []
    settings = []
    for distance in ["mean_order_dist", "mean_bp_dist"]:
        for metric in ["weighted_mi", "weighted_nmi"]:
            valid = pairs.dropna(subset=[distance, metric])
            ld_threshold, trend = distance_elbow(valid, metric, distance)
            outlier_threshold, extreme_threshold, maxima = metric_threshold(
                valid, metric, distance, ld_threshold
            )
            outliers = valid[
                (valid[distance] >= ld_threshold) & (valid[metric] > outlier_threshold)
            ]
            key = f"{metric}__{distance}"
            trend.insert(0, "analysis", key)
            trends.append(trend)
            maxima.to_csv(args.outdir / f"{key}__per_locus_non_ld_maximum.tsv", sep="\t", index=False)
            outliers.to_csv(args.outdir / f"{key}__outliers.tsv.gz", sep="\t", index=False)
            rows.append(
                {
                    "analysis": key,
                    "n_pairs_with_distance": len(valid),
                    "ld_elbow": ld_threshold,
                    "outlier_threshold": outlier_threshold,
                    "extreme_threshold": extreme_threshold,
                    "outlier_pairs": len(outliers),
                    "pairs_above_0.3_beyond_ld": int(
                        ((valid[distance] >= ld_threshold) & (valid[metric] > 0.3)).sum()
                    ),
                    "q99_beyond_ld": valid.loc[valid[distance] >= ld_threshold, metric].quantile(0.99),
                    "q999_beyond_ld": valid.loc[valid[distance] >= ld_threshold, metric].quantile(0.999),
                }
            )
            settings.append((distance, metric, ld_threshold, outlier_threshold, outliers))

    pd.DataFrame(rows).to_csv(args.outdir / "mi_nmi_distance_summary.tsv", sep="\t", index=False)
    pd.concat(trends, ignore_index=True).to_csv(
        args.outdir / "mi_nmi_distance_binned_trends.tsv", sep="\t", index=False
    )

    sns.set_theme(
        style="whitegrid",
        rc={
            "figure.facecolor": COLORS["surface"],
            "axes.facecolor": COLORS["panel"],
            "axes.edgecolor": COLORS["grid"],
            "font.family": "sans-serif",
        },
    )
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.2))
    title_map = {
        ("mean_order_dist", "weighted_mi"): "A. MI by gene order distance",
        ("mean_order_dist", "weighted_nmi"): "B. NMI by gene order distance",
        ("mean_bp_dist", "weighted_mi"): "C. MI by physical distance",
        ("mean_bp_dist", "weighted_nmi"): "D. NMI by physical distance",
    }
    axis_map = {
        ("mean_order_dist", "weighted_mi"): axes[0, 0],
        ("mean_order_dist", "weighted_nmi"): axes[0, 1],
        ("mean_bp_dist", "weighted_mi"): axes[1, 0],
        ("mean_bp_dist", "weighted_nmi"): axes[1, 1],
    }
    for distance, metric, ld_threshold, outlier_threshold, outliers in settings:
        bg = background.dropna(subset=[distance, metric])
        plot_panel(
            axis_map[(distance, metric)],
            bg,
            outliers,
            metric,
            distance,
            ld_threshold,
            outlier_threshold,
            title_map[(distance, metric)],
        )

    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.07, right=0.98, hspace=0.30, wspace=0.18)
    fig.suptitle(
        "wgMLST locus association: raw weighted MI versus normalized MI",
        x=0.07,
        y=0.97,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color=COLORS["ink"],
    )
    fig.text(
        0.07,
        0.925,
        "Jaccard reweighting threshold = 0.90; identical locus pairs and sample weights in every panel. "
        "Orange points exceed the metric-specific Tukey fence after its distance elbow.",
        ha="left",
        fontsize=9.5,
        color=COLORS["muted"],
    )
    output = args.outdir / "weighted_mi_nmi_vs_order_and_bp_distance"
    fig.savefig(output.with_suffix(".png"), dpi=260, bbox_inches="tight", facecolor=COLORS["surface"])
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor=COLORS["surface"])
    plt.close(fig)
    make_mi_direct_figure(
        pairs,
        background,
        args.outdir,
        args.gff,
        "mean_order_dist",
        "Mean gene order distance",
        "order_dist",
    )
    make_mi_direct_figure(
        pairs,
        background,
        args.outdir,
        args.gff,
        "mean_bp_dist",
        "Mean physical distance (bp)",
        "bp_dist",
    )


if __name__ == "__main__":
    main()
