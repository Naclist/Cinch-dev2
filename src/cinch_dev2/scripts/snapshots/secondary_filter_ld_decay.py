#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


METRICS = {
    "mean_order_dist": {
        "label": "Mean order distance",
        "log_x": False,
    },
    "mean_bp_dist": {
        "label": "Mean bp distance",
        "log_x": True,
    },
    "mean_phylo_dist": {
        "label": "Mean phylogenetic distance",
        "log_x": True,
    },
}


def transform_x(x: np.ndarray, log_x: bool) -> np.ndarray:
    if log_x:
        return np.log10(np.maximum(x, 0) + 1.0)
    return x.astype(float)


def make_quantile_bins(df: pd.DataFrame, x_col: str, y_col: str, n_bins: int, q: float) -> pd.DataFrame:
    work = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
    work = work[work[x_col] >= 0]
    if work.empty:
        raise ValueError(f"No usable rows for {x_col}")

    unique_n = work[x_col].nunique()
    bins = min(n_bins, max(8, unique_n))
    work["_bin"] = pd.qcut(work[x_col].rank(method="first"), q=bins, duplicates="drop")
    binned = (
        work.groupby("_bin", observed=True)
        .agg(
            x_mid=(x_col, "median"),
            x_min=(x_col, "min"),
            x_max=(x_col, "max"),
            y_median=(y_col, "median"),
        )
    )
    # Pandas named aggregation cannot directly use a dynamic quantile and a count
    # cleanly across older versions, so fill those two columns explicitly.
    grouped = work.groupby("_bin", observed=True)[y_col]
    binned["y_q"] = grouped.quantile(q)
    binned["n"] = grouped.size()
    binned = binned.reset_index(drop=True)
    return binned


def smooth_series(y: pd.Series, window: int = 7) -> pd.Series:
    window = min(window, max(3, len(y) // 5 * 2 + 1))
    return y.rolling(window=window, min_periods=1, center=True).mean()


def find_elbow(binned: pd.DataFrame, log_x: bool, trend_col: str = "y_q_smooth") -> tuple[float, pd.DataFrame]:
    x = transform_x(binned["x_mid"].to_numpy(dtype=float), log_x)
    y = binned[trend_col].to_numpy(dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    idx_map = np.where(valid)[0]
    if len(x) < 5 or np.nanmax(y) == np.nanmin(y) or np.nanmax(x) == np.nanmin(x):
        elbow_idx = idx_map[len(idx_map) // 2]
        return float(binned.loc[elbow_idx, "x_mid"]), binned

    x_norm = (x - x.min()) / (x.max() - x.min())
    y_norm = (y - y.min()) / (y.max() - y.min())

    # For a decreasing LD-decay curve, the elbow is where the trend falls farthest
    # below the straight line connecting the high-LD and baseline ends.
    baseline_line = 1.0 - x_norm
    distance = baseline_line - y_norm
    distance[:2] = -np.inf
    distance[-2:] = -np.inf
    elbow_local = int(np.nanargmax(distance))
    elbow_idx = idx_map[elbow_local]

    binned = binned.copy()
    binned["elbow_score"] = np.nan
    binned.loc[idx_map, "elbow_score"] = distance
    return float(binned.loc[elbow_idx, "x_mid"]), binned


def plot_decay(
    df: pd.DataFrame,
    binned: pd.DataFrame,
    x_col: str,
    y_col: str,
    elbow: float,
    out_base: Path,
    label: str,
    log_x: bool,
    sample_n: int,
    random_state: int,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    plot_df = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    plot_df = plot_df[plot_df[x_col] >= 0]
    if len(plot_df) > sample_n:
        plot_df = plot_df.sample(sample_n, random_state=random_state)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(
        plot_df[x_col],
        plot_df[y_col],
        s=8,
        alpha=0.08,
        color="#3b6ea8",
        linewidths=0,
        label="significant pairs",
    )
    ax.plot(
        binned["x_mid"],
        binned["y_q_smooth"],
        color="#c23b22",
        linewidth=2.5,
        label="smoothed high-quantile trend",
    )
    ax.plot(
        binned["x_mid"],
        binned["y_median_smooth"],
        color="#2a9d8f",
        linewidth=2.0,
        alpha=0.85,
        label="smoothed median trend",
    )
    ax.axvline(elbow, color="#202020", linestyle="--", linewidth=1.7, label=f"elbow = {elbow:.4g}")
    if log_x:
        ax.set_xscale("symlog", linthresh=1)
        ax.set_xlim(left=0)
    ax.set_xlabel(label)
    ax.set_ylabel(y_col)
    ax.set_title(f"{y_col} decay vs {label}")
    ax.legend(frameon=False, fontsize=10)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find distance-decay elbows and filter significant gene-gene LD-like pairs."
    )
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("Leca/secondary_filter"))
    parser.add_argument("--r2-col", default="r2_adj")
    parser.add_argument("--p-col", default="p_adj")
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--quantile", type=float, default=0.90)
    parser.add_argument("--bins", type=int, default=120)
    parser.add_argument("--sample-n", type=int, default=180000)
    parser.add_argument("--random-state", type=int, default=42)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input_tsv, sep="\t")
    for col in [args.r2_col, args.p_col, *METRICS.keys()]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    sig = df[(df[args.p_col] <= args.p_cutoff) & np.isfinite(df[args.r2_col])].copy()
    if sig.empty:
        raise ValueError(f"No significant rows found with {args.p_col} <= {args.p_cutoff}")

    thresholds = []
    binned_tables = []
    for x_col, spec in METRICS.items():
        usable = sig[np.isfinite(sig[x_col])].copy()
        binned = make_quantile_bins(usable, x_col, args.r2_col, args.bins, args.quantile)
        binned["y_median_smooth"] = smooth_series(binned["y_median"])
        binned["y_q_smooth"] = smooth_series(binned["y_q"])
        elbow, binned = find_elbow(binned, spec["log_x"])
        binned.insert(0, "metric", x_col)
        binned_tables.append(binned)
        thresholds.append(
            {
                "metric": x_col,
                "threshold": elbow,
                "direction": "keep >= threshold",
                "n_rows_used": len(usable),
                "r2_column": args.r2_col,
                "p_column": args.p_col,
                "p_cutoff": args.p_cutoff,
                "trend_quantile": args.quantile,
                "log_x_for_elbow": spec["log_x"],
            }
        )
        plot_decay(
            usable,
            binned,
            x_col,
            args.r2_col,
            elbow,
            args.outdir / f"decay_{args.r2_col}_vs_{x_col}",
            spec["label"],
            spec["log_x"],
            args.sample_n,
            args.random_state,
        )

    threshold_df = pd.DataFrame(thresholds)
    threshold_df.to_csv(args.outdir / "elbow_thresholds.tsv", sep="\t", index=False)
    pd.concat(binned_tables, ignore_index=True).to_csv(args.outdir / "binned_decay_trends.tsv", sep="\t", index=False)

    keep = sig.copy()
    for row in threshold_df.itertuples(index=False):
        keep = keep[np.isfinite(keep[row.metric]) & (keep[row.metric] >= row.threshold)]

    sig.to_csv(args.outdir / "significant_pairs_before_secondary_filter.tsv.gz", sep="\t", index=False)
    keep.to_csv(args.outdir / "epistasis_candidate_pairs_after_elbow_filter.tsv.gz", sep="\t", index=False)

    summary = pd.Series(
        {
            "input_pairs": len(df),
            "significant_pairs": len(sig),
            "after_elbow_filter_pairs": len(keep),
            "removed_pairs": len(sig) - len(keep),
            "retained_fraction_of_significant": len(keep) / len(sig) if len(sig) else np.nan,
            "r2_column": args.r2_col,
            "p_column": args.p_col,
            "p_cutoff": args.p_cutoff,
            "trend_quantile": args.quantile,
        }
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", header=False)
    print(summary.to_string())
    print("\nElbow thresholds:")
    print(threshold_df[["metric", "threshold", "direction"]].to_string(index=False))


if __name__ == "__main__":
    main()
