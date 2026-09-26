#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def signed_direction(df: pd.DataFrame, n_samples: float) -> pd.Series:
    obs_joint = df["n_joint"] / n_samples
    exp_joint = df["prev1"] * df["prev2"]
    delta = obs_joint - exp_joint
    return pd.Series(
        np.where(delta > 0, "positive", np.where(delta < 0, "negative", "zero")),
        index=df.index,
        name="raw_direction",
    )


def summarize_set(name: str, sub: pd.DataFrame) -> dict:
    return {
        "set": name,
        "n_edges": len(sub),
        "phylo_nonmissing": int(np.isfinite(sub["mean_phylo_dist"]).sum()),
        "median_phylo": sub["mean_phylo_dist"].median(),
        "q10_phylo": sub["mean_phylo_dist"].quantile(0.10),
        "q90_phylo": sub["mean_phylo_dist"].quantile(0.90),
        "median_order": sub["mean_order_dist"].median(),
        "median_bp": sub["mean_bp_dist"].median(),
        "median_r2_adj": sub["r2_adj"].median(),
        "mean_r2_adj": sub["r2_adj"].mean(),
    }


def summarize_filter(name: str, df: pd.DataFrame, mask: pd.Series) -> dict:
    sub = df.loc[mask]
    return {
        "filter": name,
        "n_edges": len(sub),
        "positive": int((sub["raw_direction"] == "positive").sum()),
        "negative": int((sub["raw_direction"] == "negative").sum()),
        "zero": int((sub["raw_direction"] == "zero").sum()),
        "median_phylo": sub["mean_phylo_dist"].median(),
        "median_order": sub["mean_order_dist"].median(),
        "median_bp": sub["mean_bp_dist"].median(),
        "median_r2_adj": sub["r2_adj"].median(),
    }


def binned_trend(df: pd.DataFrame, bins: int = 70) -> pd.DataFrame:
    rows = []
    for direction, sub in df.groupby("raw_direction"):
        if direction == "zero" or sub.empty:
            continue
        qs = np.linspace(0, 1, bins + 1)
        edges = np.unique(sub["mean_phylo_dist"].quantile(qs).to_numpy())
        if edges.size < 3:
            continue
        cats = pd.cut(sub["mean_phylo_dist"], bins=edges, include_lowest=True, duplicates="drop")
        grouped = sub.groupby(cats, observed=True)
        for _, g in grouped:
            if len(g) < 20:
                continue
            rows.append(
                {
                    "raw_direction": direction,
                    "mean_phylo_dist": g["mean_phylo_dist"].median(),
                    "r2_adj_median": g["r2_adj"].median(),
                    "r2_adj_q90": g["r2_adj"].quantile(0.90),
                    "n_edges": len(g),
                }
            )
    return pd.DataFrame(rows)


def plot_direction_decay(df: pd.DataFrame, trend: pd.DataFrame, phylo_thr: float, out_base: Path) -> None:
    rng = np.random.default_rng(42)
    plot_df = df[df["raw_direction"].isin(["positive", "negative"])].copy()
    if len(plot_df) > 200_000:
        plot_df = plot_df.iloc[rng.choice(len(plot_df), size=200_000, replace=False)]

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    palette = {"positive": "#1f77b4", "negative": "#d62728"}
    for direction, sub in plot_df.groupby("raw_direction"):
        ax.scatter(
            sub["mean_phylo_dist"],
            sub["r2_adj"],
            s=3,
            alpha=0.08,
            linewidths=0,
            color=palette[direction],
            label=f"{direction} sampled edges",
            rasterized=True,
        )
    for direction, sub in trend.groupby("raw_direction"):
        ax.plot(
            sub["mean_phylo_dist"],
            sub["r2_adj_q90"],
            color=palette[direction],
            linewidth=2.4,
            label=f"{direction} q90 trend",
        )
    ax.axvline(phylo_thr, color="black", linestyle="--", linewidth=1.5, label="phylo elbow")
    ax.set_xlabel("mean phylogenetic distance")
    ax.set_ylabel("adjusted r2")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, markerscale=3, fontsize=10)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=260)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Check sign-specific phylogenetic filtering on Cinch pair output.")
    parser.add_argument("--pairs", type=Path, default=Path("Leca/0427_1.ldlike_pairs.tsv.gz"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/sign_specific_phylo_check_0427_1"))
    parser.add_argument("--n-samples", type=float, default=1000.0)
    parser.add_argument("--p-adj", type=float, default=0.05)
    parser.add_argument("--order-threshold", type=float, default=50.09)
    parser.add_argument("--bp-threshold", type=float, default=105585.9)
    parser.add_argument("--phylo-threshold", type=float, default=0.020894)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cols = [
        "gene1",
        "gene2",
        "prev1",
        "prev2",
        "r2_adj",
        "p_adj",
        "n_joint",
        "mean_phylo_dist",
        "mean_bp_dist",
        "mean_order_dist",
    ]
    df = pd.read_csv(args.pairs, sep="\t", usecols=cols)
    for col in cols[2:]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df["p_adj"] <= args.p_adj].copy()
    df["raw_direction"] = signed_direction(df, args.n_samples)

    summary = pd.DataFrame(
        [
            summarize_set("all_sig", df),
            summarize_set("positive", df[df["raw_direction"] == "positive"]),
            summarize_set("negative", df[df["raw_direction"] == "negative"]),
        ]
    )
    summary.to_csv(args.outdir / "direction_phylo_summary.tsv", sep="\t", index=False)

    order = df["mean_order_dist"] >= args.order_threshold
    bp = df["mean_bp_dist"] >= args.bp_threshold
    pos = df["raw_direction"] == "positive"
    neg = df["raw_direction"] == "negative"
    phylo_long = df["mean_phylo_dist"] >= args.phylo_threshold
    phylo_short = df["mean_phylo_dist"] <= args.phylo_threshold

    filters = pd.DataFrame(
        [
            summarize_filter("order_only", df, order),
            summarize_filter("order_and_phylo_long_current", df, order & phylo_long),
            summarize_filter("order_sign_aware_pos_long_neg_short", df, order & ((pos & phylo_long) | (neg & phylo_short))),
            summarize_filter("order_sign_aware_pos_long_neg_any_phylo", df, order & ((pos & phylo_long) | neg)),
            summarize_filter("bp_only", df, bp),
            summarize_filter("bp_and_phylo_long_current", df, bp & phylo_long),
            summarize_filter("bp_sign_aware_pos_long_neg_short", df, bp & ((pos & phylo_long) | (neg & phylo_short))),
            summarize_filter("bp_sign_aware_pos_long_neg_any_phylo", df, bp & ((pos & phylo_long) | neg)),
        ]
    )
    filters.to_csv(args.outdir / "sign_specific_filter_counts.tsv", sep="\t", index=False)

    trend = binned_trend(df[np.isfinite(df["mean_phylo_dist"])])
    trend.to_csv(args.outdir / "direction_phylo_binned_trends.tsv", sep="\t", index=False)
    plot_direction_decay(df[np.isfinite(df["mean_phylo_dist"])], trend, args.phylo_threshold, args.outdir / "r2_adj_vs_phylo_by_direction")

    print(summary.to_string(index=False))
    print(filters.to_string(index=False))


if __name__ == "__main__":
    main()
