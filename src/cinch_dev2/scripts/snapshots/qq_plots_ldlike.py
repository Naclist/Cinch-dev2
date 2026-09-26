#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import chi2


def prepare_pvalues(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    p = p[np.isfinite(p)]
    p = p[(p >= 0) & (p <= 1)]
    if p.size == 0:
        raise ValueError("No valid p-values found.")
    positive = p[p > 0]
    floor = positive.min() / 10.0 if positive.size else np.nextafter(0, 1)
    p[p <= 0] = max(floor, np.nextafter(0, 1))
    p[p > 1] = 1.0
    return p


def qq_data(p: np.ndarray, max_points: int) -> pd.DataFrame:
    p_sorted = np.sort(p)
    n = p_sorted.size
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p_sorted)
    if n > max_points:
        idx = np.unique(np.linspace(0, n - 1, max_points).astype(int))
        expected = expected[idx]
        observed = observed[idx]
    return pd.DataFrame({"expected": expected, "observed": observed})


def nice_axis_limit(value: float, step: float = 1.0) -> float:
    if not np.isfinite(value) or value <= 0:
        return step
    return float(np.ceil(value / step) * step)


def lambda_gc(p: np.ndarray) -> float:
    p = np.clip(p, np.nextafter(0, 1), 1.0)
    chisq = chi2.isf(p, df=1)
    return float(np.nanmedian(chisq) / chi2.ppf(0.5, df=1))


def plot_qq(p: np.ndarray, title: str, out_base: Path, max_points: int) -> dict[str, float]:
    df = qq_data(p, max_points)
    expected_max = -np.log10(0.5 / p.size)
    x_limit = nice_axis_limit(expected_max * 1.02, step=0.5)
    y_limit = nice_axis_limit(df["observed"].max() * 1.02, step=1.0)

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    ax.scatter(
        df["expected"],
        df["observed"],
        s=13,
        alpha=0.35,
        linewidths=0,
        color="black",
        rasterized=True,
    )
    ax.plot([0, x_limit], [0, x_limit], color="red", linewidth=2.0)
    ax.set_xlabel(r"Expected $-\log_{10}(pvalue)$")
    ax.set_ylabel(r"Observed $-\log_{10}(pvalue)$")
    ax.set_title(title)
    ax.set_xlim(0, x_limit)
    ax.set_ylim(0, y_limit)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)

    return {
        "n_pvalues": int(p.size),
        "min_p": float(np.min(p)),
        "median_p": float(np.median(p)),
        "lambda_gc": lambda_gc(p),
        "p_lt_0.05": int(np.sum(p < 0.05)),
        "p_lt_1e_5": int(np.sum(p < 1e-5)),
        "p_lt_1e_10": int(np.sum(p < 1e-10)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="QQ plots for raw, adjusted, and distance-filtered LD-like p-values.")
    parser.add_argument("--all-pairs", type=Path, default=Path("Leca/0427_1.ldlike_pairs.tsv.gz"))
    parser.add_argument(
        "--filtered-pairs",
        type=Path,
        default=Path("Leca/secondary_filter_0427_1/epistasis_candidate_pairs_after_elbow_filter.tsv.gz"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("Leca/qq_plots_0427_1"))
    parser.add_argument("--max-points", type=int, default=250000)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_df = pd.read_csv(args.all_pairs, sep="\t", usecols=["p_raw", "p_adj"])
    filt_df = pd.read_csv(args.filtered_pairs, sep="\t", usecols=["p_adj"])

    jobs = [
        ("all_pairs_p_raw", prepare_pvalues(all_df["p_raw"]), "All pairs p_raw QQ"),
        ("all_pairs_p_adj", prepare_pvalues(all_df["p_adj"]), "All pairs p_adj QQ"),
        ("distance_filtered_p_adj", prepare_pvalues(filt_df["p_adj"]), "Distance-filtered candidate p_adj QQ"),
    ]

    rows = []
    for name, p, title in jobs:
        stats = plot_qq(p, title, args.outdir / name, args.max_points)
        rows.append({"plot": name, **stats})

    summary = pd.DataFrame(rows)
    summary.to_csv(args.outdir / "qq_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
