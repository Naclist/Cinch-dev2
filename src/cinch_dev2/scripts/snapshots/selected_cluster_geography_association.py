#!/usr/bin/env python3
"""Reproduce selected-cluster geography enrichment tables and a compact visual."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def canonical_genome(value: object) -> str:
    text = str(value).replace(".fna", "")
    match = re.search(r"(GC[AF]_\d+)", text)
    return match.group(1) if match else text


def bh_fdr(pvalues: pd.Series) -> np.ndarray:
    p = pd.to_numeric(pvalues, errors="coerce").fillna(1.0).to_numpy(float)
    n = len(p)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n, dtype=float)
    q[order] = np.minimum(adjusted, 1.0)
    return q


def read_metadata(path: Path) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("Query", "Request Failed")):
                continue
            parts = line.rstrip("\n").split("\t")
            # Historical Kp.mp.merge rows contain one more trailing field than
            # the repeated headers; Geography is stable at field 8.
            if len(parts) > 7 and re.match(r"GC[AF]_\d+", parts[0]):
                rows.append({"Query": parts[0], "Geography": parts[7]})
    if not rows:
        raise ValueError(f"No accession/geography rows found in {path}")
    metadata = pd.DataFrame(rows)
    metadata["genome_key"] = metadata["Query"].map(canonical_genome)
    return metadata.drop_duplicates("genome_key", keep="first")


def safe_level(level: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", level.lower()).strip("_")


def association_table(
    level: str,
    merged: pd.DataFrame,
    genes: list[str],
    cluster: pd.DataFrame,
    nmi: pd.DataFrame,
) -> pd.DataFrame:
    key = safe_level(level)
    other = f"non{key}"
    is_level = merged["Geography"].eq(level)
    rows: list[dict[str, object]] = []
    for gene in genes:
        present = pd.to_numeric(merged[gene], errors="coerce").fillna(0).ne(0)
        a = int((present & is_level).sum())
        b = int((~present & is_level).sum())
        c = int((present & ~is_level).sum())
        d = int((~present & ~is_level).sum())
        inside = a / (a + b) if a + b else np.nan
        outside = c / (c + d) if c + d else np.nan
        odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        rows.append(
            {
                "gene": gene,
                f"{key}_present": a,
                f"{key}_absent": b,
                f"{other}_present": c,
                f"{other}_absent": d,
                f"{key}_presence": inside,
                f"{other}_presence": outside,
                (f"presence_delta_{key}_minus_{other}" if key == "india" else f"delta_{key}_minus_{other}"): inside - outside,
                "odds_ratio": odds,
                "fisher_p": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    result["fdr_q"] = bh_fdr(result["fisher_p"])
    nmi_col = f"Geography={level}"
    if nmi_col not in nmi.columns:
        raise ValueError(f"Missing NMI column {nmi_col!r}")
    result = result.merge(nmi[nmi_col].rename(f"{key}_nmi"), left_on="gene", right_index=True, how="left")
    annotation_cols = [c for c in ["gene", "annotation_categories", "annotation_functions", "presence_fraction"] if c in cluster]
    return result.merge(cluster[annotation_cols], on="gene", how="left").sort_values(["fdr_q", "fisher_p"])


def plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    x = np.arange(len(summary))
    width = 0.34
    axes[0].bar(x - width / 2, summary["mean_level_presence"], width, label="Within geography", color="#3264a8")
    axes[0].bar(x + width / 2, summary["mean_other_presence"], width, label="Outside geography", color="white", edgecolor="#3264a8", hatch="///")
    axes[0].set_xticks(x, summary["level"])
    axes[0].set_ylabel("Mean selected-gene presence")
    axes[0].set_ylim(bottom=0)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].set_title("Presence enrichment")
    axes[1].bar(x, summary["mean_nmi"], color=["#c44e52", "#55a868"][: len(summary)])
    axes[1].scatter(x, summary["max_nmi"], marker="D", s=42, facecolor="white", edgecolor="black", label="Maximum gene NMI", zorder=3)
    axes[1].set_xticks(x, summary["level"])
    axes[1].set_ylabel("Gene–geography NMI")
    axes[1].set_ylim(bottom=0)
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].set_title("Association strength")
    fig.suptitle("Selected PCA cluster: geography associations", fontsize=14)
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        fig.savefig(outdir / f"selected_cluster_geography_summary.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-table", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--nmi", type=Path, required=True)
    parser.add_argument("--levels", nargs="+", default=["India", "Taiwan"])
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cluster = pd.read_csv(args.cluster_table, sep="\t")
    genes = cluster["gene"].astype(str).tolist()
    profile = pd.read_csv(args.profile, sep="\t", usecols=lambda c: c == "genome" or c in set(genes))
    profile["genome_key"] = profile["genome"].map(canonical_genome)
    metadata = read_metadata(args.metadata)
    merged = profile.merge(metadata[["genome_key", "Geography"]], on="genome_key", how="inner")
    nmi = pd.read_csv(args.nmi, sep="\t", index_col=0)

    summaries: list[dict[str, object]] = []
    for level in args.levels:
        key = safe_level(level)
        other = f"non{key}"
        table = association_table(level, merged, genes, cluster, nmi)
        table.to_csv(args.outdir / f"selected_cluster_{key}_association.tsv", sep="\t", index=False)
        summaries.append(
            {
                "level": level,
                "n_genomes_level": int(merged["Geography"].eq(level).sum()),
                "n_genomes_other": int(merged["Geography"].ne(level).sum()),
                "n_genes": len(table),
                "n_fdr_lt_0_05": int((table["fdr_q"] < 0.05).sum()),
                "mean_level_presence": table[f"{key}_presence"].mean(),
                "mean_other_presence": table[f"{other}_presence"].mean(),
                "mean_odds_ratio": table["odds_ratio"].replace([np.inf, -np.inf], np.nan).mean(),
                "mean_nmi": table[f"{key}_nmi"].mean(),
                "max_nmi": table[f"{key}_nmi"].max(),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.outdir / "selected_cluster_geography_summary.tsv", sep="\t", index=False)
    plot_summary(summary, args.outdir)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
