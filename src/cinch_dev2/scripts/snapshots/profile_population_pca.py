#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pointbiserialr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


UNKNOWN_VALUES = {"", "unknown", "not_available", "not available", "missing", "na", "nan", "none", "null"}


def canonical_genome_name(name: str) -> str:
    match = re.search(r"(GC[AF]_\d+)", str(name))
    return match.group(1) if match else str(name)


def clean_value(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    return np.nan if text.lower() in UNKNOWN_VALUES else text


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip", index_col=False)
    meta["genome_key"] = meta["Query"].map(canonical_genome_name)
    meta = meta[meta["genome_key"].str.match(r"GC[AF]_\d+", na=False)].drop_duplicates("genome_key")
    for col in meta.columns:
        if col != "genome_key":
            meta[col] = meta[col].map(clean_value)
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description="PCA of genome x gene presence matrix and PC association with metadata levels.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/0427_1.profile"))
    parser.add_argument("--metadata", type=Path, default=Path("Leca/Kp.mp.merge"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/profile_population_pca"))
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--metadata-col", default="Geography")
    parser.add_argument("--metadata-level", default="India")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    prof = pd.read_csv(args.profile, sep="\t")
    prof["genome_key"] = prof["genome"].map(canonical_genome_name)
    genes = [c for c in prof.columns if c not in {"genome", "genome_key"}]
    x = (prof[genes] > 0).astype(np.float32)

    meta = load_metadata(args.metadata)
    merged_meta = prof[["genome", "genome_key"]].merge(meta, on="genome_key", how="left")
    y = (merged_meta[args.metadata_col].astype(str) == args.metadata_level).astype(int).to_numpy()

    n_components = min(args.n_pcs, x.shape[0] - 1, x.shape[1])
    x_scaled = StandardScaler(with_mean=True, with_std=True).fit_transform(x.values)
    pca = PCA(n_components=n_components, random_state=42)
    scores = pca.fit_transform(x_scaled)

    pc_cols = [f"PC{i}" for i in range(1, n_components + 1)]
    scores_df = pd.DataFrame(scores, columns=pc_cols)
    scores_df.insert(0, "genome_key", prof["genome_key"])
    scores_df.insert(0, "genome", prof["genome"])
    scores_df[args.metadata_col] = merged_meta[args.metadata_col].values
    scores_df[f"{args.metadata_col}={args.metadata_level}"] = y
    scores_df.to_csv(args.outdir / "genome_presence_pca_scores.tsv", sep="\t", index=False)

    var_df = pd.DataFrame(
        {
            "PC": pc_cols,
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        }
    )
    var_df.to_csv(args.outdir / "genome_presence_pca_variance.tsv", sep="\t", index=False)

    assoc_rows = []
    for i, pc in enumerate(pc_cols):
        r, p = pointbiserialr(y, scores[:, i])
        assoc_rows.append({"PC": pc, "r_with_level": r, "p_value": p, "explained_variance_ratio": pca.explained_variance_ratio_[i]})
    assoc = pd.DataFrame(assoc_rows)
    assoc["bonferroni_p"] = np.minimum(assoc["p_value"] * len(assoc), 1.0)
    assoc.to_csv(args.outdir / f"pc_association_{args.metadata_col}_{args.metadata_level}.tsv", sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(np.arange(1, n_components + 1), pca.explained_variance_ratio_ * 100, marker="o", linewidth=1.8)
    ax.set_xlabel("PC")
    ax.set_ylabel("Explained variance (%)")
    ax.set_title("Genome presence PCA scree")
    fig.tight_layout()
    fig.savefig(args.outdir / "genome_presence_pca_scree.png", dpi=240)
    fig.savefig(args.outdir / "genome_presence_pca_scree.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(np.arange(1, n_components + 1), np.cumsum(pca.explained_variance_ratio_) * 100, marker="o", linewidth=1.8)
    ax.set_xlabel("PC")
    ax.set_ylabel("Cumulative explained variance (%)")
    ax.set_title("Genome presence PCA cumulative variance")
    fig.tight_layout()
    fig.savefig(args.outdir / "genome_presence_pca_cumulative.png", dpi=240)
    fig.savefig(args.outdir / "genome_presence_pca_cumulative.pdf")
    plt.close(fig)

    plot_df = scores_df.copy()
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.scatterplot(
        data=plot_df,
        x="PC1",
        y="PC2",
        hue=f"{args.metadata_col}={args.metadata_level}",
        palette={0: "#6b7280", 1: "#c23b22"},
        alpha=0.85,
        s=44,
        ax=ax,
    )
    ax.set_title(f"Genome presence PCA colored by {args.metadata_col}={args.metadata_level}")
    fig.tight_layout()
    fig.savefig(args.outdir / f"genome_presence_pca_{args.metadata_col}_{args.metadata_level}.png", dpi=240)
    fig.savefig(args.outdir / f"genome_presence_pca_{args.metadata_col}_{args.metadata_level}.pdf")
    plt.close(fig)

    summary = pd.Series(
        {
            "genomes": x.shape[0],
            "genes": x.shape[1],
            "pcs": n_components,
            "level_positive_count": int(y.sum()),
            "pc1_variance": pca.explained_variance_ratio_[0],
            "pc2_variance": pca.explained_variance_ratio_[1],
            "pcs_to_50pct": int(np.searchsorted(np.cumsum(pca.explained_variance_ratio_), 0.50) + 1),
            "pcs_to_80pct": int(np.searchsorted(np.cumsum(pca.explained_variance_ratio_), 0.80) + 1),
        }
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", header=False)
    print(summary.to_string())
    print("\nTop PC associations:")
    print(assoc.sort_values("p_value").head(12).to_string(index=False))


if __name__ == "__main__":
    main()
