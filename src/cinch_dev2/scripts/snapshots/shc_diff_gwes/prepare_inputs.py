#!/usr/bin/env python3
"""Prepare frozen seed, background, block, and SHC-stratified split inputs."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WORKFLOW = ROOT / "Leca/workflows/shc_diff_gwes"
INPUT = WORKFLOW / "input"


def canonical(left: str, right: str) -> str:
    return "|".join(sorted((left, right)))


def main() -> None:
    INPUT.mkdir(parents=True, exist_ok=True)
    shc = pd.read_csv(ROOT / "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction/genome_shc_assignments.tsv", sep="\t")
    rng = np.random.default_rng(20260713)
    split = []
    for _, group in shc.groupby("shc", sort=True):
        genomes = group.genome.to_numpy().copy()
        rng.shuffle(genomes)
        cut = (len(genomes) + 1) // 2
        split.extend((g, "discovery" if i < cut else "validation") for i, g in enumerate(genomes))
    pd.DataFrame(split, columns=["genome", "split"]).merge(shc, on="genome", validate="one_to_one").to_csv(INPUT / "sample_split.tsv", sep="\t", index=False)

    paths = {
        "order": ROOT / "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction/order_dist_shc_strict_direct_edges.tsv",
        "bp": ROOT / "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction/bp_dist_shc_strict_direct_edges.tsv",
    }
    rows = []
    for source, path in paths.items():
        frame = pd.read_csv(path, sep="\t")
        for row in frame.itertuples(index=False):
            rows.append({"gene1": min(row.gene1, row.gene2), "gene2": max(row.gene1, row.gene2), "source": source, "conditional_mi": row.conditional_mi})
    seeds = pd.DataFrame(rows)
    seeds["pair_id"] = [canonical(a, b) for a, b in zip(seeds.gene1, seeds.gene2)]
    seeds = seeds.groupby(["pair_id", "gene1", "gene2"], as_index=False).agg(source=("source", lambda x: "+".join(sorted(set(x)))), conditional_mi=("conditional_mi", "max"))
    seeds.to_csv(INPUT / "seed_pairs.tsv", sep="\t", index=False)

    pairs = pd.read_csv(ROOT / "Leca/source/profile/0427_1.ldlike_pairs.tsv.gz", sep="\t", compression="gzip", usecols=["gene1", "gene2"])
    backgrounds = pd.DataFrame({"gene": sorted(set(pairs.gene1) | set(pairs.gene2))})
    backgrounds.to_csv(INPUT / "background_loci.tsv", sep="\t", index=False)

    positions = pd.read_csv(ROOT / "Leca/results/gene_gene_cgmlst/flattened_gene_positions.tsv", sep="\t")
    positions["block_id"] = positions.contig.astype(str) + ":" + (positions.midpoint // 100_000).astype(int).astype(str)
    block = backgrounds.merge(positions[["gene", "contig", "midpoint", "block_id"]], on="gene", how="left", validate="one_to_one")
    missing = block.block_id.isna()
    block.loc[missing, "block_id"] = "unplaced:" + block.loc[missing, "gene"]
    block.to_csv(INPUT / "locus_blocks.tsv", sep="\t", index=False)

    summary = pd.DataFrame([{"n_seed_pairs": len(seeds), "n_background_loci": len(backgrounds), "n_placed_loci": int((~missing).sum()), "n_unplaced_loci": int(missing.sum()), "n_shc": shc.shc.nunique(), "n_discovery": int((pd.DataFrame(split)[1] == "discovery").sum()), "n_validation": int((pd.DataFrame(split)[1] == "validation").sum()), "block_bp": 100_000}])
    summary.to_csv(INPUT / "preparation_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
