#!/usr/bin/env python3
"""Add 20 known background-dependent synthetic triples to a copy of the Kp profile."""

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"
INPUT = WF / "input"


def main() -> None:
    profile = pd.read_csv(ROOT / "Leca/source/profile/0427_1.profile", sep="\t", index_col=0).fillna(0)
    split = pd.read_csv(INPUT / "sample_split.tsv", sep="\t").set_index("genome").loc[profile.index]
    rng = np.random.default_rng(20260713)
    truth = []
    synthetic_blocks = []
    toy_seeds = []
    toy_backgrounds = []

    for number in range(1, 21):
        prefix = f"toy_epi_{number:02d}"
        ga, gb, gc = f"{prefix}_A", f"{prefix}_B", f"{prefix}_C"
        a = np.zeros(len(profile), dtype=np.int16)
        b = np.zeros(len(profile), dtype=np.int16)
        c = np.zeros(len(profile), dtype=np.int16)
        for (_, _), group in split.reset_index().groupby(["shc", "split"], sort=True):
            idx = profile.index.get_indexer(group.genome)
            idx = idx.copy()
            rng.shuffle(idx)
            if len(idx) < 8:
                a[idx] = rng.integers(0, 2, len(idx))
                b[idx] = rng.integers(0, 2, len(idx))
                c[idx] = rng.integers(0, 2, len(idx))
                continue
            cut = len(idx) // 2
            zero, one = idx[:cut], idx[cut:]
            c[one] = 1
            one_pattern = np.arange(len(one)) % 2
            rng.shuffle(one_pattern)
            a[one] = one_pattern
            b[one] = one_pattern
            cells = np.tile(np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=np.int16), (int(np.ceil(len(zero) / 4)), 1))[: len(zero)]
            rng.shuffle(cells)
            a[zero], b[zero] = cells[:, 0], cells[:, 1]
        profile[ga], profile[gb], profile[gc] = a, b, c
        pair_id = "|".join(sorted((ga, gb)))
        truth.append({"pair_id": pair_id, "gene1": min(ga, gb), "gene2": max(ga, gb), "background": gc, "expected_direction": "enhancement", "injection": "B=A in C=1; balanced independent AB cells in C=0 within SHC and split"})
        toy_seeds.append({"pair_id": pair_id, "gene1": min(ga, gb), "gene2": max(ga, gb), "source": "toy_truth", "conditional_mi": np.nan})
        toy_backgrounds.append(gc)
        for label, gene in zip(("A", "B", "C"), (ga, gb, gc)):
            synthetic_blocks.append({"gene": gene, "contig": f"toy_contig_{number:02d}_{label}", "midpoint": float(number * 1_000_000 + {"A": 0, "B": 300_000, "C": 600_000}[label]), "block_id": f"toy_block_{number:02d}_{label}"})

    profile.to_csv(INPUT / "toy_profile.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame(truth).to_csv(INPUT / "toy_truth.tsv", sep="\t", index=False)
    negatives = []
    for number in range(1, 21):
        next_number = number % 20 + 1
        pair_id = "|".join(sorted((f"toy_epi_{number:02d}_A", f"toy_epi_{number:02d}_B")))
        negatives.append({"pair_id": pair_id, "background": f"toy_epi_{next_number:02d}_C", "control_type": "predeclared_negative"})
    pd.DataFrame(negatives).to_csv(INPUT / "toy_negative_controls.tsv", sep="\t", index=False)
    real_seeds = pd.read_csv(INPUT / "seed_pairs.tsv", sep="\t")
    pd.concat([real_seeds, pd.DataFrame(toy_seeds)], ignore_index=True).to_csv(INPUT / "toy_seed_pairs.tsv", sep="\t", index=False)
    real_backgrounds = pd.read_csv(INPUT / "background_loci.tsv", sep="\t")
    pd.concat([real_backgrounds, pd.DataFrame({"gene": toy_backgrounds})], ignore_index=True).to_csv(INPUT / "toy_background_loci.tsv", sep="\t", index=False)
    real_blocks = pd.read_csv(INPUT / "locus_blocks.tsv", sep="\t")
    pd.concat([real_blocks, pd.DataFrame(synthetic_blocks)], ignore_index=True).to_csv(INPUT / "toy_locus_blocks.tsv", sep="\t", index=False)
    summary = {"n_genomes": len(profile), "original_loci": profile.shape[1] - 60, "synthetic_loci": 60, "synthetic_triples": 20, "predeclared_negative_controls": 20, "seed": 20260713}
    (INPUT / "toy_injection_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(summary)


if __name__ == "__main__":
    main()
