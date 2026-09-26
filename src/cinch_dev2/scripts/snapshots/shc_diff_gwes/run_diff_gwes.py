#!/usr/bin/env python3
"""Run frozen SHC-conditioned Diff-GWES on true or injected-toy Kp profiles."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"
LN2 = np.log(2.0)


def bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def weighted_mi(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    total = w.sum()
    if total <= 0:
        return np.nan
    cells = np.array([w[(x == a) & (y == b)].sum() for a, b in [(0, 0), (0, 1), (1, 0), (1, 1)]]) / total
    px = np.array([cells[0] + cells[1], cells[2] + cells[3]])
    py = np.array([cells[0] + cells[2], cells[1] + cells[3]])
    expected = np.array([px[0] * py[0], px[0] * py[1], px[1] * py[0], px[1] * py[1]])
    valid = (cells > 0) & (expected > 0)
    return float(np.sum(cells[valid] * np.log(cells[valid] / expected[valid])))


def background_stat(
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    w: np.ndarray,
    shc: np.ndarray,
    min_group: int = 25,
    min_shc_cell: int = 3,
    min_informative: int = 2,
    cluster_indices: list[np.ndarray] | None = None,
) -> dict:
    result = {"n0": int((c == 0).sum()), "n1": int((c == 1).sum())}
    for state in (0, 1):
        numerator = 0.0
        denominator = 0.0
        informative = 0
        groups = cluster_indices if cluster_indices is not None else [np.flatnonzero(shc == cluster) for cluster in np.unique(shc)]
        for idx in groups:
            selected = idx[c[idx] == state]
            n = len(selected)
            if n < min_shc_cell:
                continue
            sx, sy = int(x[selected].sum()), int(y[selected].sum())
            if not (0 < sx < n and 0 < sy < n):
                continue
            mass = float(w[selected].sum())
            mi = weighted_mi(x[selected], y[selected], w[selected])
            if np.isfinite(mi) and mass > 0:
                numerator += mass * mi
                denominator += mass
                informative += 1
        conditional_mi = numerator / denominator if denominator > 0 else np.nan
        result[f"informative_shc_{state}"] = informative
        result[f"conditional_mi_{state}"] = conditional_mi
        result[f"epidis_{state}"] = np.sqrt(max(conditional_mi, 0.0) / LN2) if np.isfinite(conditional_mi) else np.nan
    result["eligible"] = bool(result["n0"] >= min_group and result["n1"] >= min_group and result["informative_shc_0"] >= min_informative and result["informative_shc_1"] >= min_informative)
    result["delta"] = result["epidis_1"] - result["epidis_0"] if result["eligible"] else np.nan
    return result


def stable_seed(text: str) -> int:
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def make_cluster_indices(shc: np.ndarray, minimum_size: int = 2) -> list[np.ndarray]:
    return [idx for cluster in np.unique(shc) if len(idx := np.flatnonzero(shc == cluster)) >= minimum_size]


def permuted_background(c: np.ndarray, cluster_indices: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    out = c.copy()
    for idx in cluster_indices:
        out[idx] = c[idx][rng.permutation(len(idx))]
    return out


def load_mode(mode: str):
    inp = WF / "input"
    if mode == "true":
        return ROOT / "Leca/source/profile/0427_1.profile", inp / "seed_pairs.tsv", inp / "background_loci.tsv", inp / "locus_blocks.tsv", None
    return inp / "toy_profile.tsv.gz", inp / "toy_seed_pairs.tsv", inp / "toy_background_loci.tsv", inp / "toy_locus_blocks.tsv", inp / "toy_truth.tsv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["true", "toy"], required=True)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--effect-screen", type=float, default=0.15)
    parser.add_argument("--max-backgrounds-per-seed", type=int, default=25)
    args = parser.parse_args()
    outdir = WF / "output" / args.mode
    outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s\t%(levelname)s\t%(message)s", handlers=[logging.FileHandler(outdir / "run.log", mode="w"), logging.StreamHandler()])
    log = logging.getLogger(args.mode)

    profile_path, seeds_path, backgrounds_path, blocks_path, truth_path = load_mode(args.mode)
    seeds = pd.read_csv(seeds_path, sep="\t")
    backgrounds = pd.read_csv(backgrounds_path, sep="\t").gene.tolist()
    blocks = pd.read_csv(blocks_path, sep="\t").set_index("gene")
    needed = sorted(set(backgrounds) | set(seeds.gene1) | set(seeds.gene2))
    profile = pd.read_csv(profile_path, sep="\t", usecols=lambda c: c == "genome" or c in set(needed), index_col=0).fillna(0)
    binary = (profile > 0).astype(np.uint8)
    shc_table = pd.read_csv(WF / "input/sample_split.tsv", sep="\t").set_index("genome").loc[profile.index]
    weight_table = pd.read_csv(ROOT / "Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv", sep="\t").set_index("genome").loc[profile.index]
    weights = weight_table.wg_content_hamming_weight.to_numpy(float)
    weights /= weights.sum()

    prevalence = binary[backgrounds].mean()
    backgrounds = prevalence[(prevalence >= 0.05) & (prevalence <= 0.95)].index.tolist()
    log.info("%s: %d seeds, %d prevalence-eligible backgrounds", args.mode, len(seeds), len(backgrounds))
    records = []
    discovery_mask = shc_table.split.eq("discovery").to_numpy()
    validation_mask = ~discovery_mask
    discovery_clusters = shc_table.shc.to_numpy()[discovery_mask]
    validation_clusters = shc_table.shc.to_numpy()[validation_mask]
    discovery_groups = make_cluster_indices(discovery_clusters, 3)
    validation_groups = make_cluster_indices(validation_clusters, 3)
    validation_permutation_groups = make_cluster_indices(validation_clusters, 2)
    for seed in seeds.itertuples(index=False):
        x = binary[seed.gene1].to_numpy(np.uint8)[discovery_mask]
        y = binary[seed.gene2].to_numpy(np.uint8)[discovery_mask]
        seed_blocks = {blocks.at[g, "block_id"] for g in (seed.gene1, seed.gene2) if g in blocks.index}
        for gene in backgrounds:
            if gene in (seed.gene1, seed.gene2):
                continue
            block_id = blocks.at[gene, "block_id"] if gene in blocks.index else f"unplaced:{gene}"
            if block_id in seed_blocks:
                continue
            stat = background_stat(x, y, binary[gene].to_numpy(np.uint8)[discovery_mask], weights[discovery_mask], discovery_clusters, cluster_indices=discovery_groups)
            records.append({"pair_id": seed.pair_id, "gene1": seed.gene1, "gene2": seed.gene2, "background": gene, "background_block": block_id, **{f"discovery_{k}": v for k, v in stat.items()}})
    discovery = pd.DataFrame(records)
    discovery.to_csv(outdir / "discovery_all_triples.tsv.gz", sep="\t", index=False, compression="gzip")
    candidates = discovery[discovery.discovery_eligible & (discovery.discovery_delta.abs() >= args.effect_screen)].copy()
    candidates["abs_discovery_delta"] = candidates.discovery_delta.abs()
    candidates = candidates.sort_values(["pair_id", "abs_discovery_delta"], ascending=[True, False]).groupby("pair_id", as_index=False, group_keys=False).head(args.max_backgrounds_per_seed)
    negative_keys = set()
    if truth_path is not None:
        negative_path = WF / "input/toy_negative_controls.tsv"
        negative = pd.read_csv(negative_path, sep="\t")
        negative_keys = set(zip(negative.pair_id, negative.background))
        forced = discovery[[((a, b) in negative_keys) for a, b in zip(discovery.pair_id, discovery.background)]].copy()
        forced["abs_discovery_delta"] = forced.discovery_delta.abs()
        candidates = pd.concat([candidates, forced], ignore_index=True).drop_duplicates(["pair_id", "background"])
    log.info("Selected %d validation candidates", len(candidates))

    validation_rows = []
    for row in candidates.itertuples(index=False):
        x = binary[row.gene1].to_numpy(np.uint8)[validation_mask]
        y = binary[row.gene2].to_numpy(np.uint8)[validation_mask]
        c = binary[row.background].to_numpy(np.uint8)[validation_mask]
        w = weights[validation_mask]
        clusters = validation_clusters
        observed = background_stat(x, y, c, w, clusters, cluster_indices=validation_groups)
        record = row._asdict() | {f"validation_{k}": v for k, v in observed.items()}
        if observed["eligible"]:
            null = np.empty(args.permutations)
            rng = np.random.default_rng(stable_seed(f"20260713|{row.background_block}"))
            for p in range(args.permutations):
                cp = permuted_background(c, validation_permutation_groups, rng)
                null[p] = background_stat(x, y, cp, w, clusters, cluster_indices=validation_groups)["delta"]
            null = null[np.isfinite(null)]
            record.update({"null_n": len(null), "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)), "perm_p": float((1 + np.sum(np.abs(null) >= abs(observed["delta"]))) / (len(null) + 1)), "z": float((observed["delta"] - null.mean()) / null.std(ddof=1)) if null.std(ddof=1) > 0 else np.nan})
        validation_rows.append(record)
    validation = pd.DataFrame(validation_rows)
    validation["perm_q"] = np.nan
    eligible = validation.validation_eligible.fillna(False) if len(validation) else pd.Series(dtype=bool)
    if eligible.any():
        validation.loc[eligible, "perm_q"] = bh(validation.loc[eligible, "perm_p"].to_numpy(float))
    validation["direction_replicated"] = np.sign(validation.discovery_delta) == np.sign(validation.validation_delta)
    validation["passes"] = validation.validation_eligible.fillna(False) & validation.direction_replicated.fillna(False) & (validation.perm_q < 0.05)
    if truth_path is not None:
        truth = pd.read_csv(truth_path, sep="\t")
        truth_keys = set(zip(truth.pair_id, truth.background))
        validation["is_injected_truth"] = [(a, b) in truth_keys for a, b in zip(validation.pair_id, validation.background)]
        validation["is_negative_control"] = [(a, b) in negative_keys for a, b in zip(validation.pair_id, validation.background)]
    else:
        validation["is_injected_truth"] = False
        validation["is_negative_control"] = False
    validation.to_csv(outdir / "validation_candidates.tsv", sep="\t", index=False)
    validation[validation.passes].sort_values("perm_q").to_csv(outdir / "significant_triples.tsv", sep="\t", index=False)

    summary = {"mode": args.mode, "n_seed_pairs": len(seeds), "n_backgrounds_prevalence_eligible": len(backgrounds), "n_discovery_tests": len(discovery), "n_discovery_eligible": int(discovery.discovery_eligible.sum()), "n_validation_candidates": len(validation), "n_validation_eligible": int(eligible.sum()), "n_significant": int(validation.passes.sum()), "permutations": args.permutations, "effect_screen": args.effect_screen, "max_backgrounds_per_seed": args.max_backgrounds_per_seed}
    if truth_path is not None:
        summary.update({"n_injected_truth": len(truth), "truth_selected": int(validation.is_injected_truth.sum()), "truth_significant": int((validation.is_injected_truth & validation.passes).sum()), "n_negative_controls": len(negative_keys), "negative_controls_tested": int(validation.is_negative_control.sum()), "negative_controls_significant": int((validation.is_negative_control & validation.passes).sum()), "significant_nontruth": int((~validation.is_injected_truth & validation.passes).sum())})
    pd.DataFrame([summary]).to_csv(outdir / "summary.tsv", sep="\t", index=False)
    (outdir / "run_provenance.json").write_text(json.dumps(summary | {"profile": str(profile_path), "seed_file": str(seeds_path), "background_file": str(backgrounds_path), "block_file": str(blocks_path), "sample_split_seed": 20260713, "permutation_seed_rule": "sha256(20260713|block_id)", "statistic": "sqrt(I(A;B|SHC,C=c)/ln2) difference"}, indent=2) + "\n")
    log.info("Completed: %s", summary)


if __name__ == "__main__":
    main()
