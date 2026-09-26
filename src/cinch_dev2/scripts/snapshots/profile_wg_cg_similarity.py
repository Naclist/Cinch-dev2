#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def allele_indicator(values: np.ndarray) -> sparse.csr_matrix:
    rows = []
    cols = []
    offset = 0
    for j in range(values.shape[1]):
        present = values[:, j] > 0
        if not np.any(present):
            continue
        row_idx = np.flatnonzero(present)
        alleles, codes = np.unique(values[present, j], return_inverse=True)
        rows.append(row_idx.astype(np.int32))
        cols.append((offset + codes).astype(np.int32))
        offset += len(alleles)
    if not rows:
        return sparse.csr_matrix((values.shape[0], 0), dtype=np.float32)
    row = np.concatenate(rows)
    col = np.concatenate(cols)
    data = np.ones(len(row), dtype=np.float32)
    return sparse.coo_matrix((data, (row, col)), shape=(values.shape[0], offset)).tocsr()


def similarity_matrices(values: np.ndarray, prevalence: np.ndarray, wg_min_prev: float, cg_min_prev: float):
    wg_keep = prevalence >= wg_min_prev
    cg_keep = prevalence >= cg_min_prev
    wg = values[:, wg_keep]
    cg = values[:, cg_keep]

    wg_presence = sparse.csr_matrix((wg > 0).astype(np.float32))
    wg_alleles = allele_indicator(wg)
    wg_intersection = (wg_presence @ wg_presence.T).toarray()
    wg_counts = np.asarray(wg_presence.sum(axis=1)).ravel()
    wg_union = wg_counts[:, None] + wg_counts[None, :] - wg_intersection
    wg_same_allele = (wg_alleles @ wg_alleles.T).toarray()
    wg_content_jaccard = np.divide(
        wg_intersection,
        wg_union,
        out=np.zeros_like(wg_intersection),
        where=wg_union > 0,
    )
    wg_allele_union = np.divide(
        wg_same_allele,
        wg_union,
        out=np.zeros_like(wg_same_allele),
        where=wg_union > 0,
    )
    wg_xor = wg_counts[:, None] + wg_counts[None, :] - 2.0 * wg_intersection
    wg_content_hamming_similarity = 1.0 - wg_xor / wg.shape[1]

    cg_presence = sparse.csr_matrix((cg > 0).astype(np.float32))
    cg_alleles = allele_indicator(cg)
    cg_callable = (cg_presence @ cg_presence.T).toarray()
    cg_same_allele = (cg_alleles @ cg_alleles.T).toarray()
    cg_allele_callable = np.divide(
        cg_same_allele,
        cg_callable,
        out=np.zeros_like(cg_same_allele),
        where=cg_callable > 0,
    )
    cgav_similarity_pseudocount = (cg_same_allele + 0.5) / (cg_callable + 1.0)
    cgav_distance = -np.log(cgav_similarity_pseudocount)
    return (
        wg_keep,
        cg_keep,
        wg_content_jaccard,
        wg_content_hamming_similarity,
        wg_allele_union,
        cg_allele_callable,
        cgav_similarity_pseudocount,
        cgav_distance,
    )


def weights_from_similarity(similarity: np.ndarray, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    neighbours = (similarity >= threshold).sum(axis=1).astype(np.int32)
    return neighbours, 1.0 / neighbours


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate wgMLST and cgMLST isolate similarities directly from an allele profile.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("Leca/results/profile_similarity/wg_cg_profile_0427_1"),
    )
    parser.add_argument("--wg-min-prevalence", type=float, default=0.01)
    parser.add_argument("--cg-min-prevalence", type=float, default=0.95)
    parser.add_argument("--reweighting-similarity", type=float, default=0.90)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    profile = pd.read_csv(args.profile, sep="\t", index_col=0).fillna(0)
    values = profile.to_numpy(dtype=np.int32)
    prevalence = (values > 0).mean(axis=0)
    (
        wg_keep,
        cg_keep,
        wg_content_jaccard,
        wg_content_hamming_similarity,
        wg_allele_union,
        cg_allele_callable,
        cgav_similarity_pseudocount,
        cgav_distance,
    ) = similarity_matrices(values, prevalence, args.wg_min_prevalence, args.cg_min_prevalence)

    np.savez_compressed(
        args.outdir / "wg_cg_similarity_matrices.npz",
        genomes=profile.index.to_numpy(dtype=str),
        wg_gene_content_jaccard=wg_content_jaccard,
        wg_gene_content_hamming_similarity=wg_content_hamming_similarity,
        wg_allele_concordance_union=wg_allele_union,
        cg_allele_concordance_callable=cg_allele_callable,
        cgav_similarity_pseudocount=cgav_similarity_pseudocount,
        cgav_distance=cgav_distance,
    )
    pd.DataFrame(
        {
            "gene": profile.columns,
            "prevalence": prevalence,
            "included_wg": wg_keep,
            "included_cg": cg_keep,
        }
    ).to_csv(args.outdir / "locus_sets.tsv", sep="\t", index=False)

    iu = np.triu_indices(len(profile), k=1)
    pairs = pd.DataFrame(
        {
            "genome1": profile.index.to_numpy()[iu[0]],
            "genome2": profile.index.to_numpy()[iu[1]],
            "wg_gene_content_jaccard": wg_content_jaccard[iu],
            "wg_gene_content_hamming_similarity": wg_content_hamming_similarity[iu],
            "wg_allele_concordance_union": wg_allele_union[iu],
            "cg_allele_concordance_callable": cg_allele_callable[iu],
            "cgav_similarity_pseudocount": cgav_similarity_pseudocount[iu],
            "cgav_distance": cgav_distance[iu],
        }
    )
    pairs.to_csv(args.outdir / "pairwise_wg_cg_similarities.tsv.gz", sep="\t", index=False)

    weights = pd.DataFrame({"genome": profile.index})
    effective_sizes = []
    for name, matrix in [
        ("wg_content", wg_content_jaccard),
        ("wg_content_hamming", wg_content_hamming_similarity),
        ("wg_allele", wg_allele_union),
        ("cg_allele", cg_allele_callable),
        ("cgav_pseudocount", cgav_similarity_pseudocount),
    ]:
        neighbours, sample_weights = weights_from_similarity(matrix, args.reweighting_similarity)
        weights[f"{name}_neighbours_ge_{args.reweighting_similarity:g}"] = neighbours
        weights[f"{name}_weight"] = sample_weights
        effective_sizes.append(
            {
                "similarity": name,
                "threshold": args.reweighting_similarity,
                "effective_sample_size_sum_weights": sample_weights.sum(),
                "min_neighbours": neighbours.min(),
                "median_neighbours": np.median(neighbours),
                "max_neighbours": neighbours.max(),
            }
        )
    weights.to_csv(args.outdir / "spydrpick_style_sample_weights.tsv", sep="\t", index=False)
    pd.DataFrame(effective_sizes).to_csv(args.outdir / "sample_weight_summary.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_genomes": len(profile),
                "profile_loci": profile.shape[1],
                "wg_loci_prevalence_ge_cutoff": int(wg_keep.sum()),
                "wg_min_prevalence": args.wg_min_prevalence,
                "cg_loci_prevalence_ge_cutoff": int(cg_keep.sum()),
                "cg_min_prevalence": args.cg_min_prevalence,
                "reweighting_similarity": args.reweighting_similarity,
                "wg_content_pair_median": float(np.median(wg_content_jaccard[iu])),
                "wg_content_hamming_pair_median": float(np.median(wg_content_hamming_similarity[iu])),
                "wg_allele_pair_median": float(np.median(wg_allele_union[iu])),
                "cg_allele_pair_median": float(np.median(cg_allele_callable[iu])),
                "cgav_distance_pair_median": float(np.median(cgav_distance[iu])),
            }
        ]
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    print("\nWeights:")
    print(pd.DataFrame(effective_sizes).to_string(index=False))


if __name__ == "__main__":
    main()
