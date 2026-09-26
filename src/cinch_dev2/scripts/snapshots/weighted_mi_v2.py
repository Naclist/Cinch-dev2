#!/usr/bin/env python3
"""Compute thesis-style sample-reweighted binary MI on the frozen Kp pair universe."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def file_sha256(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def weighted_binary_information(
    c11: np.ndarray,
    m1: np.ndarray,
    m2: np.ndarray,
    total_weight: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return MI (nats), symmetric NMI, H(X), and H(Y) without pseudocounts."""
    counts = np.column_stack(
        [
            total_weight - m1 - m2 + c11,
            m2 - c11,
            m1 - c11,
            c11,
        ]
    )
    counts = np.clip(counts, 0.0, total_weight)
    probabilities = counts / total_weight
    px = np.column_stack(
        [probabilities[:, 0] + probabilities[:, 1], probabilities[:, 2] + probabilities[:, 3]]
    )
    py = np.column_stack(
        [probabilities[:, 0] + probabilities[:, 2], probabilities[:, 1] + probabilities[:, 3]]
    )
    expected = np.column_stack(
        [px[:, 0] * py[:, 0], px[:, 0] * py[:, 1], px[:, 1] * py[:, 0], px[:, 1] * py[:, 1]]
    )
    terms = np.zeros_like(probabilities)
    valid = (probabilities > 0.0) & (expected > 0.0)
    terms[valid] = probabilities[valid] * np.log(probabilities[valid] / expected[valid])
    mi = terms.sum(axis=1)

    hx_terms = np.zeros_like(px)
    hy_terms = np.zeros_like(py)
    px_valid = px > 0.0
    py_valid = py > 0.0
    hx_terms[px_valid] = -px[px_valid] * np.log(px[px_valid])
    hy_terms[py_valid] = -py[py_valid] * np.log(py[py_valid])
    hx = hx_terms.sum(axis=1)
    hy = hy_terms.sum(axis=1)
    nmi = np.divide(2.0 * mi, hx + hy, out=np.zeros_like(mi), where=(hx + hy) > 0.0)
    return mi, nmi, hx, hy


def rank_correlation(x: pd.Series, y: pd.Series) -> float:
    return float(x.rank(method="average").corr(y.rank(method="average"), method="pearson"))


def run_self_tests() -> None:
    mi, nmi, _, _ = weighted_binary_information(
        np.array([0.5]), np.array([0.5]), np.array([0.5]), 1.0
    )
    assert np.isclose(mi[0], np.log(2.0))
    assert np.isclose(nmi[0], 1.0)

    mi, nmi, _, _ = weighted_binary_information(
        np.array([0.25]), np.array([0.5]), np.array([0.5]), 1.0
    )
    assert np.isclose(mi[0], 0.0, atol=1e-15)
    assert np.isclose(nmi[0], 0.0, atol=1e-15)


def write_field_dictionary(path: Path) -> None:
    path.write_text(
        """# weighted MI v2 field dictionary

- `mi`: sample-reweighted mutual information of binary locus presence states, in natural-log units (nats), estimated without a pseudocount.
- `nmi`: symmetric normalized MI, `2 * mi / (H(X) + H(Y))`, using the same weighted empirical probabilities.
- `gene1_frequency_0`, `gene1_frequency_1`: thesis equation (2-11) `beta` and `alpha` for `gene1`; weighted absence and presence proportions after sample weights are normalized to sum one.
- `gene2_frequency_0`, `gene2_frequency_1`: corresponding weighted state proportions for `gene2`.
- `frequency_weight`: compatibility scalar fixed to `1.0`. Thesis equation (2-11) supplies a two-component marginal mixing vector, not a second scalar multiplier for MI.
- `frequency_weighted_MI`: equal to `mi`. The alpha/beta locus frequencies are already inside the joint and marginal probabilities defining MI; multiplying MI again would double-weight frequency and would not reproduce the audited formula.
- `sample_weight_scheme`: provenance label for the sample weights.
- `sample_weight_tau`: Hamming-distance neighbourhood cutoff. Similarity >=0.90 is distance <=0.10 in the frozen upstream calculation.
- `sample_weight_sum_raw`: sum of unnormalized inverse-neighbour weights; reported as an effective-size diagnostic.
- `sample_weight_sum_normalized`: sum after thesis equation (2-7) normalization; expected to be 1.

The retained source fields (`gene1`, `gene2`, `r2_adj`, distances, and related columns) keep their definitions from `kp_universe_spec.md`. This is weighted MI v2, not EpiDis.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument(
        "--pairs", type=Path, default=Path("Leca/source/profile/0427_1.ldlike_pairs.tsv.gz")
    )
    parser.add_argument(
        "--weights",
        type=Path,
        default=Path(
            "Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv"
        ),
    )
    parser.add_argument(
        "--old-mi",
        type=Path,
        default=Path(
            "Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/"
            "wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("Leca/results/gene_gene_wgmlst/weighted_mi_v2"),
    )
    parser.add_argument("--weight-column", default="wg_content_hamming_weight")
    parser.add_argument("--tau", type=float, default=0.10)
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()

    run_self_tests()
    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger("weighted_mi_v2")
    log.info("Starting weighted MI v2")
    log.info("Command: %s", " ".join(sys.argv))

    pairs = pd.read_csv(args.pairs, sep="\t", compression="gzip")
    genes = sorted(set(pairs["gene1"]) | set(pairs["gene2"]))
    log.info("Frozen pair universe: %d pairs, %d loci", len(pairs), len(genes))
    if len(pairs) != 2_609_249 or len(genes) != 3_279:
        raise ValueError("Pair universe does not match kp_universe_spec.md freeze")
    if pairs.duplicated(["gene1", "gene2"]).any():
        raise ValueError("Duplicate ordered gene pairs found in frozen pair table")

    profile_header = pd.read_csv(args.profile, sep="\t", nrows=0).columns.tolist()
    missing_genes = sorted(set(genes) - set(profile_header[1:]))
    if missing_genes:
        raise ValueError(f"{len(missing_genes)} pair-universe loci are absent from profile")
    profile = pd.read_csv(
        args.profile, sep="\t", usecols=[profile_header[0], *genes], index_col=0
    ).fillna(0)
    if profile.index.duplicated().any():
        raise ValueError("Duplicate genome identifiers in profile")
    binary = (profile[genes].to_numpy() > 0).astype(np.float64)
    log.info("Loaded profile subset: %d genomes x %d loci", *binary.shape)

    weight_table = pd.read_csv(args.weights, sep="\t").set_index("genome")
    missing_samples = profile.index.difference(weight_table.index)
    extra_samples = weight_table.index.difference(profile.index)
    if len(missing_samples) or len(extra_samples):
        raise ValueError(
            f"Sample-ID mismatch: {len(missing_samples)} missing and {len(extra_samples)} extra in weights"
        )
    if args.weight_column not in weight_table.columns:
        raise ValueError(f"Missing weight column: {args.weight_column}")
    weights_raw = weight_table.loc[profile.index, args.weight_column].to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(weights_raw)) or np.any(weights_raw <= 0.0):
        raise ValueError("Sample weights must be finite and positive")
    raw_sum = float(weights_raw.sum())
    weights = weights_raw / raw_sum
    log.info(
        "Sample weights: %s, raw sum %.12g, normalized sum %.12g, min %.12g, max %.12g",
        args.weight_column,
        raw_sum,
        weights.sum(),
        weights.min(),
        weights.max(),
    )

    weighted_marginal = weights @ binary
    weighted_joint = (binary * weights[:, None]).T @ binary
    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    idx1 = pairs["gene1"].map(gene_index).to_numpy(dtype=np.int64)
    idx2 = pairs["gene2"].map(gene_index).to_numpy(dtype=np.int64)
    mi, nmi, hx, hy = weighted_binary_information(
        weighted_joint[idx1, idx2],
        weighted_marginal[idx1],
        weighted_marginal[idx2],
        1.0,
    )
    pairs["mi"] = mi
    pairs["nmi"] = nmi
    pairs["gene1_frequency_0"] = 1.0 - weighted_marginal[idx1]
    pairs["gene1_frequency_1"] = weighted_marginal[idx1]
    pairs["gene2_frequency_0"] = 1.0 - weighted_marginal[idx2]
    pairs["gene2_frequency_1"] = weighted_marginal[idx2]
    pairs["frequency_weight"] = 1.0
    pairs["frequency_weighted_MI"] = mi
    pairs["sample_weight_scheme"] = "wg_presence_hamming_tau_0.10_inverse_neighbour_normalized"
    pairs["sample_weight_tau"] = args.tau
    pairs["sample_weight_sum_raw"] = raw_sum
    pairs["sample_weight_sum_normalized"] = float(weights.sum())

    if np.any(mi < -1e-12) or np.any(nmi < -1e-12) or np.any(nmi > 1.0 + 1e-10):
        raise ValueError("MI/NMI numerical range check failed")
    log.info("Calculated MI/NMI; MI range %.12g to %.12g", mi.min(), mi.max())

    old = pd.read_csv(
        args.old_mi,
        sep="\t",
        compression="gzip",
        usecols=["gene1", "gene2", "weighted_mi", "weighted_nmi"],
    )
    if len(old) != len(pairs):
        raise ValueError("Old weighted-MI table row count differs from frozen pair universe")
    compare = pairs[["gene1", "gene2", "mi", "nmi"]].merge(
        old, on=["gene1", "gene2"], how="left", validate="one_to_one"
    )
    if compare[["weighted_mi", "weighted_nmi"]].isna().any().any():
        raise ValueError("Old weighted-MI table does not cover the frozen pair universe")

    correlation = pd.DataFrame(
        [
            {
                "comparison": "mi_v2_vs_old_weighted_mi",
                "n_pairs": len(compare),
                "pearson": compare["mi"].corr(compare["weighted_mi"]),
                "spearman": rank_correlation(compare["mi"], compare["weighted_mi"]),
                "mean_absolute_difference": float(
                    np.mean(np.abs(compare["mi"] - compare["weighted_mi"]))
                ),
            },
            {
                "comparison": "nmi_v2_vs_old_weighted_nmi",
                "n_pairs": len(compare),
                "pearson": compare["nmi"].corr(compare["weighted_nmi"]),
                "spearman": rank_correlation(compare["nmi"], compare["weighted_nmi"]),
                "mean_absolute_difference": float(
                    np.mean(np.abs(compare["nmi"] - compare["weighted_nmi"]))
                ),
            },
        ]
    )

    metrics = {
        "mi": pairs["mi"],
        "nmi": pairs["nmi"],
        "frequency_weight": pairs["frequency_weight"],
        "frequency_weighted_MI": pairs["frequency_weighted_MI"],
    }
    quantiles = [0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.995, 0.999, 1.0]
    summary_rows = []
    for name, values in metrics.items():
        q = values.quantile(quantiles)
        for probability, value in q.items():
            summary_rows.append(
                {
                    "metric": name,
                    "statistic": f"quantile_{probability:g}",
                    "value": float(value),
                }
            )
        summary_rows.extend(
            [
                {"metric": name, "statistic": "mean", "value": float(values.mean())},
                {"metric": name, "statistic": "std", "value": float(values.std())},
                {"metric": name, "statistic": "n", "value": int(values.notna().sum())},
            ]
        )
    distribution = pd.DataFrame(summary_rows)

    full_path = args.outdir / "kp_weighted_mi_v2_all_pairs.tsv.gz"
    pairs.to_csv(full_path, sep="\t", index=False, compression="gzip")
    correlation.to_csv(args.outdir / "old_weighted_mi_correlation.tsv", sep="\t", index=False)
    distribution.to_csv(args.outdir / "distribution_summary.tsv", sep="\t", index=False)
    pairs.nlargest(args.top_n, "mi").to_csv(
        args.outdir / "top_pairs_by_mi.tsv", sep="\t", index=False
    )
    pairs.nlargest(args.top_n, "nmi").to_csv(
        args.outdir / "top_pairs_by_nmi.tsv", sep="\t", index=False
    )
    write_field_dictionary(args.outdir / "field_definitions.md")

    provenance = {
        "method_name": "weighted MI v2",
        "not_epidis": True,
        "profile": str(args.profile),
        "profile_sha256": file_sha256(args.profile),
        "pairs": str(args.pairs),
        "pairs_sha256": file_sha256(args.pairs),
        "weights": str(args.weights),
        "weights_sha256": file_sha256(args.weights),
        "old_mi": str(args.old_mi),
        "old_mi_sha256": file_sha256(args.old_mi),
        "n_genomes": int(binary.shape[0]),
        "n_loci": int(binary.shape[1]),
        "n_pairs": int(len(pairs)),
        "binary_encoding": "profile value > 0 is 1; otherwise 0",
        "sample_weight_column": args.weight_column,
        "sample_weight_tau_hamming_distance": args.tau,
        "sample_weight_raw_sum": raw_sum,
        "sample_weight_normalized_sum": float(weights.sum()),
        "pseudocount": 0.0,
        "mi_log_base": "natural",
        "frequency_weight_decision": "fixed 1.0; alpha/beta already embedded in MI probabilities",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (args.outdir / "run_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    log.info("Wrote full table: %s", full_path)
    log.info("Correlation summary:\n%s", correlation.to_string(index=False))
    log.info("Completed weighted MI v2")


if __name__ == "__main__":
    main()
