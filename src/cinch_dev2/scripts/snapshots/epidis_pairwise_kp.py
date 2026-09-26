#!/usr/bin/env python3
"""Compute thesis-defined pairwise EpiDis on the frozen Kp pair universe."""

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


EPSILON = 1e-27


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def xlog2_ratio(x: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    valid = (x > 0.0) & (denominator > 0.0)
    out[valid] = x[valid] * np.log2(x[valid] / denominator[valid])
    return out


def directional_epidis(
    joint_11: np.ndarray,
    conditioning_1: np.ndarray,
    target_1: np.ndarray,
    epsilon: float = EPSILON,
) -> np.ndarray:
    """Implement thesis equations 2-8 through 2-14 literally, without post-epsilon renormalization."""
    alpha = conditioning_1
    beta = 1.0 - alpha
    valid = (alpha > 0.0) & (beta > 0.0)
    result = np.full(alpha.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result

    a = alpha[valid]
    b = beta[valid]
    c11 = joint_11[valid]
    t1 = target_1[valid]
    p1 = np.clip(c11 / a, 0.0, 1.0)
    q1 = np.clip((t1 - c11) / b, 0.0, 1.0)
    p = np.column_stack((p1 + epsilon, 1.0 - p1 + epsilon))
    q = np.column_stack((q1 + epsilon, 1.0 - q1 + epsilon))
    mixture = a[:, None] * p + b[:, None] * q
    jsd = a * xlog2_ratio(p, mixture).sum(axis=1) + b * xlog2_ratio(q, mixture).sum(axis=1)
    jsd = np.maximum(jsd, 0.0)
    result[valid] = np.sqrt(jsd)
    return result


def rank_correlation(x: pd.Series, y: pd.Series) -> float:
    return float(x.rank(method="average").corr(y.rank(method="average")))


def correlations(frame: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("epidis_vs_mi_v2", "epidis", "mi_v2"),
        ("epidis_squared_vs_mi_v2_bits", "epidis_squared", "mi_v2_bits"),
        ("epidis_vs_old_weighted_mi", "epidis", "old_weighted_mi"),
        ("epidis_vs_old_weighted_nmi", "epidis", "old_weighted_nmi"),
    ]
    rows = []
    for name, left, right in comparisons:
        valid = frame[[left, right]].dropna()
        rows.append(
            {
                "comparison": name,
                "n_pairs": len(valid),
                "pearson": valid[left].corr(valid[right]),
                "spearman": rank_correlation(valid[left], valid[right]),
                "max_absolute_difference": float(np.max(np.abs(valid[left] - valid[right]))),
            }
        )
    return pd.DataFrame(rows)


def validate_formula() -> pd.DataFrame:
    cases = pd.DataFrame(
        [
            ("perfect_association", 0.5, 0.5, 0.5, 1.0),
            ("independence", 0.25, 0.5, 0.5, 0.0),
            ("asymmetric_marginals", 0.20, 0.25, 0.70, np.nan),
            ("constant_target", 0.0, 0.4, 0.0, 0.0),
        ],
        columns=["case", "joint_11", "x1", "y1", "expected_epidis"],
    )
    xy = directional_epidis(cases.joint_11.to_numpy(), cases.x1.to_numpy(), cases.y1.to_numpy())
    yx = directional_epidis(cases.joint_11.to_numpy(), cases.y1.to_numpy(), cases.x1.to_numpy())
    cases["epidis_x_to_y"] = xy
    cases["epidis_y_to_x"] = yx
    cases["direction_absolute_difference"] = np.abs(xy - yx)
    checkable = cases.expected_epidis.notna()
    cases["expected_pass"] = True
    cases.loc[checkable, "expected_pass"] = np.isclose(
        cases.loc[checkable, "epidis_x_to_y"], cases.loc[checkable, "expected_epidis"], atol=1e-12
    )
    if not cases.expected_pass.all() or cases.direction_absolute_difference.max() > 1e-12:
        raise AssertionError("Synthetic EpiDis validation failed")
    return cases


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=Path("Leca/source/profile/0427_1.profile"))
    parser.add_argument("--pairs", type=Path, default=Path("Leca/source/profile/0427_1.ldlike_pairs.tsv.gz"))
    parser.add_argument("--weights", type=Path, default=Path("Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv"))
    parser.add_argument("--mi-v2", type=Path, default=Path("Leca/results/gene_gene_wgmlst/weighted_mi_v2/kp_weighted_mi_v2_all_pairs.tsv.gz"))
    parser.add_argument("--old-mi", type=Path, default=Path("Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"))
    parser.add_argument("--outdir", type=Path, default=Path("Leca/results/gene_gene_wgmlst/epidis_pairwise_kp"))
    parser.add_argument("--weight-column", default="wg_content_hamming_weight")
    parser.add_argument("--tau", type=float, default=0.10)
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s\t%(levelname)s\t%(message)s", handlers=[logging.FileHandler(args.outdir / "run.log", mode="w"), logging.StreamHandler(sys.stdout)])
    log = logging.getLogger("epidis")
    log.info("Starting thesis-literal pairwise EpiDis")
    validation = validate_formula()
    validation.to_csv(args.outdir / "small_scale_validation.tsv", sep="\t", index=False)

    pairs = pd.read_csv(args.pairs, sep="\t", compression="gzip")
    genes = sorted(set(pairs.gene1) | set(pairs.gene2))
    if len(pairs) != 2_609_249 or len(genes) != 3_279:
        raise ValueError("Frozen Kp pair universe mismatch")
    if pairs.duplicated(["gene1", "gene2"]).any():
        raise ValueError("Duplicate pair keys")

    header = pd.read_csv(args.profile, sep="\t", nrows=0).columns.tolist()
    if not set(genes).issubset(header[1:]):
        raise ValueError("Pair loci missing from profile")
    profile = pd.read_csv(args.profile, sep="\t", usecols=[header[0], *genes], index_col=0).fillna(0)
    binary = (profile[genes].to_numpy() > 0).astype(np.float64)
    weight_table = pd.read_csv(args.weights, sep="\t").set_index("genome")
    if len(profile.index.difference(weight_table.index)) or len(weight_table.index.difference(profile.index)):
        raise ValueError("Profile and weight sample IDs differ")
    raw_weights = weight_table.loc[profile.index, args.weight_column].to_numpy(dtype=float)
    weights = raw_weights / raw_weights.sum()
    marginal = weights @ binary
    joint = (binary * weights[:, None]).T @ binary
    index = {gene: i for i, gene in enumerate(genes)}
    i = pairs.gene1.map(index).to_numpy(dtype=np.int64)
    j = pairs.gene2.map(index).to_numpy(dtype=np.int64)
    c11 = joint[i, j]
    forward = directional_epidis(c11, marginal[i], marginal[j])
    reverse = directional_epidis(c11, marginal[j], marginal[i])
    direction_diff = np.abs(forward - reverse)
    max_diff_index = int(np.nanargmax(direction_diff))
    log.info(
        "Directional audit: max abs diff %.12g at %s / %s (%.12g vs %.12g)",
        direction_diff[max_diff_index],
        pairs.iloc[max_diff_index].gene1,
        pairs.iloc[max_diff_index].gene2,
        forward[max_diff_index],
        reverse[max_diff_index],
    )
    if np.nanmax(direction_diff) > 1e-9:
        raise ValueError("Directional EpiDis discrepancy exceeds tolerance")
    pairs["epidis_i_to_j"] = forward
    pairs["epidis_j_to_i"] = reverse
    pairs["epidis_direction_abs_diff"] = direction_diff
    pairs["epidis"] = (forward + reverse) / 2.0
    pairs["epidis_squared"] = pairs.epidis**2
    pairs["gene1_weighted_presence"] = marginal[i]
    pairs["gene2_weighted_presence"] = marginal[j]

    mi_v2 = pd.read_csv(args.mi_v2, sep="\t", compression="gzip", usecols=["gene1", "gene2", "mi"] ).rename(columns={"mi": "mi_v2"})
    old = pd.read_csv(args.old_mi, sep="\t", compression="gzip", usecols=["gene1", "gene2", "weighted_mi", "weighted_nmi"]).rename(columns={"weighted_mi": "old_weighted_mi", "weighted_nmi": "old_weighted_nmi"})
    pairs = pairs.merge(mi_v2, on=["gene1", "gene2"], validate="one_to_one").merge(old, on=["gene1", "gene2"], validate="one_to_one")
    pairs["mi_v2_bits"] = pairs.mi_v2 / np.log(2.0)
    valid_scores = pairs.epidis.dropna()
    q1, q3 = valid_scores.quantile([0.25, 0.75])
    iqr = q3 - q1
    threshold = q3 + 3.0 * iqr
    pairs["epidis_iqr_outlier"] = pairs.epidis > threshold

    corr = correlations(pairs)
    corr.to_csv(args.outdir / "correlations_with_mi.tsv", sep="\t", index=False)
    threshold_table = pd.DataFrame([{"n_valid_pairs": len(valid_scores), "q1": q1, "q3": q3, "iqr": iqr, "k": 2, "iqr_multiplier": 3.0, "threshold": threshold, "n_outliers": int(pairs.epidis_iqr_outlier.sum()), "outlier_fraction": float(pairs.epidis_iqr_outlier.mean())}])
    threshold_table.to_csv(args.outdir / "epidis_iqr_threshold.tsv", sep="\t", index=False)
    pairs.nlargest(args.top_n, "epidis").to_csv(args.outdir / "top_epidis_pairs.tsv", sep="\t", index=False)
    pairs.to_csv(args.outdir / "kp_pairwise_epidis_all_pairs.tsv.gz", sep="\t", index=False, compression="gzip")

    summary = {
        "n_genomes": int(binary.shape[0]), "n_loci": len(genes), "n_pairs": len(pairs),
        "n_valid_epidis": int(valid_scores.size), "n_undefined_epidis": int(pairs.epidis.isna().sum()),
        "epidis_min": float(valid_scores.min()), "epidis_median": float(valid_scores.median()),
        "epidis_mean": float(valid_scores.mean()), "epidis_max": float(valid_scores.max()),
        "max_direction_absolute_difference": float(np.nanmax(direction_diff)),
        "iqr_threshold": float(threshold), "n_iqr_outliers": int(pairs.epidis_iqr_outlier.sum()),
        "iqr_outlier_fraction": float(pairs.epidis_iqr_outlier.mean()),
    }
    pd.DataFrame([summary]).to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    provenance = {
        "method": "thesis-literal EpiDis equations 2-8 through 2-14",
        "profile": str(args.profile), "profile_sha256": sha256(args.profile),
        "pairs": str(args.pairs), "pairs_sha256": sha256(args.pairs),
        "weights": str(args.weights), "weights_sha256": sha256(args.weights),
        "binary_encoding": "profile value > 0 is presence 1; all other values are absence/non-call 0",
        "pair_universe": "frozen 3,279 loci and 2,609,249 pairs",
        "sample_weight_column": args.weight_column, "sample_weight_tau": args.tau,
        "sample_weight_raw_sum": float(raw_weights.sum()), "sample_weight_normalized_sum": float(weights.sum()),
        "epsilon": EPSILON, "epsilon_post_addition_renormalized": False,
        "log_base": 2, "direction_rule": "compute both directions; report mean after tolerance validation",
        "threshold": "all valid pair scores; Q3 + 2*1.5*IQR = Q3 + 3*IQR",
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
    }
    (args.outdir / "run_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (args.outdir / "implementation_notes.md").write_text(
        "# Kp pairwise EpiDis implementation notes\n\n"
        "The score directly implements thesis equations 2-8 through 2-14. It does not call MI or JSD library functions. "
        "The public repository has no pairwise EpiDis implementation, so there is no code path to reproduce.\n\n"
        "Kp decisions: presence is profile value >0; all other values are collapsed to 0; thesis-style Hamming-neighbour "
        "weights at tau=0.10 are normalized to sum one; epsilon=1e-27 is added literally without post-addition "
        "renormalization; both directions are retained and their mean is reported only after numerical equality testing; "
        "the comparison universe is the existing 3,279 loci and 2,609,249 pairs.\n\n"
        "Differences from source code: the audited GitHub repository only implements major-state binary encoding and does not "
        "implement weights, EpiDis, thresholds, or pairwise scanning. Differences from the thesis are Kp-specific handling of "
        "non-positive calls as absence/non-call, the selected tau=0.10 within the thesis range, and explicit skipping of a "
        "zero-mass conditioning background.\n"
    )
    log.info("Completed: %s", json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
