"""Exact-equivalent statistical primitives extracted from the local Cinch workflow.

These functions deliberately preserve the frozen definitions. Engineering changes
must be tested against them before replacing any production kernel.
"""

from __future__ import annotations

import numpy as np


EPSILON = 1e-27


def weighted_binary_information(
    joint_11: np.ndarray,
    marginal_x1: np.ndarray,
    marginal_y1: np.ndarray,
    total_weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return weighted MI in nats, symmetric NMI, H(X), and H(Y).

    Inputs are weighted masses rather than raw sample arrays. No pseudocount is
    applied. This is the frozen weighted-MI-v2 definition.
    """

    c11 = np.asarray(joint_11, dtype=np.float64)
    m1 = np.asarray(marginal_x1, dtype=np.float64)
    m2 = np.asarray(marginal_y1, dtype=np.float64)
    counts = np.column_stack(
        [total_weight - m1 - m2 + c11, m2 - c11, m1 - c11, c11]
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


def _xlog2_ratio(x: np.ndarray, denominator: np.ndarray) -> np.ndarray:
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
    """Evaluate the frozen EpiDis equations in one conditioning direction."""

    joint_11 = np.asarray(joint_11, dtype=np.float64)
    alpha = np.asarray(conditioning_1, dtype=np.float64)
    target_1 = np.asarray(target_1, dtype=np.float64)
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
    jsd = a * _xlog2_ratio(p, mixture).sum(axis=1) + b * _xlog2_ratio(q, mixture).sum(axis=1)
    result[valid] = np.sqrt(np.maximum(jsd, 0.0))
    return result


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment matching the frozen SHC implementation."""

    pvalues = np.asarray(pvalues, dtype=np.float64)
    if pvalues.ndim != 1:
        raise ValueError("pvalues must be one-dimensional")
    if np.any(~np.isfinite(pvalues)) or np.any((pvalues < 0.0) | (pvalues > 1.0)):
        raise ValueError("pvalues must be finite values in [0, 1]")
    if pvalues.size == 0:
        return pvalues.copy()
    order = np.argsort(pvalues, kind="mergesort")
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output

