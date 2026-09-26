import numpy as np

from cinch_dev2.core.statistics import bh_adjust, directional_epidis, weighted_binary_information


def test_binary_mi_perfect_association():
    mi, nmi, _, _ = weighted_binary_information(
        np.array([0.5]), np.array([0.5]), np.array([0.5])
    )
    assert np.isclose(mi[0], np.log(2.0))
    assert np.isclose(nmi[0], 1.0)


def test_binary_mi_independence():
    mi, nmi, _, _ = weighted_binary_information(
        np.array([0.25]), np.array([0.5]), np.array([0.5])
    )
    assert np.isclose(mi[0], 0.0, atol=1e-15)
    assert np.isclose(nmi[0], 0.0, atol=1e-15)


def test_epidis_identity_and_symmetry():
    c11 = np.array([0.5, 0.25, 0.20])
    x1 = np.array([0.5, 0.5, 0.25])
    y1 = np.array([0.5, 0.5, 0.70])
    forward = directional_epidis(c11, x1, y1)
    reverse = directional_epidis(c11, y1, x1)
    mi, _, _, _ = weighted_binary_information(c11, x1, y1)
    assert np.allclose(forward, reverse, atol=1e-12)
    assert np.allclose(forward**2, mi / np.log(2.0), atol=1e-12)


def test_bh_adjustment():
    observed = bh_adjust(np.array([0.01, 0.04, 0.03, 0.002]))
    expected = np.array([0.02, 0.04, 0.04, 0.008])
    assert np.allclose(observed, expected)


def test_bh_rejects_invalid_values():
    try:
        bh_adjust(np.array([0.2, np.nan]))
    except ValueError:
        pass
    else:
        raise AssertionError("invalid p-values were accepted")

