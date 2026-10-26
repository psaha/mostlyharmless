"""Certificates and contract checks for the rotation (F4-F8, V4-V7)."""

import math

import numpy as np
import pytest

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.wigner import (
    rotate_coefficients,
    rotation_matrix,
    wigner_d,
    wigner_D,
)

TOL = 1e-11
LMAX = 3


def random_coeffs(lmax=LMAX, seed=11):
    rng = np.random.default_rng(seed)
    return Coeffs.from_vector(rng.normal(size=(lmax + 1) ** 2))


def full_range_norm(c):
    """|a_0|^2 + 2 sum_{m>0} |a_m|^2, which a unitary rotation preserves."""
    total = 0.0
    for l in range(c.lmax + 1):
        block = c.blocks[l]
        total += abs(block[0]) ** 2 + 2.0 * sum(
            abs(block[m]) ** 2 for m in range(1, l + 1))
    return total


# ------------------------------------------------------------------ V4 to V7
def test_v4_d_matrix_against_the_exponential():
    assert dv.check_wigner_d_against_exponential(l_max=6) < TOL


def test_v5_d_matrix_properties():
    assert dv.check_wigner_d_properties(l_max=6) < TOL


def test_v6_rotation_against_geometry():
    assert dv.check_rotation_against_geometry(l_max=4) < TOL


def test_v6_rotation_nonnegative_form():
    # the transform as the manual states it, no negative-order coefficients
    assert dv.check_rotation_nonnegative_form(l_max=4) < TOL


def test_v7_rotation_composition():
    assert dv.check_rotation_composition(l_max=3) < TOL


# ------------------------------------------------------------------ contracts
def test_d_matrix_is_real():
    for l in range(4):
        assert np.allclose(wigner_d(l, 0.7).imag, 0.0)


def test_D_matrix_is_unitary():
    for l in range(4):
        D = wigner_D(l, 0.3, 1.2, -0.7)
        assert np.allclose(D.conj().T @ D, np.eye(2 * l + 1), atol=TOL)


def test_rotation_matrix_is_a_proper_rotation():
    R = rotation_matrix(0.4, 1.1, 2.3)
    assert np.allclose(R.T @ R, np.eye(3), atol=TOL)
    assert np.isclose(np.linalg.det(R), 1.0, atol=TOL)


def test_identity_viewing_angles():
    # alpha = 0, beta = 0, PA = 0 gives R = R_z(0) R_y(0) R_z(0) = I
    assert np.allclose(rotation_matrix(0.0, 0.0, 0.0), np.eye(3), atol=TOL)

    c = random_coeffs()
    rotated = rotate_coefficients(c, 0.0, 0.0, 0.0)
    for l in range(c.lmax + 1):
        assert np.allclose(rotated.blocks[l], c.blocks[l], atol=TOL)


def test_rotation_preserves_the_full_range_norm():
    c = random_coeffs()
    for alpha, beta, pa in ((0.3, 0.9, 1.4), (2.1, 0.4, -0.6)):
        assert full_range_norm(rotate_coefficients(c, alpha, beta, pa)) == (
            pytest.approx(full_range_norm(c), abs=1e-12))


def test_rotation_returns_coeffs_of_the_same_degree():
    c = random_coeffs()
    r = rotate_coefficients(c, 0.3, 0.9, 1.4)
    assert isinstance(r, Coeffs)
    assert r.lmax == c.lmax
    assert r.to_vector().size == c.to_vector().size


def test_rotation_of_the_zero_map_is_the_zero_map():
    r = rotate_coefficients(Coeffs.zeros(3), 0.7, 1.3, 2.9)
    assert np.all(r.to_vector() == 0)


def test_rotation_leaves_the_uniform_component_alone():
    # D16: maps are normalised to a_00 = 1, and the l = 0 block is invariant
    # under any rotation, since d at l = 0 is the identity.
    rng = np.random.default_rng(5)
    vector = rng.normal(size=16)
    vector[0] = 1.0  # a_00; in the (l+1)^2 packing this is the first entry
    c = Coeffs.from_vector(vector)
    assert c[0, 0] == pytest.approx(1.0)
    for alpha, beta, pa in ((0.7, 1.1, 2.4), (3.0, 0.2, -0.5)):
        assert rotate_coefficients(c, alpha, beta, pa)[0, 0] == pytest.approx(
            1.0, abs=1e-14)


def test_the_cached_rotation_matches_a_fresh_d_matrix():
    """The d-matrix cache must not change a rotation, only its cost.

    It is keyed on beta alone, so the risk is a stale matrix surviving a change of
    inclination --- which would rotate the map by the previous angle.
    """
    import numpy as np

    from mostlyharmless.visylm.coefficients import Coeffs
    from mostlyharmless.visylm import wigner

    rng = np.random.default_rng(3)
    blocks = []
    for l in range(4):
        block = rng.standard_normal(l + 1) + 1j * rng.standard_normal(l + 1)
        block[0] = block[0].real
        blocks.append(block)
    coefficients = Coeffs(tuple(blocks))

    # two rotations at different inclinations, then back to the first: the second call at
    # the first angle must not reuse the matrix built for the second
    first = wigner.rotate_coefficients(coefficients, 0.3, 0.7, 0.1)
    second = wigner.rotate_coefficients(coefficients, 0.3, 1.2, 0.1)
    again = wigner.rotate_coefficients(coefficients, 0.3, 0.7, 0.1)

    print(f"\nrepeated rotation agrees to {np.max(np.abs(first.to_real()[3] - again.to_real()[3])):.1e}")
    assert np.allclose(first.to_real()[3], again.to_real()[3], atol=1e-12)
    assert not np.allclose(first.to_real()[3], second.to_real()[3], atol=1e-3)
