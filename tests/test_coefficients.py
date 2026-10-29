"""Tests for the Coeffs container.

Not certificates: these are contract checks.  The point of the container is that
a real map is real *by construction*, so most of what is tested here is that the
construction holds, and that the two representations of the same map agree.
"""

import math

import numpy as np
import pytest

from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.spherical_harmonics import real_sph_harm, sph_harm

LMAX = 3
ANGLES = ((0.4, 0.2), (1.7, -1.1), (2.6, 2.9))


def random_coeffs(lmax=LMAX, seed=7):
    rng = np.random.default_rng(seed)
    return Coeffs.from_vector(rng.normal(size=(lmax + 1) ** 2))


def assert_close(a, b, tol=1e-14):
    """Compare two Coeffs allowing for the rounding in a sqrt(2) round trip."""
    assert a.lmax == b.lmax
    for block_a, block_b in zip(a.blocks, b.blocks):
        assert np.allclose(block_a, block_b, rtol=0, atol=tol)


# ------------------------------------------------------------------ basics
def test_zeros_is_a_zero_map():
    c = Coeffs.zeros(2)
    assert c.lmax == 2
    assert all(np.all(block == 0) for block in c.blocks)
    assert np.all(c.to_vector() == 0)


def test_vector_length_is_lmax_plus_one_squared():
    for lmax in range(4):
        assert Coeffs.zeros(lmax).to_vector().size == (lmax + 1) ** 2


def test_blocks_have_the_stored_half_length():
    c = random_coeffs()
    assert [len(b) for b in c.blocks] == [l + 1 for l in range(LMAX + 1)]


def test_accessor_derives_the_negative_orders():
    c = random_coeffs()
    for l in range(c.lmax + 1):
        for m in range(1, l + 1):
            assert c[l, -m] == pytest.approx((-1) ** m * np.conjugate(c[l, m]))


def test_accessor_rejects_orders_beyond_l():
    c = random_coeffs()
    with pytest.raises(ValueError):
        c[2, 3]


# ------------------------------------------------------- the map is real
def test_reconstructed_map_is_real_and_matches_the_real_basis():
    """The two representations of the same map must agree, and be real."""
    c = random_coeffs()
    real = c.to_real()
    for theta, phi in ANGLES:
        from_complex = 0.0
        for l in range(c.lmax + 1):
            for m in range(-l, l + 1):
                from_complex += c[l, m] * sph_harm(l, m, theta, phi)
        from_real = 0.0
        for l in range(c.lmax + 1):
            for m in range(-l, l + 1):
                from_real += real[l][l + m] * real_sph_harm(l, m, theta, phi)

        assert abs(from_complex.imag) < 1e-12
        assert from_complex.real == pytest.approx(from_real, abs=1e-12)


def test_the_worked_example_from_the_discussion():
    """A pure cos(phi) pattern at l = 1.

    Real coefficients p_1 = 1, q_-1 = 0 must give a_1 = +1/sqrt(2) and
    a_-1 = -1/sqrt(2), and describe the real map -sqrt(3/4pi) sin(theta) cos(phi).
    """
    real = [np.zeros(1), np.array([0.0, 0.0, 1.0])]  # index m + l
    c = Coeffs.from_real(real)

    assert c[1, 1] == pytest.approx(1 / math.sqrt(2))
    assert c[1, -1] == pytest.approx(-1 / math.sqrt(2))

    for theta, phi in ANGLES:
        expected = -math.sqrt(3.0) * math.sin(theta) * math.cos(phi)
        got = (c[1, 1] * sph_harm(1, 1, theta, phi)
               + c[1, -1] * sph_harm(1, -1, theta, phi))
        assert got.imag == pytest.approx(0.0, abs=1e-14)
        assert got.real == pytest.approx(expected, abs=1e-12)


# ------------------------------------------------------------- round trips
def test_real_round_trip():
    c = random_coeffs()
    assert_close(Coeffs.from_real(c.to_real()), c)


def test_vector_round_trip():
    c = random_coeffs()
    assert_close(Coeffs.from_vector(c.to_vector()), c)


def test_unfolded_round_trip():
    # exact: unfolding builds the negative orders from the same formula the
    # fold checks against, so nothing is rounded on the way
    c = random_coeffs()
    assert Coeffs.from_unfolded(c.unfolded()) == c


def test_unfolded_blocks_have_the_full_length():
    c = random_coeffs()
    assert [len(b) for b in c.unfolded()] == [2 * l + 1 for l in range(LMAX + 1)]


# -------------------------------------------------------------- validation
def test_from_unfolded_catches_a_reality_violation():
    c = random_coeffs()
    blocks = [b.copy() for b in c.unfolded()]
    blocks[2][0] += 0.5  # corrupt the order -2 entry of degree 2
    with pytest.raises(ValueError, match="reality condition"):
        Coeffs.from_unfolded(blocks)


def test_wrong_block_lengths_are_rejected():
    with pytest.raises(ValueError, match="expected"):
        Coeffs((np.zeros(3),))  # degree 0 must have 1 entry
    with pytest.raises(ValueError, match="expected"):
        Coeffs.from_real((np.zeros(1), np.zeros(2)))  # degree 1 needs 3
    with pytest.raises(ValueError, match=r"\(lmax\+1\)\^2"):
        Coeffs.from_vector(np.zeros(5))


def test_a_complex_a_l0_is_rejected():
    with pytest.raises(ValueError, match="a_l0 must be real"):
        Coeffs((np.zeros(1), np.array([0.3 + 0.4j, 0.0])))


def test_blocks_are_read_only():
    c = random_coeffs()
    with pytest.raises(ValueError):
        c.blocks[0][0] = 1.0


def test_the_flat_vector_puts_a_l0_at_l_squared_plus_l():
    """A degree's row has length 2l+1 with (l, 0) in the middle, so the flat index
    of a_l0 is l*l + l, not l*l.

    Writing a zonal map at l*l puts its coefficients in the sin(l*phi) slots, which
    drew a banana where a zonal pattern was meant, and that is how this was noticed.
    The certificates do not go through the flat vector --- it is the fitting interface
    --- so nothing had pinned it.
    """
    vector = np.arange(1.0, 10.0)          # (lmax + 1)^2 = 9, so lmax = 2
    coeffs = Coeffs.from_vector(vector)
    for l in range(3):
        assert coeffs[l, 0] == pytest.approx(vector[l * l + l])
    assert coeffs.to_vector() == pytest.approx(vector)
