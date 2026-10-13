"""Certificates and contract checks for the odd family (F12-F13, V11-V13, V18)."""

import math

import pytest
from scipy.special import spherical_jn

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.spherical_harmonics import normalisation
from mostlyharmless.visylm.visibility import (
    radial_integral,
    visibility_mode,
)

TOL = 1e-9


def odd_modes(l_max):
    """Every (l, m) with l + m odd, m >= 0."""
    for l in range(1, l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                yield l, m


# ------------------------------------------------------------------ V11, V13, V18
def test_v11_odd_collapse_is_exact():
    assert dv.check_odd_collapse_exact(l_max=8) == 0


def test_v18_k_identity():
    assert dv.check_k_identity(l_max=6) < 1e-10


def test_v13_odd_closed_form_against_quadrature():
    assert dv.check_odd_closed_form(l_max=5) < TOL


# ------------------------------------------------------------------ the limits
def test_the_closed_form_limit_at_zero_baseline():
    """j_l(a)/a -> 1/3 for l = 1 and 0 for l >= 2, so the mode visibility at
    q = 0 is 2 pi N_lm times that limit."""
    for l, m in odd_modes(5):
        expected_limit = 1.0 / 3.0 if l == 1 else 0.0
        expected = 2.0 * math.pi * normalisation(l, m) * expected_limit
        assert complex(visibility_mode(l, m, 0.0, 0.0)) == pytest.approx(
            expected, abs=1e-14)


def test_the_lowest_odd_mode_at_zero_carries_one_third():
    # l = 1, m = 0: the one case where the limit is not zero
    expected = 2.0 * math.pi * normalisation(1, 0) / 3.0
    assert complex(visibility_mode(1, 0, 0.0, 0.0)).real == pytest.approx(
        expected)


# ------------------------------------------------------------------ structure
def test_the_odd_modes_match_their_closed_expression():
    """The manual's closed form, evaluated independently of the dispatch."""
    for l, m in odd_modes(5):
        coefficient = ((-1) ** m
                       * math.prod(range(l + m, 0, -2))
                       / math.prod(range(l - m - 1, 0, -2) or [1]))
        for q in (0.0, 0.7, 3.0):
            a = 2.0 * math.pi * q
            if a == 0.0:
                bessel = (1.0 / 3.0) if l == 1 else 0.0
            else:
                bessel = float(spherical_jn(l, a)) / a
            expected = (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
                        * coefficient * bessel)
            assert complex(visibility_mode(l, m, q, 0.0)) == pytest.approx(
                expected, abs=1e-12)


def test_the_even_modes_still_use_the_quadrature():
    # until stage 4 supplies their closed form
    for l in range(1, 5):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                continue
            for q in (0.5, 2.0):
                expected = (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
                            * complex(radial_integral(l, m, q)))
                assert complex(visibility_mode(l, m, q, 0.0)) == pytest.approx(
                    expected, abs=1e-12)


def test_visibility_mode_rejects_a_negative_order():
    with pytest.raises(ValueError):
        visibility_mode(2, -1, 1.0, 0.0)
