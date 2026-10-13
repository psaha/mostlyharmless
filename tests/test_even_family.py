"""Certificates and contract checks for the even family (F14-F17, V8-V10, V12, V14)."""

import math
from pathlib import Path

import pytest

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.symbolic import emit
from mostlyharmless.visylm.spherical_harmonics import normalisation
from mostlyharmless.visylm.visibility import (
    radial_integral,
    visibility_mode,
)
from mostlyharmless.visylm.zernike import (
    zernike_coefficients,
    zernike_polynomial,
    zernike_radial,
)

TOL = 1e-9


def even_modes(l_max):
    """Every (l, m) with l + m even, m >= 0."""
    for l in range(l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 0:
                yield l, m


# ------------------------------------------------------------------ V8, V9, V10, V12, V14
def test_v8_jacobi_form_of_the_zernike_polynomials():
    assert dv.check_zernike_jacobi(n_max=8) == 0


def test_v9_zernike_orthogonality():
    assert dv.check_zernike_orthogonality(n_max=8) == 0


def test_v10_zernike_hankel_identity():
    assert dv.check_zernike_hankel(n_max=6, order=6) == 0


def test_v12_projection_coefficients_equal_coefficient_matching():
    assert dv.check_zernike_projection(l_max=6) == 0


def test_v14_even_closed_form_against_quadrature():
    assert dv.check_even_closed_form(l_max=6) < TOL


# ------------------------------------------------------------------ the limits
def test_the_even_closed_form_limit_at_zero_baseline():
    """At q = 0 only the term with m + 2k = 0 survives, so the limit is c_0/2 for
    m = 0 and zero otherwise."""
    for l, m in even_modes(6):
        expected_limit = (float(zernike_coefficients(l, m)[0]) / 2.0
                          if m == 0 else 0.0)
        expected = 2.0 * math.pi * normalisation(l, m) * expected_limit
        # the value at q = 0 is real and signed: c_0 is negative for some l
        assert complex(visibility_mode(l, m, 0.0, 0.0)) == pytest.approx(
            expected, abs=1e-13)


def test_the_uniform_disk_is_the_lowest_even_mode():
    # l = m = 0 is the uniform disk: c_0 = 1 and R_00 = J_1(2 pi q)/(2 pi q)
    assert zernike_coefficients(0, 0) == (1,)
    for q in (0.0, 0.5, 3.0):
        from scipy.special import jv
        a = 2.0 * math.pi * q
        expected = (2.0 * math.pi * normalisation(0, 0)
                    * (0.5 if a == 0.0 else float(jv(1, a)) / a))
        assert complex(visibility_mode(0, 0, q, 0.0)).real == pytest.approx(
            expected, abs=1e-12)


# ------------------------------------------------------------------ contracts
def test_the_jacobi_evaluation_matches_the_polynomial():
    for n, m in ((0, 0), (2, 0), (3, 1), (4, 2), (6, 2)):
        expression = zernike_polynomial(n, m)
        # (0, 0) is the constant 1, so there may be no free symbol to substitute
        symbols = sorted(expression.free_symbols, key=str)
        for value in (0.0, 0.25, 0.5, 1.0):
            substituted = (expression.subs(symbols[0], value) if symbols
                           else expression)
            assert float(zernike_radial(n, m, value)) == pytest.approx(
                float(substituted), abs=1e-12)


def test_the_coefficients_are_exact_rationals():
    from fractions import Fraction
    for l, m in even_modes(6):
        for value in zernike_coefficients(l, m):
            assert Fraction(str(value)) == Fraction(value)


def test_the_even_family_rejects_odd_l_plus_m():
    with pytest.raises(ValueError):
        zernike_coefficients(2, 1)
    with pytest.raises(ValueError):
        zernike_coefficients(1, 0)


def test_the_even_modes_no_longer_need_the_quadrature():
    """After stage 4 every mode has a closed form, so the quadrature agrees with
    the dispatch but is no longer reached by it."""
    for l, m in even_modes(5):
        for q in (0.4, 2.0):
            expected = (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
                        * complex(radial_integral(l, m, q)))
            assert complex(visibility_mode(l, m, q, 0.0)) == pytest.approx(
                expected, abs=1e-12)


# ------------------------------------------------------------------ the emitter
def test_the_generated_coefficient_table_is_current():
    path = (Path(__file__).resolve().parents[1] / "docs" / "generated"
            / "coefficients.tex")
    if not path.exists():
        pytest.skip(f"{path} not generated yet; run `make data`")
    assert path.read_text() == emit.coefficients_table()
