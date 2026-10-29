"""Certificates and implementation checks for the spherical harmonics (F1-F3).

Tags, rungs and scopes follow ROADMAP.md section 4 and decision D10.  A test
named after a tag is a certificate; the rest are contract or anchor checks.
"""

import math

import numpy as np
import pytest
from scipy.special import sph_harm_y

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.spherical_harmonics import (
    associated_legendre,
    factorial_ratio,
    m_values,
    normalisation,
    sph_harm,
)

TOL = 1e-11
HIGH_PRECISION_TOL = 1e-25


# ---------------------------------------------------------------- V1, exact
def test_v1_legendre_orthogonality_exact():
    assert dv.check_legendre_orthogonality_exact(l_max=5) == 0


# ------------------------------------------------------------ V1, numerical
def test_v1_matches_sympy_ynm():
    assert dv.check_ylm_matches_sympy_ynm(l_max=3) < HIGH_PRECISION_TOL


def test_v1_orthonormality_of_the_implementation():
    assert dv.check_orthonormality_numeric(l_max=4) < TOL


# ---------------------------------------------------------------- V2, exact
def test_v2_reality_exact():
    assert dv.check_reality_exact(l_max=4) == 0


def test_v2_reality_of_the_implementation():
    for l in range(5):
        for m in range(1, l + 1):
            for theta, phi in ((0.4, 0.2), (1.9, -1.1), (2.8, 3.0)):
                lhs = sph_harm(l, -m, theta, phi)
                rhs = (-1) ** m * np.conjugate(sph_harm(l, m, theta, phi))
                assert lhs == pytest.approx(rhs, abs=1e-12)


# ------------------------------------------------------------ V3, numerical
def test_v3_addition_theorem_complex_form():
    assert dv.check_addition_theorem_numeric(l_max=5) < TOL


def test_v3_addition_theorem_real_form():
    # the form the manual states; every index non-negative
    assert dv.check_addition_theorem_real_numeric(l_max=5) < TOL


# ------------------------------------------------------------------ anchors
def test_y00_is_the_constant():
    assert sph_harm(0, 0, 0.7, 1.3) == pytest.approx(1 + 0j)   # Y_00 = 1 in this convention


def test_y10_is_the_z_harmonic():
    theta = 0.9
    assert sph_harm(1, 0, theta, 2.0) == pytest.approx(
        math.sqrt(3.0) * math.cos(theta) + 0j)


def test_matches_scipy_sph_harm_y():
    # scipy 1.15 renamed this and made theta the polar angle, which is our
    # convention.  The old scipy.special.sph_harm had theta as the azimuth.
    for l in range(4):
        for m in range(-l, l + 1):
            for theta in (0.3, 1.1, 2.4):
                for phi in (0.0, 1.7, -2.2):
                    # scipy is orthonormal; ours carries sqrt(4 pi) relative to it
                    assert sph_harm(l, m, theta, phi) == pytest.approx(
                        math.sqrt(4 * math.pi) * sph_harm_y(l, m, theta, phi),
                        abs=1e-12)


def test_negative_order_relation_for_legendre():
    for l in range(4):
        for m in range(1, l + 1):
            for x in (-0.7, 0.0, 0.35):
                lhs = associated_legendre(l, -m, x)
                rhs = ((-1) ** m * factorial_ratio(l, m)
                       * associated_legendre(l, m, x))
                assert lhs == pytest.approx(rhs, abs=1e-12)


def test_real_basis_is_orthonormal():
    assert dv.check_real_basis_orthonormality_numeric(l_max=4) < TOL


def test_real_basis_matches_the_documented_convention():
    assert dv.check_real_basis_matches_complex(l_max=4) < TOL


def test_m_values_are_ascending_from_minus_l():
    assert list(m_values(3)) == [-3, -2, -1, 0, 1, 2, 3]


def test_normalisation_uses_the_signed_m():
    # The trap this guards: with |m| the negative orders come out unnormalised.
    assert normalisation(1, -1) == pytest.approx(2.0 * normalisation(1, 1))


def test_orders_beyond_l_are_rejected():
    with pytest.raises(ValueError):
        factorial_ratio(2, 3)
    with pytest.raises(ValueError):
        associated_legendre(2, 3, 0.5)
