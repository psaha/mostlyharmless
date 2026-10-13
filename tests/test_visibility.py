"""Certificates and contract checks for the master reduction (F9-F11, V16-V17)."""

import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import jv

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.spherical_harmonics import (
    associated_legendre,
    normalisation,
)
from mostlyharmless.visylm.visibility import (
    radial_integral,
    visibility,
    visibility_mode,
    visibility_uniform_disk,
)

TOL = 1e-9


def uniform_coeffs(l_max=0):
    """The uniform map: a_00 = 1 and nothing else."""
    vector = np.zeros((l_max + 1) ** 2)
    vector[0] = 1.0
    return Coeffs.from_vector(vector)


# ------------------------------------------------------------------ V16, V17
def test_v17_jacobi_anger_factors():
    assert dv.check_jacobi_anger() < 1e-11


def test_v16_uniform_disk():
    assert dv.check_uniform_disk_visibility() < TOL


# ------------------------------------------------------------------ contracts
def test_radial_integral_matches_an_independent_quadrature():
    # scipy's adaptive quad, which shares nothing with Gauss-Legendre
    for l in range(4):
        for m in range(l + 1):
            for q in (0.0, 0.7, 2.0):
                expected, _ = quad(
                    lambda rho: (associated_legendre(l, m, math.sqrt(1 - rho ** 2))
                                 * jv(m, 2.0 * math.pi * q * rho) * rho),
                    0.0, 1.0, limit=200)
                assert radial_integral(l, m, q) == pytest.approx(expected, abs=1e-10)


def test_radial_integral_vanishes_for_m_positive_at_the_origin():
    # J_m(0) = 0 for m >= 1, and there is no division by q to rescue it
    for l in range(1, 5):
        for m in range(1, l + 1):
            assert radial_integral(l, m, 0.0) == pytest.approx(0.0, abs=1e-14)


def test_contribution_matches_equation_35():
    """A map holding one stored coefficient a_lm contributes

        4 pi (-i)^m N_lm R_lm(q) Re(a_lm exp(i m psi)) .

    This is the statement that the partner term, with V_l^{-m} = conj(V_l^m),
    has been folded in correctly.
    """
    rng = np.random.default_rng(21)
    for l in range(1, 4):
        for m in range(1, l + 1):
            a = complex(rng.normal(), rng.normal())
            blocks = [np.zeros(k + 1, dtype=complex) for k in range(l + 1)]
            blocks[l][m] = a
            coeffs = Coeffs(tuple(blocks))
            for psi in (0.4, 1.3, 2.7):
                for q in (0.0, 1.5):
                    u, v = q * math.cos(psi), q * math.sin(psi)
                    expected = (4.0 * math.pi * (-1j) ** m * normalisation(l, m)
                                * radial_integral(l, m, q)
                                * np.real(a * np.exp(1j * m * psi)))
                    assert visibility(coeffs, u, v) == pytest.approx(
                        expected, abs=1e-10)


def test_a_single_mode_map_reproduces_the_mode_visibility():
    blocks = [np.zeros(k + 1, dtype=complex) for k in range(3)]
    blocks[2][1] = 2.0 - 0.5j
    coeffs = Coeffs(tuple(blocks))
    for u, v in ((0.7, 0.3), (2.0, 0.0)):
        expected = (blocks[2][1] * visibility_mode(2, 1, u, v)
                    + (-1) ** 1 * np.conjugate(blocks[2][1])
                    * np.conjugate(visibility_mode(2, 1, u, v)))
        assert visibility(coeffs, u, v) == pytest.approx(expected, abs=1e-12)


def test_visibility_is_linear_in_the_coefficients():
    rng = np.random.default_rng(3)
    coeffs = Coeffs.from_vector(rng.normal(size=9))
    for scale in (2.0, -0.5):
        scaled = Coeffs.from_vector(scale * coeffs.to_vector())
        assert visibility(scaled, 1.2, 0.7) == pytest.approx(
            scale * visibility(coeffs, 1.2, 0.7), abs=1e-12)


def test_uniform_disk_anchor_is_one_at_the_origin_and_depends_only_on_q():
    assert visibility_uniform_disk(0.0, 0.0, 1.0) == pytest.approx(1.0)
    for theta in (0.3, 1.1, 2.7):
        u, v = 1.7 * math.cos(theta), 1.7 * math.sin(theta)
        assert visibility_uniform_disk(u, v, 1.0) == pytest.approx(
            visibility_uniform_disk(1.7, 0.0, 1.0), abs=1e-12)


def test_negative_orders_are_rejected():
    with pytest.raises(ValueError):
        visibility_mode(2, -1, 1.0, 0.0)
    with pytest.raises(ValueError):
        radial_integral(2, -1, 1.0)
    with pytest.raises(ValueError):
        visibility_mode(2, 3, 1.0, 0.0)


def test_the_uniform_map_flux_is_the_analytic_value():
    # a_00 = 1 gives mean intensity 1/sqrt(4 pi) and flux sqrt(pi)/2
    assert float(np.real(visibility(uniform_coeffs(0), 0.0, 0.0))) == pytest.approx(
        math.sqrt(math.pi) / 2.0)
