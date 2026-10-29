"""The visibility transform: from a map's coefficients to V(u, v).

Conventions are section 01 of the manual.  Two things are worth stating here
because both are easy to get wrong and both are certified:

* the measure is planar over the projected disk, so the Jacobi-Anger step turns
  the azimuthal integral into a Bessel function and leaves

      V_lm = 2 pi (-i)^m e^{i m psi} N_lm int_0^1 P_{lm}(sqrt(1-rho^2))
                                          J_m(2 pi q rho) rho drho ,

  with q and psi the polar form of (u, v);

* modes are labelled by non-negative orders only.  Negative orders are never
  computed: they are related to the non-negative ones by

      V_{l,-m} = conj(V_{lm}) ,

  which is the analogue of the harmonic extension Y_{l,-m} = (-1)^m conj(Y_{lm}),
  with the sign absent because the (-1)^m is cancelled when the baseline is
  reversed.  That relation is what lets a stored coefficient pair up with its
  partner, and it is the only place a negative order is needed at all.

The radial integral is evaluated by quadrature rather than in closed form.  That
is deliberate: it is the reference implementation, and the closed forms of the odd
and even families are checked against it.  It must never dispatch to a closed form
itself, or the check would become a comparison of the formula with itself.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import jv, spherical_jn

from mostlyharmless.visylm.spherical_harmonics import (
    associated_legendre,
    normalisation,
)
from mostlyharmless.visylm.zernike import zernike_coefficients

__all__ = [
    "radial_integral",
    "visibility_mode",
    "visibility",
    "visibility_uniform_disk",
]


def _uv_polar(u, v):
    """The radial baseline and its azimuth."""
    return np.hypot(u, v), np.arctan2(v, u)


def _double_factorial(n):
    """n!!, with 0!! = (-1)!! = 1.  Small non-negative arguments only."""
    if n <= 0:
        return 1
    result = 1
    for k in range(n, 0, -2):
        result *= k
    return result


def _spherical_bessel_over_argument(l, a):
    """j_l(a)/a, including its value at a = 0.

    The singularity is removable.  Since j_l(a) ~ a^l/(2l+1)!!, the limit is 1/3
    for l = 1 and 0 for l >= 2.  For every a != 0 plain division is accurate ---
    scipy evaluates j_l directly, so there is no cancellation to lose --- and only
    a = 0 needs a case of its own.  A small-but-nonzero a is deliberately not
    special-cased: a series would introduce a truncation error where division
    introduces none.
    """
    a = np.asarray(a, dtype=float)
    limit = 1.0 / 3.0 if l == 1 else 0.0
    return np.where(a == 0.0, limit,
                    spherical_jn(l, a) / np.where(a == 0.0, 1.0, a))


def _radial_integral_odd(l, m, q):
    """The odd-family closed form, for l + m odd.

        R_lm(q) = (l+m)!!/(l-m-1)!! j_l(2 pi q)/(2 pi q)

    Kept private, and reached through :func:`visibility_mode`.  The quadrature in
    :func:`radial_integral` is the reference this is checked against, so that
    function must never dispatch here.
    """
    a = 2.0 * math.pi * np.asarray(q, dtype=float)
    coefficient = _double_factorial(l + m) / _double_factorial(l - m - 1)
    return coefficient * _spherical_bessel_over_argument(l, a)


def _bessel_over_argument(order, a):
    """J_order(a)/a, including its value at a = 0.

    Since J_nu(a) ~ a^nu, the limit is 1/2 for order 1 and 0 for order above 1.
    As with the spherical case, division is accurate for every a != 0 and only
    a = 0 needs a case of its own.
    """
    a = np.asarray(a, dtype=float)
    limit = 0.5 if order == 1 else 0.0
    return np.where(a == 0.0, limit,
                    jv(order, a) / np.where(a == 0.0, 1.0, a))


def _radial_integral_even(l, m, q):
    """The even-family closed form, for l + m even.

        R_lm(q) = sum_k c_k (-1)^k J_{m+2k+1}(2 pi q)/(2 pi q)

    with the exact coefficients c_k of :func:`zernike_coefficients`.  At q = 0 only
    the term with m + 2k = 0 survives, so the limit is c_0/2 for m = 0 and zero
    otherwise.  For l = m = 0 this is J_1(2 pi q)/(2 pi q), the uniform disk.
    """
    a = 2.0 * math.pi * np.asarray(q, dtype=float)
    total = 0.0
    for k, coefficient in enumerate(zernike_coefficients(l, m)):
        total = total + (float(coefficient) * (-1) ** k
                         * _bessel_over_argument(m + 2 * k + 1, a))
    return total


def _n_nodes(q):
    """Gauss-Legendre nodes enough for the oscillation of the Bessel factor.

    The integrand oscillates about 2 pi q times across the disk, so the count
    grows with q.  Both the count and the rule are part of what makes this a
    reliable reference, so they are stated rather than tuned silently.
    """
    q_max = float(np.max(np.abs(np.atleast_1d(q)))) if np.size(q) else 0.0
    return max(64, 16 * (int(2.0 * math.pi * q_max) + 2))


def radial_integral(l, m, q, n_nodes=None):
    """int_0^1 P_{lm}(sqrt(1-rho^2)) J_m(2 pi q rho) rho drho.

    The intermediate of the master reduction, without the harmonic
    normalisation.

    The quadrature is done in theta, with rho = sin(theta), so that the integrand
    is analytic on the interval:

        int_0^{pi/2} P_{lm}(cos theta) J_m(2 pi q sin theta)
                     sin(theta) cos(theta) dtheta .

    Integrating in rho instead would carry the factor sqrt(1-rho^2) to the end
    point, where its derivative diverges.  That is an integrable endpoint
    singularity, so the integral still converges, but only algebraically --- 64
    Gauss-Legendre nodes leave an error near 1e-6 for the odd-parity modes.  This
    was found by the contract check against scipy's adaptive quadrature, which is
    why that check exists.

    There is no division by q, so q = 0 needs no special case.  This function is
    the reference the closed forms of stages 3 and 4 are checked against.
    """
    m = int(m)
    if m < 0:
        raise ValueError(f"orders are non-negative in this package, got m={m}")
    if m > l:
        raise ValueError(f"need m <= l, got l={l}, m={m}")

    q = np.asarray(q, dtype=float)
    q_flat = np.atleast_1d(q).ravel()
    if n_nodes is None:
        n_nodes = _n_nodes(q_flat)

    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)
    theta = 0.25 * math.pi * (nodes + 1.0)          # [-1, 1] -> [0, pi/2]
    quarter = 0.25 * math.pi * weights

    sin_theta = np.sin(theta)
    legendre = associated_legendre(l, m, np.cos(theta))
    argument = 2.0 * math.pi * q_flat[:, None] * sin_theta[None, :]
    integrand = (legendre[None, :] * jv(m, argument)
                 * sin_theta[None, :] * np.cos(theta)[None, :] * quarter[None, :])
    values = np.sum(integrand, axis=1)
    return values.reshape(q.shape)


def visibility_mode(l, m, u, v):
    """The complex visibility of the harmonic Y_{lm}, for m >= 0.

        V_{lm} = 2 pi (-i)^m N_lm exp(i m psi) R_lm(q)

    For ``l + m`` odd the radial integral is the closed form of section 5, a single
    spherical Bessel function; for ``l + m`` even it is the closed form of section
    6, a finite sum of ordinary Bessel functions through the Zernike expansion.
    Both are dispatched here, and the quadrature of :func:`radial_integral` is
    never used in production --- it is the reference the two are checked against.

    The negative orders are not computed here.  They are given by
    ``V_{l,-m} = conj(V_{lm})``, which is what :func:`visibility` uses to pair a
    stored coefficient with its partner.
    """
    m = int(m)
    if m < 0:
        raise ValueError(f"orders are non-negative in this package, got m={m}")
    if m > l:
        raise ValueError(f"need m <= l, got l={l}, m={m}")
    q, psi = _uv_polar(u, v)
    if (l + m) % 2 == 1:
        radial = _radial_integral_odd(l, m, q)
    else:
        radial = _radial_integral_even(l, m, q)
    return (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
            * np.exp(1j * m * psi) * radial)


def visibility(coeffs, u, v):
    """The complex visibility of a map, by linearity over the harmonics.

    A stored coefficient a_lm with m > 0 contributes

        a_lm V_{lm} + (-1)^m conj(a_lm) conj(V_{lm}) ,

    the second term being the partner that the reality of the map determines.
    Their sum is 4 pi (-i)^m N_lm R_lm(q) Re(a_lm exp(i m psi)).  At m = 0 there
    is a single term, a_l0 V_{l0}.
    """
    total = 0.0
    for l, block in enumerate(coeffs.blocks):
        total = total + block[0] * visibility_mode(l, 0, u, v)
        for m in range(1, l + 1):
            mode = visibility_mode(l, m, u, v)
            total = total + (block[m] * mode
                             + (-1) ** m * np.conjugate(block[m])
                             * np.conjugate(mode))
    return total


def visibility_uniform_disk(u, v, theta_d):
    """The normalised visibility of a uniform disk, 2 J_1(x)/x, x = 2 pi q theta_d.

    Normalised means divided by its value at the origin, which is the total flux.
    This is the classical form the uniform map must reproduce.
    """
    q, _ = _uv_polar(u, v)
    x = 2.0 * math.pi * q * theta_d
    safe = np.where(x == 0.0, 1.0, x)
    return np.where(x == 0.0, 1.0, 2.0 * jv(1, safe) / safe)
