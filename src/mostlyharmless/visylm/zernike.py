"""Zernike radial polynomials, and the even-family coefficients.

Two jobs live here.  :func:`zernike_radial` evaluates R_n^m(rho) for given numbers,
through the Jacobi polynomial.  :func:`mode_polynomial` and
:func:`zernike_coefficients` do the exact rational algebra the even family needs:
the first writes P_l^m(sqrt(1-rho^2)) as rho^m A(rho^2), the second expands that in
Zernike radial polynomials.

The convention is the standard one, R_n^m(1) = 1, defined for n - m even:

    R_n^m(rho) = sum_{s=0}^{(n-m)/2} (-1)^s (n-s)!
                 / (s! ((n-m)/2 + m - s)! ((n-m)/2 - s)!) rho^{n-2s} .

Two identities are used later and certified rather than assumed: the Jacobi form

    R_n^m(rho) = (-1)^{(n-m)/2} rho^m P_{(n-m)/2}^{(m,0)}(1 - 2 rho^2) ,

and the orthogonality

    int_0^1 R_n^m(rho) R_n'^m(rho) rho drho = delta_{nn'} / (2n + 2) .

The exact side uses SymPy, which is already a dependency, and every result is
memoised: a fit asks for the same coefficients many times, and they never change.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.special import eval_jacobi

__all__ = [
    "zernike_radial",
    "zernike_polynomial",
    "mode_polynomial",
    "zernike_coefficients",
]

_RHO = sp.Symbol("rho", nonnegative=True)
_X = sp.Symbol("x")


def _check_degree(n, m):
    if m < 0:
        raise ValueError(f"orders are non-negative in this package, got m={m}")
    if n < m:
        raise ValueError(f"need m <= n, got n={n}, m={m}")
    if (n - m) % 2 == 1:
        raise ValueError(f"need n - m even, got n={n}, m={m}")


def _check_even_family(l, m):
    if m < 0:
        raise ValueError(f"orders are non-negative in this package, got m={m}")
    if m > l:
        raise ValueError(f"need m <= l, got l={l}, m={m}")
    if (l + m) % 2 == 1:
        raise ValueError(f"the even family needs l + m even, got l={l}, m={m}")


def zernike_radial(n, m, rho):
    """R_n^m(rho), numeric, through the Jacobi polynomial."""
    _check_degree(n, m)
    rho = np.asarray(rho, dtype=float)
    half = (n - m) // 2
    return (-1) ** half * rho ** m * eval_jacobi(half, m, 0, 1.0 - 2.0 * rho ** 2)


@lru_cache(maxsize=None)
def zernike_polynomial(n, m):
    """R_n^m as an exact expression in rho."""
    _check_degree(n, m)
    half = (n - m) // 2
    total = sp.Integer(0)
    for s in range(half + 1):
        coefficient = (sp.Integer((-1) ** s) * sp.factorial(n - s)
                       / (sp.factorial(s) * sp.factorial(half + m - s)
                          * sp.factorial(half - s)))
        total = total + coefficient * _RHO ** (n - 2 * s)
    return sp.expand(total)


@lru_cache(maxsize=None)
def mode_polynomial(l, m):
    """P_l^m(sqrt(1-rho^2)) as rho^m A(rho^2), exactly, for l + m even.

    Uses P_l^m(x) = (-1)^m (1-x^2)^{m/2} d^m P_l / dx^m with the Condon-Shortley
    phase of section 1.  For l + m even the m-th derivative is an even function of
    x, so after x^2 -> 1 - rho^2 it is a polynomial in rho^2.
    """
    l, m = int(l), int(m)
    _check_even_family(l, m)
    derivative = sp.Poly(sp.diff(sp.legendre(l, _X), _X, m), _X)
    inside = sp.Integer(0)
    for (power,), coefficient in derivative.terms():
        if power % 2 == 1:
            raise ValueError("the derivative should be even when l + m is even")
        inside = inside + coefficient * (1 - _RHO ** 2) ** (power // 2)
    return sp.expand(sp.Integer((-1) ** m) * _RHO ** m * inside)


@lru_cache(maxsize=None)
def zernike_coefficients(l, m):
    """The even-family coefficients c_k, exact, for l + m even.

    With P_l^m(sqrt(1-rho^2)) = rho^m A(rho^2) = sum_k c_k R_{m+2k}^m(rho), the
    orthogonality gives

        c_k = (2(m + 2k) + 2) int_0^1 rho^m A(rho^2) R_{m+2k}^m(rho) rho drho .

    Returned as a tuple, one exact rational per k = 0, ..., (l-m)/2, and memoised.
    """
    l, m = int(l), int(m)
    _check_even_family(l, m)
    polynomial = mode_polynomial(l, m)
    values = []
    for k in range((l - m) // 2 + 1):
        n = m + 2 * k
        integral = sp.integrate(polynomial * zernike_polynomial(n, m) * _RHO,
                                (_RHO, 0, 1))
        values.append(sp.Rational(sp.simplify(integral)) * (2 * n + 2))
    return tuple(values)
