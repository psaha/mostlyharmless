"""Spherical harmonics, in the conventions of section 01 of the manual.

    Y_l^m(theta, phi) = N_lm P_l^m(cos theta) exp(i m phi)

    N_lm = sqrt((2l+1)/(4 pi) * (l-m)! / (l+m)!)      <- the SIGNED m

with the Condon-Shortley phase carried inside P_l^m, theta the polar angle and
phi the azimuthal angle.

Two details are easy to get wrong, so they are done explicitly here rather than
left to a library's convention:

1. The factorials in N_lm take the signed m.  With ``abs(m)`` the negative
   orders come out unnormalised, by a factor [(l-|m|)!/(l+|m|)!]^2.
2. Negative orders go through the standard relation
   P_l^{-m} = (-1)^m (l-m)!/(l+m)! P_l^m for m >= 0, applied here rather than
   relying on how scipy treats a negative order.

The real basis used by the public interface is cos(m phi) for m >= 0 and
sin(-m phi) for m < 0, with an extra sqrt(2) for m != 0 so that it is
orthonormal on the sphere.

Every convention above is asserted by the certificates in
``mostlyharmless.symbolic.derivations``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import lpmv

__all__ = [
    "m_values",
    "factorial_ratio",
    "normalisation",
    "associated_legendre",
    "sph_harm",
    "real_sph_harm",
    "sph_harm_addition",
]


def m_values(l):
    """The m values for degree l, ascending.  Index i corresponds to m = i - l."""
    return np.arange(-l, l + 1)


def factorial_ratio(l, m):
    """(l-m)! / (l+m)! for -l <= m <= l, with the signed m."""
    if abs(m) > l:
        raise ValueError(f"need |m| <= l, got m={m}, l={l}")
    return math.factorial(l - m) / math.factorial(l + m)


def normalisation(l, m):
    """N_lm of equation (ylm) in section 01.  Note the signed m."""
    return math.sqrt((2 * l + 1) / (4 * math.pi) * factorial_ratio(l, m))


def associated_legendre(l, m, x):
    """P_l^m(x) with the Condon-Shortley phase, for -l <= m <= l.

    The negative-order branch uses the standard relation, written for m >= 0,

        P_l^{-m} = (-1)^m (l-m)!/(l+m)! P_l^m,

    applied explicitly rather than trusting scipy's convention for a negative
    order.  This is the one place in the package where a magnitude is needed at
    all, because lpmv takes only a non-negative order.
    """
    if abs(m) > l:
        raise ValueError(f"need -l <= m <= l, got m={m}, l={l}")
    order = abs(int(m))
    p = lpmv(order, l, x)
    if m < 0:
        p = (-1) ** order * factorial_ratio(l, order) * p
    return p


def sph_harm(l, m, theta, phi):
    """Y_l^m(theta, phi), complex, with theta the polar angle."""
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    return (normalisation(l, m)
            * associated_legendre(l, m, np.cos(theta))
            * np.exp(1j * m * phi))


def real_sph_harm(l, m, theta, phi):
    """The real basis, orthonormal on the sphere.

        m > 0 :  sqrt(2) Re Y_l^m
        m = 0 :  Y_l^0
        m < 0 :  sqrt(2) Im Y_l^{-m}

    The order -m is written out rather than as a magnitude, since m is negative
    in that branch.  No (-1)^m is inserted; several other real-basis conventions
    include one, and the difference matters when comparing to them.
    """
    if m > 0:
        return math.sqrt(2.0) * np.real(sph_harm(l, m, theta, phi))
    if m == 0:
        return np.real(sph_harm(l, 0, theta, phi))
    return math.sqrt(2.0) * np.imag(sph_harm(l, -m, theta, phi))


def sph_harm_addition(l, theta1, phi1, theta2, phi2):
    """sum_m Y_l^m(theta1,phi1) conj(Y_l^m(theta2,phi2)).

    By the addition theorem this equals (2l+1)/(4 pi) P_l(cos gamma), with

        cos gamma = cos(theta1) cos(theta2)
                  + sin(theta1) sin(theta2) cos(phi1 - phi2).
    """
    total = 0.0
    for m in range(-l, l + 1):
        total = total + sph_harm(l, m, theta1, phi1) * np.conjugate(
            sph_harm(l, m, theta2, phi2))
    return total
