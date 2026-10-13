"""Machine-checked lemmas backing the manual.

Every ``check_*`` function returns a residual: zero means the claim holds.  The
tag, the rung and the scope are recorded in ROADMAP.md section 4 and quoted in
the manual next to the claim they support.

Two independent constructions are used throughout, which is what decision D10
asks for.  The exact checks build Y_l^m from SymPy's ``assoc_legendre``; the
implementation in ``mostlyharmless.visylm`` goes through SciPy's ``lpmv``.
Nothing is shared but the conventions, and those are written out explicitly in
both places.

Rungs, restated because they are easy to blur:

    exact          residual zero with the parameters free
    range-limited  residual exactly zero for every case in a stated finite set
    numerical      agreement to a stated tolerance
"""

from __future__ import annotations

import math

import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.special import eval_legendre

from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.spherical_harmonics import (
    real_sph_harm,
    sph_harm,
    sph_harm_addition,
)
from mostlyharmless.visylm.wigner import (
    rotate_coefficients,
    rotation_matrix,
    wigner_d,
    wigner_D,
)

__all__ = [
    "Y",
    "Ybar",
    "check_ylm_matches_sympy_ynm",
    "check_legendre_orthogonality_exact",
    "check_orthonormality_numeric",
    "check_reality_exact",
    "check_addition_theorem_numeric",
    "check_addition_theorem_real_numeric",
    "check_real_basis_orthonormality_numeric",
    "check_real_basis_matches_complex",
    "check_wigner_d_against_exponential",
    "check_wigner_d_properties",
    "check_rotation_against_geometry",
    "check_rotation_composition",
    "check_rotation_nonnegative_form",
]

THETA, PHI = sp.symbols("theta phi", real=True)

_ANGLE_PAIRS = [
    (sp.Rational(1, 3), sp.Rational(2, 5)),
    (sp.Rational(4, 5), -sp.Rational(7, 3)),
    (sp.Rational(9, 10), sp.Rational(11, 6)),
]


# --------------------------------------------------------------------------
# Y_l^m in SymPy, for the exact checks
# --------------------------------------------------------------------------

def _coefficient(l, m, theta=THETA):
    """The theta-dependent coefficient of exp(i m phi) in Y_l^m.

    Built from SymPy's assoc_legendre, i.e. independently of the lpmv route in
    the implementation, with the negative-order relation applied explicitly.
    The order -m is written out where it is needed, and no magnitude appears.
    """
    n = sp.sqrt((2 * l + 1) / (4 * sp.pi)
                * sp.factorial(l - m) / sp.factorial(l + m))
    if m >= 0:
        return n * sp.assoc_legendre(l, m, sp.cos(theta))
    p = (sp.factorial(l + m) / sp.factorial(l - m)
         * sp.assoc_legendre(l, -m, sp.cos(theta)))
    return n * (-1) ** (-m) * p


def Y(l, m, theta=THETA, phi=PHI):
    """Y_l^m in the conventions of the manual, as a SymPy expression."""
    return _coefficient(l, m, theta) * sp.exp(sp.I * m * phi)


def Ybar(l, m, theta=THETA, phi=PHI):
    """The complex conjugate of Y_l^m, written out explicitly.

    For real theta the associated Legendre factor is real, so conjugating only
    touches the exponential.  Doing it by construction keeps SymPy from having
    to prove that assoc_legendre is real, which it will not do.
    """
    return _coefficient(l, m, theta) * sp.exp(-sp.I * m * phi)


# --------------------------------------------------------------------------
# V1 -- orthonormality
# --------------------------------------------------------------------------

def _legendre_poly_part(l, order, x):
    """Q with P_l^order(x) = (1-x^2)^(order/2) Q(x), for order >= 0.

    The Rodrigues form, so that SymPy only ever integrates polynomials.
    """
    q = sp.legendre(l, x)
    for _ in range(order):
        q = sp.diff(q, x)
    return sp.expand((-1) ** order * q)


def check_legendre_orthogonality_exact(l_max=5):
    """V1, exact rung -- the theta factor.

        int_{-1}^{1} P_l^m(x) P_l'^m(x) dx
            = delta_ll' 2/(2l+1) (l+m)!/(l-m)!

    with the SIGNED m.  Both sides are polynomials in x, so SymPy is exact and
    fast.  Returns the sum of the absolute residuals over every case: zero means
    all of them hold.
    """
    x = sp.Symbol("x", real=True)
    total = sp.Integer(0)
    for order in range(0, l_max + 1):
        qs = {l: _legendre_poly_part(l, order, x)
              for l in range(order, l_max + 1)}
        for sign in (+1, -1):
            m = sign * order
            for l in range(order, l_max + 1):
                for lp in range(order, l_max + 1):
                    if sign > 0:
                        weight = sp.Integer(1)
                    else:
                        weight = (sp.factorial(l - order) / sp.factorial(l + order)
                                  * sp.factorial(lp - order)
                                  / sp.factorial(lp + order))
                    integrand = sp.expand((1 - x ** 2) ** order
                                          * qs[l] * qs[lp])
                    value = sp.integrate(weight * integrand, (x, -1, 1))
                    if l == lp:
                        expected = (sp.Rational(2, 2 * l + 1)
                                    * sp.factorial(l + m) / sp.factorial(l - m))
                    else:
                        expected = sp.Integer(0)
                    total += sp.Abs(sp.simplify(value - expected))
    return sp.simplify(total)


def check_ylm_matches_sympy_ynm(l_max=3, digits=30):
    """V1, numerical rung -- our conventions against SymPy's own Ynm.

    Compared at 30 digits rather than symbolically: the two are algebraically
    equal but not identically written (nested square roots), so a symbolic zero
    would test the simplifier more than the claim.
    """
    worst = 0.0
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            for th, ph in _ANGLE_PAIRS:
                ours = sp.N(Y(l, m, th, ph), digits)
                theirs = sp.N(sp.Ynm(l, m, th, ph), digits)
                worst = max(worst, abs(complex(ours - theirs)))
    return worst


def check_orthonormality_numeric(l_max=4, n_theta=None, n_phi=None):
    """V1, numerical rung -- orthonormality of the implementation end to end.

    Gauss-Legendre in cos(theta) and the uniform Fourier rule in phi.  Both are
    exact for this integrand, so what remains is floating point.
    """
    if n_theta is None:
        n_theta = 2 * l_max + 2
    if n_phi is None:
        n_phi = 4 * l_max + 4
    nodes, weights = np.polynomial.legendre.leggauss(n_theta)
    theta = np.arccos(nodes)
    phi = 2 * np.pi * np.arange(n_phi) / n_phi
    worst = 0.0
    for l in range(l_max + 1):
        for lp in range(l_max + 1):
            for m in range(-l, l + 1):
                for mp in range(-lp, lp + 1):
                    a = sph_harm(l, m, theta[:, None], phi[None, :])
                    b = sph_harm(lp, mp, theta[:, None], phi[None, :])
                    value = np.sum(weights[:, None] * np.conjugate(b) * a)
                    value = float(np.real(value)) * (2 * np.pi / n_phi)
                    expected = 1.0 if (l == lp and m == mp) else 0.0
                    worst = max(worst, abs(value - expected))
    return worst


# --------------------------------------------------------------------------
# V2 -- reality, and the negative orders
# --------------------------------------------------------------------------

def check_reality_exact(l_max=4):
    """V2, exact rung -- Y_l^{-m} = (-1)^m conj(Y_l^m).

    Reduced to the theta coefficients, so what is being tested is the
    cancellation between the factorial ratio of the normalisation and the one
    introduced by the negative-order relation.  If either is written with |m|,
    this is where it shows.
    """
    total = sp.Integer(0)
    for l in range(l_max + 1):
        for m in range(1, l + 1):
            diff = _coefficient(l, -m) - (-1) ** m * _coefficient(l, m)
            total += sp.Abs(sp.simplify(sp.powsimp(sp.expand(diff), force=True)))
    return sp.simplify(total)


# --------------------------------------------------------------------------
# V3 -- the addition theorem
# --------------------------------------------------------------------------

def check_addition_theorem_numeric(l_max=6, digits=40, seed=20261012):
    """V3, numerical rung -- the addition theorem.

        sum_m Y_l^m(1) conj(Y_l^m(2)) = (2l+1)/(4 pi) P_l(cos gamma)

    The theorem itself is classical and cited as such.  What this checks is our
    implementation of it, in two ways.  First, SymPy's evaluation of our own
    Y_l^m at exact rational angles, to `digits` digits.  Second, the float
    implementation against SymPy's Legendre, so that the reference is not the
    limiting error.

    The exact rung is deliberately not attempted: a symbolic route has to
    reconcile sqrt(1-cos^2 theta) with sin(theta), and the sign of sin(theta)
    turns that into a statement about assumptions rather than about the theorem.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    for l in range(l_max + 1):
        # (a) exact rational angles, evaluated to `digits` digits.
        #     theta is kept inside (0, pi), where sin(theta) > 0, so that the
        #     principal square root in assoc_legendre is the right branch.
        for _ in range(4):
            a1 = sp.Rational(int(rng.integers(1, 30)), 10)
            a2 = sp.Rational(int(rng.integers(1, 30)), 10)
            b1 = sp.Rational(int(rng.integers(-30, 31)), 10)
            b2 = sp.Rational(int(rng.integers(-30, 31)), 10)
            lhs = sum(Y(l, m, a1, b1) * Ybar(l, m, a2, b2)
                      for m in range(-l, l + 1))
            cos_gamma = (sp.cos(a1) * sp.cos(a2)
                         + sp.sin(a1) * sp.sin(a2) * sp.cos(b1 - b2))
            rhs = (2 * l + 1) / (4 * sp.pi) * sp.legendre(l, cos_gamma)
            worst = max(worst, abs(complex(sp.N(lhs - rhs, digits))))

        # (b) the implementation itself, in floating point.
        for _ in range(4):
            t1 = rng.uniform(0.1, np.pi - 0.1)
            t2 = rng.uniform(0.1, np.pi - 0.1)
            p1 = rng.uniform(-np.pi, np.pi)
            p2 = rng.uniform(-np.pi, np.pi)
            cos_gamma = (math.cos(t1) * math.cos(t2)
                         + math.sin(t1) * math.sin(t2) * math.cos(p1 - p2))
            rhs = (2 * l + 1) / (4 * math.pi) * float(
                sp.legendre(l, sp.Float(cos_gamma, 30)))
            lhs = float(np.real(sph_harm_addition(l, t1, p1, t2, p2)))
            worst = max(worst, abs(lhs - rhs))
    return worst


# --------------------------------------------------------------------------
# The real basis
# --------------------------------------------------------------------------

def check_addition_theorem_real_numeric(l_max=6, seed=20261012):
    """V3, numerical rung -- the addition theorem in the real basis.

        R_l^0(1) R_l^0(2)
          + sum_{m>0} [ R_l^{m,c}(1) R_l^{m,c}(2) + R_l^{m,s}(1) R_l^{m,s}(2) ]
        = (2l+1)/(4 pi) P_l(cos gamma)

    This is the form the manual states, and every index in it is non-negative.
    In the implementation's calling convention the cosine branch of order m is
    reached as ``real_sph_harm(l, m, ...)`` and the sine branch as
    ``real_sph_harm(l, -m, ...)``; that labelling is item 6 of the roadmap's open
    list, and the correspondence is what matters here.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    for l in range(l_max + 1):
        for _ in range(6):
            t1 = rng.uniform(0.1, np.pi - 0.1)
            t2 = rng.uniform(0.1, np.pi - 0.1)
            p1 = rng.uniform(-np.pi, np.pi)
            p2 = rng.uniform(-np.pi, np.pi)
            cos_gamma = (math.cos(t1) * math.cos(t2)
                         + math.sin(t1) * math.sin(t2) * math.cos(p1 - p2))
            total = real_sph_harm(l, 0, t1, p1) * real_sph_harm(l, 0, t2, p2)
            for m in range(1, l + 1):
                total += (real_sph_harm(l, m, t1, p1) * real_sph_harm(l, m, t2, p2)
                          + real_sph_harm(l, -m, t1, p1)
                          * real_sph_harm(l, -m, t2, p2))
            expected = (2 * l + 1) / (4 * math.pi) * float(
                sp.legendre(l, sp.Float(cos_gamma, 30)))
            worst = max(worst, abs(total - expected))
    return worst


def check_real_basis_orthonormality_numeric(l_max=4, n_theta=None, n_phi=None):
    """The interface basis is orthonormal on the sphere.

    Real functions, so no conjugation is needed.
    """
    if n_theta is None:
        n_theta = 2 * l_max + 2
    if n_phi is None:
        n_phi = 4 * l_max + 4
    nodes, weights = np.polynomial.legendre.leggauss(n_theta)
    theta = np.arccos(nodes)
    phi = 2 * np.pi * np.arange(n_phi) / n_phi
    worst = 0.0
    for l in range(l_max + 1):
        for lp in range(l_max + 1):
            for m in range(-l, l + 1):
                for mp in range(-lp, lp + 1):
                    f = (real_sph_harm(l, m, theta[:, None], phi[None, :])
                         * real_sph_harm(lp, mp, theta[:, None], phi[None, :]))
                    value = float(np.sum(weights[:, None] * f)) * (2 * np.pi / n_phi)
                    expected = 1.0 if (l == lp and m == mp) else 0.0
                    worst = max(worst, abs(value - expected))
    return worst


def check_real_basis_matches_complex(l_max=4):
    """The real basis agrees with the convention stated in spherical_harmonics.

    A contract check rather than a certificate: it pins the code to the
    documented definition, sqrt(2) Re or Im of the complex harmonic.
    """
    worst = 0.0
    thetas = np.linspace(0.05, np.pi - 0.05, 7)
    phis = np.linspace(-3.0, 3.0, 7)
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            for theta in thetas:
                for phi in phis:
                    if m > 0:
                        expected = math.sqrt(2.0) * np.real(
                            sph_harm(l, m, theta, phi))
                    elif m == 0:
                        expected = np.real(sph_harm(l, 0, theta, phi))
                    else:
                        expected = math.sqrt(2.0) * np.imag(
                            sph_harm(l, -m, theta, phi))
                    worst = max(worst, abs(real_sph_harm(l, m, theta, phi) - expected))
    return worst


# --------------------------------------------------------------------------
# V4-V7 -- the rotation
# --------------------------------------------------------------------------

def _jy_matrix(l):
    """(J_y)_{m'm} in the |l m> basis with hbar = 1, index [m' + l, m + l]."""
    Jp = np.zeros((2 * l + 1, 2 * l + 1))
    Jm = np.zeros((2 * l + 1, 2 * l + 1))
    for m in range(-l, l + 1):
        if m + 1 <= l:
            Jp[m + 1 + l, m + l] = math.sqrt((l - m) * (l + m + 1))
        if m - 1 >= -l:
            Jm[m - 1 + l, m + l] = math.sqrt((l + m) * (l - m + 1))
    return (Jp - Jm) / 2j


def check_wigner_d_against_exponential(l_max=6, n_beta=6):
    """V4, numerical rung -- the explicit d sum against the matrix exponential.

        d^l_{m'm}(beta) = <l m'| exp(-i beta J_y) |l m>

    Two independent constructions of the same object: a factorial sum in
    half-angles, and the exponential of the angular momentum matrix.
    """
    worst = 0.0
    for l in range(l_max + 1):
        Jy = _jy_matrix(l)
        for beta in np.linspace(0.1, math.pi - 0.1, n_beta):
            reference = expm(-1j * beta * Jy)
            worst = max(worst, float(np.max(np.abs(wigner_d(l, beta) - reference))))
    return worst


def check_wigner_d_properties(l_max=6):
    """V5, numerical rung -- the elementary properties of d.

    Reality and orthogonality, d(0) = I, d^l_{00} = P_l(cos beta), and the
    symmetry d^l_{m'm} = (-1)^(m-m') d^l_{mm'}.
    """
    worst = 0.0
    for l in range(l_max + 1):
        size = 2 * l + 1
        for beta in (0.0, 0.3, 1.1, 2.5, math.pi):
            d = wigner_d(l, beta)
            worst = max(worst, float(np.max(np.abs(d.T @ d - np.eye(size)))))
            if beta == 0.0:
                worst = max(worst, float(np.max(np.abs(d - np.eye(size)))))
            worst = max(worst, abs(d[l, l] - eval_legendre(l, math.cos(beta))))
            for i in range(size):
                for j in range(size):
                    worst = max(worst, abs(d[i, j] - (-1) ** (j - i) * d[j, i]))
    return worst


def _map_value(coeffs, theta, phi):
    """The real map from the stored non-negative half, by the manual's formula."""
    total = 0.0
    for l, block in enumerate(coeffs.blocks):
        total += block[0] * sph_harm(l, 0, theta, phi)
        for m in range(1, l + 1):
            total += 2.0 * np.real(block[m] * sph_harm(l, m, theta, phi))
    return float(np.real(total))


def _angles_from(v):
    """Polar and azimuthal angles of a 3-vector."""
    return (math.acos(max(-1.0, min(1.0, v[2]))), math.atan2(v[1], v[0]))


def _random_direction(rng):
    v = rng.normal(size=3)
    return v / np.linalg.norm(v)


def check_rotation_against_geometry(l_max=4, seed=20261012, n_rot=3, n_dir=3):
    """V6, numerical rung -- the coefficient transform against direct geometry.

    Rotating the coefficients must equal evaluating the original map at the
    inversely rotated direction,

        I_rotated(n) = I_body(R^T n) .

    The reference shares nothing with the algebra.
    """
    rng = np.random.default_rng(seed)
    coeffs = Coeffs.from_vector(rng.normal(size=(l_max + 1) ** 2))
    worst = 0.0
    for _ in range(n_rot):
        alpha = rng.uniform(0.0, 2 * math.pi)
        beta = rng.uniform(0.0, math.pi)
        pa = rng.uniform(0.0, 2 * math.pi)
        rotated = rotate_coefficients(coeffs, alpha, beta, pa)
        R = rotation_matrix(alpha, beta, pa)
        for _ in range(n_dir):
            v = _random_direction(rng)
            theta_r, phi_r = _angles_from(v)
            theta_b, phi_b = _angles_from(R.T @ v)
            worst = max(worst, abs(_map_value(rotated, theta_r, phi_r)
                                   - _map_value(coeffs, theta_b, phi_b)))
    return worst


def check_rotation_composition(l_max=3, seed=20261012, n_dir=4):
    """V7, numerical rung -- two rotations against the composite direction.

    Applying R1 and then R2 must equal the single rotation R2 R1, which the
    geometric reference checks without needing any Euler angles.
    """
    rng = np.random.default_rng(seed)
    coeffs = Coeffs.from_vector(rng.normal(size=(l_max + 1) ** 2))
    first = (rng.uniform(0.0, 2 * math.pi), rng.uniform(0.0, math.pi),
             rng.uniform(0.0, 2 * math.pi))
    second = (rng.uniform(0.0, 2 * math.pi), rng.uniform(0.0, math.pi),
              rng.uniform(0.0, 2 * math.pi))
    twice = rotate_coefficients(rotate_coefficients(coeffs, *first), *second)
    R = rotation_matrix(*second) @ rotation_matrix(*first)
    worst = 0.0
    for _ in range(n_dir):
        v = _random_direction(rng)
        theta_r, phi_r = _angles_from(v)
        theta_b, phi_b = _angles_from(R.T @ v)
        worst = max(worst, abs(_map_value(twice, theta_r, phi_r)
                               - _map_value(coeffs, theta_b, phi_b)))
    return worst


def check_rotation_nonnegative_form(l_max=4, seed=20261012):
    """V6, part two -- the manual's transform, with no negative-order coefficient.

        b_{lm'} = D^l_{m'0} a_{l0}
                  + sum_{m=1..l} [ D^l_{m'm} a_{lm}
                                   + (-1)^m D^l_{m',-m} conj(a_{lm}) ] ,  m' >= 0

    checked against the full-range transform that ``rotate_coefficients`` uses.
    The m = 0 term appears once and not twice, since a_{l0} is real.
    """
    rng = np.random.default_rng(seed)
    coeffs = Coeffs.from_vector(rng.normal(size=(l_max + 1) ** 2))
    alpha = rng.uniform(0.0, 2 * math.pi)
    beta = rng.uniform(0.0, math.pi)
    pa = rng.uniform(0.0, 2 * math.pi)
    rotated = rotate_coefficients(coeffs, alpha, beta, pa)
    worst = 0.0
    for l in range(l_max + 1):
        # D^l_{m'm}(R) for the viewing angles is wigner_D with the zyz angles
        # (pi/2 - PA, beta, alpha)
        D = wigner_D(l, math.pi / 2.0 - pa, beta, alpha)
        for mp in range(l + 1):
            total = D[mp + l, l] * coeffs[l, 0].real
            for m in range(1, l + 1):
                total += (D[mp + l, m + l] * coeffs[l, m]
                          + (-1) ** m * D[mp + l, l - m]
                          * np.conjugate(coeffs[l, m]))
            worst = max(worst, abs(total - rotated[l, mp]))
    return worst
