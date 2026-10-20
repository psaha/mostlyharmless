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
from scipy.special import eval_legendre, jv, spherical_jn

from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.spherical_harmonics import (
    associated_legendre,
    normalisation,
    real_sph_harm,
    sph_harm,
    sph_harm_addition,
)
from mostlyharmless.visylm.visibility import (
    radial_integral,
    visibility,
    visibility_mode,
    visibility_uniform_disk,
)
from mostlyharmless.visylm.wigner import (
    rotate_coefficients,
    rotation_matrix,
    wigner_d,
    wigner_D,
)
from mostlyharmless.visylm.zernike import (
    mode_polynomial,
    zernike_coefficients,
    zernike_polynomial,
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
    "check_jacobi_anger",
    "check_uniform_disk_visibility",
    "check_odd_collapse_exact",
    "check_k_identity",
    "check_odd_closed_form",
    "check_zernike_jacobi",
    "check_zernike_orthogonality",
    "check_zernike_hankel",
    "check_zernike_projection",
    "check_even_closed_form",
    "check_full_map_fft",
    "check_deformation_first_order",
    "check_level_surface",
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
        # D^l_{m'm}(R) for the viewing angles, with R = R_z(PA) R_y(beta) R_z(alpha):
        # the row angle is PA and the column angle is alpha
        D = wigner_D(l, pa, beta, alpha)
        for mp in range(l + 1):
            total = D[mp + l, l] * coeffs[l, 0].real
            for m in range(1, l + 1):
                total += (D[mp + l, m + l] * coeffs[l, m]
                          + (-1) ** m * D[mp + l, l - m]
                          * np.conjugate(coeffs[l, m]))
            worst = max(worst, abs(total - rotated[l, mp]))
    return worst


# --------------------------------------------------------------------------
# V16, V17 -- the master reduction
# --------------------------------------------------------------------------

def check_jacobi_anger(m_max=3, n_phi=512):
    """V17, numerical rung -- the azimuthal factors of the reduction.

        int_0^{2 pi} cos(m phi) exp(-i z cos(phi - psi)) dphi
            = 2 pi (-i)^m J_m(z) cos(m psi) ,

    and the same with sin.  This is the step that produces the master reduction;
    it is a statement about Bessel functions, independent of any convention of
    ours.  The integrand is periodic, so the uniform rule converges fast and the
    residual is at the roundoff level.
    """
    phi = 2.0 * math.pi * np.arange(n_phi) / n_phi
    worst = 0.0
    for m in range(m_max + 1):
        for z in (0.5, 2.3, 7.0):
            for psi in (0.0, 0.7, 2.5):
                kernel = np.exp(-1j * z * np.cos(phi - psi))
                cosine = np.sum(np.cos(m * phi) * kernel) * (2.0 * math.pi / n_phi)
                sine = np.sum(np.sin(m * phi) * kernel) * (2.0 * math.pi / n_phi)
                common = 2.0 * math.pi * (-1j) ** m * jv(m, z)
                worst = max(worst,
                            abs(cosine - common * np.cos(m * psi)),
                            abs(sine - common * np.sin(m * psi)))
    return worst


def check_uniform_disk_visibility(q_values=(0.0, 0.3, 1.0, 2.5, 5.0, 9.0), l_max=0):
    """V16, numerical rung -- the uniform map against the Airy form.

    A map with a_00 = 1 and every other coefficient zero is the uniform
    brightness case, whose normalised visibility is 2 J_1(x)/x.  Two things are
    checked: the flux at the origin, which must be sqrt(pi)/2, and the normalised
    visibilities at a range of baselines.

    This is the end-to-end anchor for the whole convention stack: the harmonic
    normalisation, the placement of 2 pi, and the planar measure all enter, and
    the classical result fixes them jointly.
    """
    vector = np.zeros((l_max + 1) ** 2)
    vector[0] = 1.0
    coeffs = Coeffs.from_vector(vector)

    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))
    worst = abs(flux - math.sqrt(math.pi) / 2.0)

    for q in q_values:
        value = np.real(visibility(coeffs, q, 0.0)) / flux
        worst = max(worst, abs(value - visibility_uniform_disk(q, 0.0, 1.0)))
    return worst


# --------------------------------------------------------------------------
# V11, V13, V18 -- the odd family
# --------------------------------------------------------------------------

def _factorial2(n):
    """n!! as an exact integer, with 0!! = (-1)!! = 1."""
    if n <= 0:
        return sp.Integer(1)
    return sp.factorial2(n)


def check_odd_collapse_exact(l_max=8):
    """V11, exact rung -- the two coefficients of the collapse coincide.

    With the K identity substituted into the doubled integral, the odd-family
    result is

        (1/(2l+1)) [ (l-m+1) K_{l+1} + (l+m) K_{l-1} ]
            -> (2(-1)^m/(2l+1)) [ C_R j_{l+1} + C_L j_{l-1} ] ,

    and the claim is that C_R and C_L, once the double factorials are tidied, are
    the same number.  That is what lets the spherical Bessel recurrence fuse the
    two terms into one.  Exact rational arithmetic for every valid (l, m), so this
    is a proof of the finite statement rather than evidence for it.
    """
    worst = sp.Integer(0)
    for l in range(1, l_max + 1):
        for m in range(0, l + 1):
            if (l + m) % 2 == 0:
                continue
            common = (2 * (-1) ** m * _factorial2(l + m)
                      / ((2 * l + 1) * _factorial2(l - m - 1)))
            right = (2 * (-1) ** m * (l - m + 1) * _factorial2(l + m)
                     / ((2 * l + 1) * _factorial2(l - m + 1)))
            left = (2 * (-1) ** m * (l + m) * _factorial2(l + m - 2)
                    / ((2 * l + 1) * _factorial2(l - m - 1)))
            worst = sp.Max(worst, sp.Abs(right - common), sp.Abs(left - common))
    return sp.simplify(worst)


def check_k_identity(l_max=6, n_theta=200):
    """V18, numerical rung -- the K identity that the odd derivation cites.

        int_{-1}^{1} P_l^m(z) J_m(a sqrt(1-z^2)) dz
            = (-1)^m 2 (l+m-1)!!/(l-m)!! j_l(a) ,     l + m even,

    in our Condon-Shortley convention.  The sign is the one thing the citation
    does not settle for us, so it is checked rather than assumed.  Gauss-Legendre
    in theta, where the integrand is analytic.
    """
    nodes, weights = np.polynomial.legendre.leggauss(n_theta)
    theta = 0.5 * math.pi * (nodes + 1.0)          # [-1, 1] -> [0, pi]
    weight = 0.5 * math.pi * weights
    sine = np.sin(theta)
    worst = 0.0
    for l in range(l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                continue
            for a in (0.4, 2.0, 5.5):
                integrand = (associated_legendre(l, m, np.cos(theta))
                             * jv(m, a * sine))
                value = float(np.sum(weight * integrand * sine))
                expected = float((-1) ** m * 2 * _factorial2(l + m - 1)
                                 / _factorial2(l - m)) * float(spherical_jn(l, a))
                worst = max(worst, abs(value - expected))
    return worst


def check_odd_closed_form(l_max=5, q_values=(0.0, 0.3, 1.0, 2.5, 5.0, 9.0)):
    """V13, numerical rung -- the odd closed form against the quadrature.

    Compared through ``visibility_mode``, so the dispatch is checked as well as
    the closed form.  The quadrature is the reference and is never itself
    dispatched to a closed form, which is what keeps this an independent
    comparison rather than a tautology.  The grid includes q = 0, where the closed
    form needs its limit.
    """
    worst = 0.0
    for l in range(1, l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 0:
                continue
            for q in q_values:
                mode = complex(visibility_mode(l, m, q, 0.0))
                reference = (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
                             * complex(radial_integral(l, m, q)))
                worst = max(worst, abs(mode - reference))
    return worst


# --------------------------------------------------------------------------
# V8, V9, V10, V12, V14 -- the even family
# --------------------------------------------------------------------------

def check_zernike_jacobi(n_max=8):
    """V8, exact rung -- the Jacobi form of the Zernike radial polynomial.

        R_n^m(rho) = (-1)^k rho^m P_k^(m,0)(1 - 2 rho^2) ,   k = (n-m)/2 ,

    compared as polynomials by expansion.  Exact over the stated range.
    """
    rho = sp.Symbol("rho", nonnegative=True)
    for n in range(n_max + 1):
        for m in range(n + 1):
            if (n - m) % 2 == 1:
                continue
            k = (n - m) // 2
            jacobi = (-1) ** k * rho ** m * sp.jacobi(k, m, 0, 1 - 2 * rho ** 2)
            difference = sp.expand(zernike_polynomial(n, m) - jacobi)
            if difference != 0:
                return difference
    return sp.Integer(0)


def check_zernike_orthogonality(n_max=8):
    """V9, exact rung -- int_0^1 R_n^m R_n'^m rho drho = delta_{nn'}/(2n+2)."""
    rho = sp.Symbol("rho", nonnegative=True)
    worst = sp.Integer(0)
    for m in range(n_max + 1):
        orders = [n for n in range(m, n_max + 1) if (n - m) % 2 == 0]
        for n in orders:
            for other in orders:
                integral = sp.integrate(zernike_polynomial(n, m)
                                        * zernike_polynomial(other, m) * rho,
                                        (rho, 0, 1))
                expected = (sp.Rational(1, 2 * n + 2) if n == other
                            else sp.Integer(0))
                worst = sp.Max(worst, sp.Abs(sp.simplify(integral - expected)))
    return sp.simplify(worst)


def check_zernike_hankel(n_max=6, order=6):
    """V10, exact rung -- the Zernike-Hankel identity, through its power series.

        int_0^1 R_n^m(rho) J_m(a rho) rho drho = (-1)^k J_{n+1}(a)/a ,
        k = (n-m)/2 .

    Both sides are expanded in a: the left from the exact moments of R_n^m and the
    Bessel series, the right from the Bessel series alone.  Every power up to
    a^(n + 2*order) is compared, so the rung is exact with a stated series scope
    rather than a proof of the identity for all a.  The identity itself is cited;
    what this fixes is its sign and normalisation in our conventions.
    """
    rho = sp.Symbol("rho", nonnegative=True)
    worst = sp.Integer(0)
    for n in range(n_max + 1):
        for m in range(n + 1):
            if (n - m) % 2 == 1:
                continue
            k = (n - m) // 2
            left = {}
            for j in range(k + order + 1):
                moment = sp.integrate(
                    zernike_polynomial(n, m) * rho ** (m + 2 * j + 1),
                    (rho, 0, 1))
                power = m + 2 * j
                term = (sp.Rational((-1) ** j, 2 ** (m + 2 * j)) * moment
                        / (sp.factorial(j) * sp.factorial(m + j)))
                left[power] = left.get(power, sp.Integer(0)) + term
            right = {}
            for j in range(order + 1):
                power = n + 2 * j
                right[power] = (sp.Rational((-1) ** (k + j), 2 ** (n + 1 + 2 * j))
                                / (sp.factorial(j) * sp.factorial(n + 1 + j)))
            for power in sorted(set(left) | set(right)):
                worst = sp.Max(worst, sp.Abs(left.get(power, sp.Integer(0))
                                             - right.get(power, sp.Integer(0))))
    return sp.simplify(worst)


def check_zernike_projection(l_max=6):
    """V12, exact rung -- the projection coefficients equal coefficient matching.

    The same c_k obtained two ways: from the orthogonality integral, and by
    matching powers of rho from the top down.  The highest power in
    P_l^m = sum_k c_k R_{m+2k}^m occurs in R_l^m alone, so the system is triangular
    once solved downwards, and the leftover polynomial must vanish identically.
    Exact rational arithmetic throughout.
    """
    rho = sp.Symbol("rho", nonnegative=True)
    worst = sp.Integer(0)
    for l in range(l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                continue
            order = (l - m) // 2
            residual = sp.expand(mode_polynomial(l, m))
            matched = [None] * (order + 1)
            for k in range(order, -1, -1):
                n = m + 2 * k
                leading = sp.binomial(n, (n - m) // 2)
                coefficient = (sp.Poly(residual, rho).coeff_monomial(rho ** n)
                               / leading)
                matched[k] = coefficient
                residual = sp.expand(residual
                                     - coefficient * zernike_polynomial(n, m))
            worst = sp.Max(worst, sp.Abs(sp.simplify(residual)))
            for projected, by_matching in zip(zernike_coefficients(l, m), matched):
                worst = sp.Max(worst,
                               sp.Abs(sp.simplify(projected - by_matching)))
    return sp.simplify(worst)


def check_even_closed_form(l_max=6, q_values=(0.0, 0.3, 1.0, 2.5, 5.0, 9.0)):
    """V14, numerical rung -- the even closed form against the quadrature.

    Through ``visibility_mode``, so the dispatch is checked as well, and the
    quadrature is the reference.  Includes q = 0, where the limit is c_0/2 for
    m = 0 and zero otherwise.  At l = m = 0 this repeats the uniform-disk anchor
    by a different route.
    """
    worst = 0.0
    for l in range(l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                continue
            for q in q_values:
                mode = complex(visibility_mode(l, m, q, 0.0))
                reference = (2.0 * math.pi * (-1j) ** m * normalisation(l, m)
                             * complex(radial_integral(l, m, q)))
                worst = max(worst, abs(mode - reference))
    return worst


# --------------------------------------------------------------------------
# V15 -- the full map against a 2-D FFT
# --------------------------------------------------------------------------

def _zonal_starspot(l_max, half_angle=0.9, depth=0.6):
    """A starspot on a uniformly bright surface, as zonal coefficients with a_00 = 1.

    The surface has brightness 1 outside the cap and 1 - depth inside it.  Built
    here rather than imported from the figure module, so that the certificate does
    not depend on the plotting code.
    """
    cos_edge = math.cos(half_angle)
    nodes, weights = np.polynomial.legendre.leggauss(200)
    mu = 0.5 * (nodes + 1.0) * (1.0 - cos_edge) + cos_edge
    weight = 0.5 * (1.0 - cos_edge) * weights
    moments = [float(np.sum(weight * eval_legendre(l, mu))) for l in range(l_max + 1)]
    zeroth = 2.0 - depth * moments[0]
    blocks = []
    for l in range(l_max + 1):
        inside = (2.0 if l == 0 else 0.0) - depth * moments[l]
        block = np.zeros(l + 1, dtype=complex)
        block[0] = (inside / zeroth) * math.sqrt(2 * l + 1)
        blocks.append(block)
    return Coeffs(tuple(blocks))


def _fft_visibility(coeffs, n_side=129):
    """The visibility by a 2-D FFT of the sampled image, in our planar measure.

    The image is sampled at x, y = j h on a grid of odd length, h = 2/n_side,
    with the origin at the centre.  An odd length matters: it makes the sample points
    symmetric about the origin, so the discrete transform needs no phase correction.
    The measure enters as the pixel area h^2, and the ifftshift/fftshift pair puts
    the zero frequency at the centre of the array.
    """
    half = (n_side - 1) // 2
    step = 2.0 / n_side
    axis = step * np.arange(-half, half + 1)
    x_east, x_north = np.meshgrid(axis, axis)
    rho_squared = x_east ** 2 + x_north ** 2
    z = np.sqrt(np.clip(1.0 - rho_squared, 0.0, None))
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    # By D24 the body x axis is North, so the harmonics' azimuth is measured from
    # North: the same frame that psi = arctan2(v, u) lives in.
    phi = np.arctan2(x_east, x_north)

    image = np.zeros_like(x_east)
    real = coeffs.to_real()
    for l, row in enumerate(real):
        image = image + row[l] * real_sph_harm(l, 0, theta, phi)
        for m in range(1, l + 1):
            image = image + (row[l + m] * real_sph_harm(l, m, theta, phi)
                             + row[l - m] * real_sph_harm(l, -m, theta, phi))
    image = np.where(rho_squared <= 1.0, image, 0.0)

    transformed = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(image)))
    frequencies = np.fft.fftshift(np.fft.fftfreq(n_side, d=step))
    return step ** 2 * transformed, frequencies


def check_full_map_fft(l_max=4, n_side=129, band=16):
    """V15, numerical rung -- the whole map against a 2-D FFT of its image.

    A starspot, rotated so that it has m != 0 content, is fed through the whole chain
    --- real coefficients, both parity families, the rotation --- and its visibility
    compared with a direct transform of the rendered image.  This is the only check
    that would catch an error common to both families, and it is the end-to-end
    anchor for the map.

    Only the lowest frequencies are compared.  The disk's rim is a discontinuity, so
    its Gibbs ringing contaminates the high frequencies at any resolution; the
    comparison covers |p| <= band, well below Nyquist.  The residual is relative to
    the flux and includes the zero-frequency component, which must equal it.
    """
    coeffs = rotate_coefficients(_zonal_starspot(l_max), 0.0, 0.9, 0.4)
    values, frequencies = _fft_visibility(coeffs, n_side)
    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))

    u_grid, v_grid = np.meshgrid(frequencies, frequencies)
    # _fft_visibility renders with x_east as the first meshgrid output, so the FFT's
    # first index belongs to East, which D24 makes v.
    exact = visibility(coeffs, v_grid, u_grid)

    centre = n_side // 2
    worst = abs(values[centre, centre] - exact[centre, centre]) / abs(flux)
    window = slice(centre - band, centre + band + 1)
    worst = max(worst, float(np.max(np.abs(values[window, window]
                                           - exact[window, window]))) / abs(flux))
    return worst


# --------------------------------------------------------------------------
# V19, V20 -- the deformation field
# --------------------------------------------------------------------------

def _sphere_directions(n):
    """n directions spread over the sphere, by the golden-angle spiral."""
    index = np.arange(n) + 0.5
    mu = 1.0 - 2.0 * index / n
    phi = np.pi * (1.0 + 5.0 ** 0.5) * index
    sine = np.sqrt(1.0 - mu ** 2)
    return np.stack([sine * np.cos(phi), sine * np.sin(phi), mu], axis=-1)


def check_deformation_first_order(amplitudes=(1e-4, 1e-3), n_dir=24):
    """V19, numerical rung -- the deformation field against its first-order limit.

    On a level surface delta r = -phi_1/g_0 and delta g = d phi_1/d r, so for small
    amplitudes the solver must reproduce those, with the residual going as eps^2.
    Both quadrupoles are present, perpendicular, and the residual is reported divided
    by eps^2 so that the two amplitudes can share a tolerance.
    """
    worst = 0.0
    for amplitude in amplitudes:
        field = Deformation(tide=amplitude, tide_direction=(1.0, 0.0, 0.0),
                            centrifugal=amplitude, spin_direction=(0.0, 0.0, 1.0))
        for direction in _sphere_directions(n_dir):
            radius = field.level_radius(direction)
            gravity = field.surface_gravity(direction)
            delta_r, delta_g = field.first_order_radii(direction)
            worst = max(worst,
                        abs(radius - (1.0 + delta_r)) / amplitude ** 2,
                        abs(gravity - (1.0 + delta_g)) / amplitude ** 2)
    return worst


def check_level_surface(tide=0.03, centrifugal=0.01, n_dir=32):
    """V20, numerical rung -- the surface is level, and the gradient is analytic.

    Two properties that hold at any amplitude, and so are invisible to V19.  First,
    Phi takes the same value at every point of the surface the solver returns.
    Second, the analytic gradient agrees with a central difference, which is what
    makes g_eff = |grad Phi| trustworthy rather than merely asserted.

    The amplitudes are a tide-dominant pair, since eps_tide + 1.5 eps_cent is what
    has to stay below about 0.059 for the surface to close at all.
    """
    field = Deformation(tide=tide, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=centrifugal, spin_direction=(0.0, 0.0, 1.0))
    worst = 0.0
    for direction in _sphere_directions(n_dir):
        radius = field.level_radius(direction)
        point = radius * direction
        worst = max(worst, abs(float(field.potential(point)) + 1.0))

        step = 1e-6
        difference = np.array([
            (float(field.potential(point + step * axis))
             - float(field.potential(point - step * axis))) / (2.0 * step)
            for axis in np.eye(3)])
        worst = max(worst, float(np.max(np.abs(difference - field.gradient(point)))))
    return worst
