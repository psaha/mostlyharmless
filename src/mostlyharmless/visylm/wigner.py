"""Wigner rotation matrices, and the rotation of a map's coefficients.

Conventions, from section 01 of the manual and D14 of the roadmap:

* The rotation is active and right-handed and is built as

      R(alpha, beta, PA) = R_z(pi/2 - PA) R_y(beta) R_z(alpha) ,

  so that PA is the position angle of the spin axis measured North towards East,
  and beta the inclination measured from the line of sight.

* The harmonic rotation law is

      Y_{lm}(R^-1 n) = sum_{m'} D^l_{m'm}(R) Y_{lm'}(n) ,

  hence a map's coefficients transform as

      b_{lm'} = sum_m D^l_{m'm}(R) a_{lm} .

* Negative orders are indices of the representation, never coefficients (D14).
  The rotation therefore unfolds the stored non-negative half internally, applies
  D, and folds the result back, checking the reality condition as it folds.

The d matrix is real; only the two diagonal phase factors of D are complex, and
they carry the spin phase and the position angle rather than the inclination.
"""

from __future__ import annotations

import math

import numpy as np

from mostlyharmless.visylm.coefficients import Coeffs

__all__ = ["wigner_d", "wigner_D", "rotation_matrix", "rotate_coefficients"]


def wigner_d(l, beta):
    """The real Wigner matrix d^l_{m'm}(beta) = <l m'| exp(-i beta J_y) |l m>.

    Returned with index ``[m' + l, m + l]``.  The sum is the standard explicit
    one, evaluated in floating point; certificate V4 checks it against the matrix
    exponential.
    """
    d = np.zeros((2 * l + 1, 2 * l + 1))
    cb = math.cos(beta / 2.0)
    sb = math.sin(beta / 2.0)
    for i in range(2 * l + 1):
        mp = i - l
        for j in range(2 * l + 1):
            m = j - l
            total = 0.0
            for s in range(max(0, m - mp), min(l + m, l - mp) + 1):
                total += ((-1) ** (mp - m + s)
                          / (math.factorial(l + m - s)
                             * math.factorial(s)
                             * math.factorial(mp - m + s)
                             * math.factorial(l - mp - s))
                          * cb ** (2 * l + m - mp - 2 * s)
                          * sb ** (mp - m + 2 * s))
            d[i, j] = (math.sqrt(math.factorial(l + mp) * math.factorial(l - mp)
                                 * math.factorial(l + m) * math.factorial(l - m))
                       * total)
    return d


def wigner_D(l, alpha, beta, gamma):
    """D^l_{m'm}(alpha, beta, gamma) = exp(-i m' alpha) d(beta) exp(-i m gamma).

    The standard zyz convention, in which the rotation is
    R = R_z(alpha) R_y(beta) R_z(gamma).  Note which angle multiplies which index:
    ``alpha`` goes with the row index m' and ``gamma`` with the column index m.
    For the viewing angles of this package that makes the body-to-observer matrix
    ``wigner_D(l, pi/2 - PA, beta, alpha)``.
    """
    m = np.arange(-l, l + 1)
    return (np.exp(-1j * alpha * m)[:, None]
            * wigner_d(l, beta)
            * np.exp(-1j * gamma * m)[None, :])


def rotation_matrix(alpha, beta, pa):
    """The 3x3 active right-handed body-to-observer rotation R(alpha, beta, PA).

    By D24 the observer frame is (x, y, z) = (North, East, toward the observer), so
    PA is measured from North through East directly and no offset is needed:

        R = R_z(PA) R_y(beta) R_z(alpha) .

    Check it: the body z axis, which is the spin, maps to
    ``(sin(beta) cos(PA), sin(beta) sin(PA), cos(beta))``, so it tilts toward the
    azimuth PA and its distance from the sky plane is ``cos(beta)``.

    Orthogonal, so its transpose is its inverse; the geometric certificates use
    it as the independent reference for the coefficient transform.
    """
    ca, sa = math.cos(alpha), math.sin(alpha)
    cb, sb = math.cos(beta), math.sin(beta)
    cg, sg = math.cos(pa), math.sin(pa)
    Rz_a = np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])
    Ry_b = np.array([[cb, 0.0, sb], [0.0, 1.0, 0.0], [-sb, 0.0, cb]])
    Rz_g = np.array([[cg, -sg, 0.0], [sg, cg, 0.0], [0.0, 0.0, 1.0]])
    return Rz_g @ Ry_b @ Rz_a


_D_MATRICES = {"beta": None, "matrices": None}


def d_matrices(l_max, beta):
    """``wigner_d(l, beta)`` for ``l = 0 .. l_max``, remembered for the current beta.

    The two z-rotations are already phases, so the d-matrix is nearly all of a rotation's
    cost, and it does not depend on alpha or PA.  A single-entry cache is enough for the
    interactive case: a slider drag changes one angle per event, so a drag in alpha reuses
    what the previous event built, and only a drag in beta rebuilds.
    """
    beta = float(beta)
    cached = _D_MATRICES["matrices"]
    if (_D_MATRICES["beta"] != beta or cached is None or len(cached) <= l_max):
        cached = [wigner_d(l, beta) for l in range(l_max + 1)]
        _D_MATRICES["beta"] = beta
        _D_MATRICES["matrices"] = cached
    return cached


def rotate_coefficients(coeffs, alpha, beta, pa):
    """Rotate a map's coefficients from the body frame to the observer frame.

    Implements ``b_{lm'} = sum_m D^l_{m'm}(R) a_{lm}`` with R = R(alpha, beta, PA)
    from :func:`rotation_matrix`, i.e. ``R_z(PA) R_y(beta) R_z(alpha)``, which is
    ``wigner_D(l, PA, beta, alpha)``.  Negative orders appear only as indices here;
    the result is folded back onto the stored non-negative half, and the fold checks
    the reality condition.
    """
    gamma = pa
    inner = alpha
    unfolded = list(coeffs.unfolded())
    matrices = d_matrices(len(unfolded) - 1, beta)
    rotated_blocks = []
    for l, block in enumerate(unfolded):
        m = np.arange(-l, l + 1)
        # the standard d matrix is quoted in the Condon-Shortley basis, in which our
        # Y_lm is (-1)^m times the standard one for m > 0 and equal to it for m <= 0;
        # the extension is plain conjugation, so the change of basis is one-sided
        sign = np.where((m > 0) & (m % 2 == 1), -1.0, 1.0)
        shifted = np.exp(-1j * inner * m) * sign * block
        rotated = matrices[l] @ shifted
        rotated_blocks.append(np.exp(-1j * gamma * m) * sign * rotated)
    return Coeffs.from_unfolded(rotated_blocks)
