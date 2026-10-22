"""The fitted rung: Lambda from the level surface, coefficients fitted.

Manual section 8.8.  The model class is the affine row's --- ``V = sum_lm A_lm
V_lm(Lambda^T u)`` with ``Lambda`` the projected stretch --- but the coefficients are
chosen by least squares in the image plane instead of being taken from the brightness
pattern's projection.  At fixed ``Lambda`` the residual is the truncation floor *given*
the ellipsoidal shape, which is what the model's fitted row reports: comparing it with
the affine row's error separates the two errors that row conflates.

The whole fit reduces to sums over the samples the renderer already produces.  No image
grid is built and nothing is interpolated: the basis's Gram matrix carries no data and is
computed once, and the right-hand side is one pass over the samples.
"""

__all__ = ["disk_harmonic", "design", "gram", "block_condition", "Fit", "fit_image",
           "residual_of", "model_image"]

import functools
from dataclasses import dataclass

import numpy as np

from ..visylm.coefficients import Coeffs
from ..visylm.spherical_harmonics import sph_harm


def _lift(east_coordinate, north_coordinate):
    """The visible direction at disk points, in the observer frame of (77).

    Sky coordinates are ``(east, north)`` --- the order Lambda is stored in --- and the
    azimuth is measured from North towards East, as the harmonics of section 2 are.
    """
    rho_squared = east_coordinate**2 + north_coordinate**2
    z = np.sqrt(np.clip(1.0 - rho_squared, 0.0, None))
    return np.arccos(np.clip(z, -1.0, 1.0)), np.arctan2(east_coordinate, north_coordinate)


def disk_harmonic(l, m, east_coordinate, north_coordinate):
    """``Phi_lm``, the projection of ``Y_l^m`` onto the disk, in closed form (eq. 89).

    The projection *is* the visible surface value: the opaque convention absorbs the
    viewing cosine into the area element, so there is none here.  For a uniform star
    ``Phi_00`` is constant, the image is the uniform disk of flux ``pi``, and the
    visibility is the classical anchor --- which is what ``test_the_fitted_row_matches_the
    _uniform_disk`` checks.

    Off the unit disk it is zero, and that must be explicit: the clip confines
    ``arccos`` to a real angle, but it does not make the harmonic vanish there --- at
    ``rho > 1`` the lifted direction is the limb direction, and ``Y_l^m`` of it is a
    perfectly ordinary number.  The star has no support beyond the silhouette, so the
    mask is part of the definition and not a convenience.
    """
    theta, phi = _lift(east_coordinate, north_coordinate)
    inside = east_coordinate**2 + north_coordinate**2 <= 1.0
    return np.where(inside, sph_harm(l, m, theta, phi), 0.0)


def design(l_max=4):
    """The real design functions, grouped by order.

    For ``l_max = 4`` this is fifteen parameters, matching ``1 + 5 + 9``: one per degree
    at ``m = 0``, and two per degree for each ``m > 0``, carrying the real and imaginary
    parts of ``A_lm``.  Odd degrees are absent because the parity argument of section 8.4
    forbids them, not for convenience.

    Entries are ``(l, m, part)`` with ``part`` in ``{"re", "im"}``.  The list is
    contiguous in ``m``, so the Gram matrix's blocks are slices of it.
    """
    even = list(range(0, l_max + 1, 2))
    plan = [(l, 0, "re") for l in even]
    for m in range(1, l_max + 1):
        for l in even:
            if l >= m:
                plan.append((l, m, "re"))
                plan.append((l, m, "im"))
    return plan


def _evaluate(plan, east_coordinate, north_coordinate):
    """The design functions at disk points, shape ``(len(plan),) + shape``.

    The pairing is the harmonics' own, ``Y_l^{-m} = (-1)^m conj(Y_l^m)``, so the
    contribution of order ``m`` to the real image is ``2 Re(A_lm Phi_lm)``: the two real
    columns are ``2 Re Phi`` and ``-2 Im Phi``.  That is *not* the visibility's pairing of
    section 10, and the two must not be confused.
    """
    columns = []
    for l, m, part in plan:
        phi = disk_harmonic(l, m, east_coordinate, north_coordinate)
        if m == 0:
            columns.append(np.asarray(phi.real, dtype=float))
        elif part == "re":
            columns.append(2.0 * np.asarray(phi.real, dtype=float))
        else:
            columns.append(-2.0 * np.asarray(phi.imag, dtype=float))
    return np.array(columns)


@functools.lru_cache(maxsize=None)
def _gram_and_blocks(l_max=4, n_theta=96, n_phi=128):
    """The design functions' Gram matrix on the disk, and the per-order index blocks."""
    plan = design(l_max)
    nodes, node_weights = np.polynomial.legendre.leggauss(n_theta)
    theta = 0.25 * np.pi * (nodes + 1.0)                    # [0, pi/2]
    w_theta = 0.25 * np.pi * node_weights
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi            # periodic: trapezoid is exact

    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    rho = np.sin(theta_grid)
    east_coordinate = rho * np.cos(phi_grid)
    north_coordinate = rho * np.sin(phi_grid)
    measure = np.outer(w_theta * np.sin(theta) * np.cos(theta),
                       np.full(n_phi, 2.0 * np.pi / n_phi))

    columns = _evaluate(plan, east_coordinate, north_coordinate)
    matrix = np.einsum("ipq,jpq,pq->ij", columns, columns, measure)

    blocks = {}
    for index, (_, m, _) in enumerate(plan):
        blocks.setdefault(m, []).append(index)
    return matrix, tuple((m, np.array(indices)) for m, indices in sorted(blocks.items()))


def gram(l_max=4):
    """The design functions' Gram matrix; it carries no data, so it is computed once."""
    return _gram_and_blocks(l_max)[0]


def block_condition(l_max=4):
    """The per-order condition numbers of the Gram matrix (section 8.8)."""
    matrix = gram(l_max)
    return tuple(float(np.linalg.cond(matrix[np.ix_(indices, indices)]))
                 for _, indices in _gram_and_blocks(l_max)[1])


@dataclass(frozen=True)
class Fit:
    """What a fitted row reports: coefficients, residual, and the conditioning."""

    coefficients: Coeffs
    residual: float
    residual_direct: float
    norm: float
    condition: tuple


def _as_coeffs(plan, values, l_max):
    """The fitted real parameters as the standard half-complex coefficient blocks."""
    blocks = [np.zeros(l + 1, dtype=complex) for l in range(l_max + 1)]
    for (l, m, part), value in zip(plan, values):
        if m == 0:
            blocks[l][0] = complex(value, 0.0)
        elif part == "re":
            blocks[l][m] = complex(value, blocks[l][m].imag)
        else:
            blocks[l][m] = complex(blocks[l][m].real, value)
    return Coeffs(tuple(blocks))


def fit_image(samples, stretch, east, north, l_max=4):
    """Fit the coefficients to rendered samples (manual section 8.8).

    ``samples`` is the triple from :func:`render.surface_samples`, and ``stretch`` the
    2x2 ``Lambda`` in the sky's ``(east, north)`` order.  The right-hand side is a sum
    over the samples themselves, so the cost is one 2x2 solve and fifteen mode
    evaluations per sample.

    The residual comes by two routes.  ``residual`` eliminates the normal equations;
    ``residual_direct`` does not.  They must agree --- a slip of a factor of
    ``|det Lambda|`` would not show in the coefficients and would show here --- and
    :func:`check_fitted_round_trip` holds them to that at machine precision.
    """
    points, brightness, weights = samples
    sky = np.stack([points @ east, points @ north], axis=-1)
    disk = sky @ np.linalg.inv(stretch).T                  # z = Lambda^-1 x

    plan = design(l_max)
    columns = _evaluate(plan, disk[:, 0], disk[:, 1])
    weighted = weights * brightness
    rhs = columns @ weighted

    matrix = gram(l_max)
    values = np.linalg.solve(matrix, rhs)

    scale = abs(float(np.linalg.det(stretch)))
    norm = scale * float(np.sum(weights * brightness**2))
    residual = norm - float(rhs @ values)
    residual_direct = (norm - 2.0 * float(rhs @ values)
                       + float(values @ matrix @ values))
    return Fit(_as_coeffs(plan, values, l_max), residual, residual_direct, norm,
               block_condition(l_max))


def model_image(coefficients, stretch, half_width=1.5, n=241):
    """The fitted model's image, indexed ``[north, east]`` like ``sky_image``.

    ``I_model(x) = F(Lambda^-1 x) / |det Lambda|`` with ``F`` the fitted series, which is
    eq. (88).  Evaluating it on a grid is all that the display needs: the visibility
    itself comes from the mode sum, not from this array.
    """
    axis = np.linspace(-half_width, half_width, n)
    east_coordinate = np.broadcast_to(axis[None, :], (n, n))
    north_coordinate = np.broadcast_to(axis[:, None], (n, n))

    inverse = np.linalg.inv(stretch)
    disk_east = inverse[0, 0] * east_coordinate + inverse[0, 1] * north_coordinate
    disk_north = inverse[1, 0] * east_coordinate + inverse[1, 1] * north_coordinate

    image = np.zeros((n, n))
    for l, block in enumerate(coefficients.blocks):
        image = image + block[0].real * disk_harmonic(l, 0, disk_east, disk_north).real
        for m in range(1, l + 1):
            phi = disk_harmonic(l, m, disk_east, disk_north)
            image = image + 2.0 * (block[m] * phi).real
    image = image / abs(float(np.linalg.det(stretch)))

    # Outside the silhouette, NaN rather than zero: that is what ``sky_image`` returns, and
    # the display divides by the *mean of the finite pixels*.  A zero rim would be counted,
    # diluting the mean by the disk-to-square area ratio and amplifying the contrast by
    # roughly three --- which looks like a plausible pattern in the wrong colours.
    silhouette = disk_east**2 + disk_north**2 > 1.0
    return np.where(silhouette, np.nan, image)


def residual_of(coefficients, samples, stretch, east, north, l_max=4):
    """The image-space residual of *any* coefficient set, in the fitted row's metric.

    Needed to compare the fitted row with a projected one on equal terms: the fit is
    optimal in this norm, so its residual must not exceed the projection's when both live
    in the same span.  A coefficient set of lower truncation is padded with zeros, which is
    what "the same span" means.  Returns ``(residual, norm)``.
    """
    points, brightness, weights = samples
    sky = np.stack([points @ east, points @ north], axis=-1)
    disk = sky @ np.linalg.inv(stretch).T

    plan = design(l_max)
    values = _parameters_of(coefficients, plan)
    columns = _evaluate(plan, disk[:, 0], disk[:, 1])
    rhs = columns @ (weights * brightness)

    matrix = gram(l_max)
    scale = abs(float(np.linalg.det(stretch)))
    norm = scale * float(np.sum(weights * brightness**2))
    residual = (norm - 2.0 * float(values @ rhs)
                + float(values @ matrix @ values))
    return residual, norm


def _parameters_of(coefficients, plan):
    """The real parameter vector a coefficient set corresponds to, in the design order."""
    values = []
    for l, m, part in plan:
        block = 0.0
        if l < len(coefficients.blocks) and m < len(coefficients.blocks[l]):
            block = coefficients.blocks[l][m]
        values.append(float(block.real if (m == 0 or part == "re") else block.imag))
    return np.array(values)
