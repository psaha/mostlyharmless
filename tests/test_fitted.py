"""The fitted rung: the disk harmonics, the Gram matrix, and the fit itself.

The row is a diagnostic, so the tests are about *what the numbers mean*: that under
truncation the fit is exact where it should be, that its residual is the truncation floor
rather than a defect, and that the two routes to that residual agree.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.symbolic import fitted, render
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk

EYE = np.array([np.sin(np.radians(60.0)), 0.0, np.cos(np.radians(60.0))])
# No deformation: the level surface is the unit sphere, so Lambda is the identity and the
# fitted coefficients are the pattern's own -- which is what the round trip is about.
SPHERE = Deformation(tide=0.0, tide_direction=(1.0, 0.0, 0.0),
                     centrifugal=0.0, spin_direction=(0.0, 0.0, 1.0))
TIDE = 0.096
SPIN = 0.032


def _fit(field, eye=EYE, samples=20000, l_max=4):
    east, north = render.sky_basis(eye)
    star = render.surface_samples(field, eye, n_samples=samples)
    from mostlyharmless.symbolic import affine
    stretch = affine.projected_stretch(field, eye)
    return fitted.fit_image(star, stretch, east, north, l_max=l_max), stretch, (east, north)


def test_the_disk_harmonic_of_a_uniform_star_is_constant():
    """``Phi_00`` is 1 on the disk and vanishes outside it: the uniform disk."""
    east_coordinate = np.array([0.0, 0.3, -0.7, 0.99])
    north_coordinate = np.array([0.0, -0.3, 0.5, 0.0])
    values = fitted.disk_harmonic(0, 0, east_coordinate, north_coordinate)
    assert np.allclose(values, 1.0)

    outside = fitted.disk_harmonic(0, 0, np.array([1.2]), np.array([0.0]))
    assert np.allclose(outside, 0.0)


def test_the_gram_matrix_is_anchored_and_block_diagonal():
    """The ``(0, 0)`` entry is the disk's area, and orders do not mix."""
    matrix = fitted.gram(4)
    assert abs(matrix[0, 0] - np.pi) < 1e-12   # int over the disk of Phi_00^2 = 1

    orders = np.array([m for _, m, _ in fitted.design(4)])
    between = matrix[np.ix_(orders != 0, orders == 0)]
    assert np.max(np.abs(between)) < 1e-12

    condition = fitted.block_condition(4)
    assert len(condition) == 5
    assert all(np.isfinite(c) and c > 0.0 for c in condition)


def test_the_fitted_row_recovers_a_uniform_disk():
    """With no deformation the fit must give ``a_00 = 1`` and nothing else.

    That is the same statement as the uniform-disk anchor: the fitted visibility must be
    the classical ``2 J_1(x)/x`` with flux ``pi``.
    """
    fit, stretch, (east, north) = _fit(SPHERE)
    assert np.allclose(stretch, np.eye(2), atol=1e-12)

    coefficients = fit.coefficients
    assert abs(coefficients.blocks[0][0] - 1.0) < 2e-3
    for l in range(0, 5, 2):
        for m in range(l + 1):
            if (l, m) != (0, 0):
                assert abs(coefficients.blocks[l][m]) < 2e-3

    for u, v in ((0.3, 0.0), (0.0, 0.4), (0.5, 0.5)):
        modelled = visibility(coefficients, u, v)
        anchor = visibility_uniform_disk(u, v, 1.0)
        assert abs(modelled / visibility(coefficients, 0.0, 0.0)
                   - anchor) < 3e-3


def test_the_two_residual_routes_agree():
    """``||I'||^2 - b^dag G^-1 b`` against the direct expansion of the same quadratic.

    A slip of a factor of ``|det Lambda|`` would leave the coefficients plausible and show
    up only here.
    """
    field = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
    fit, _, _ = _fit(field)
    assert fit.residual > 0.0
    assert abs(fit.residual - fit.residual_direct) < 1e-10 * fit.norm


def test_the_residual_falls_with_the_mode_count():
    """The residual is truncation: adding ``l = 4`` must reduce it, on the real star."""
    field = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
    two, _, _ = _fit(field, l_max=2)
    four, _, _ = _fit(field, l_max=4)
    print(f"\nfitted residual: l<=2 {two.residual:.6f}, l<=4 {four.residual:.6f}")
    assert four.residual < two.residual


def test_the_fitted_image_carries_the_right_flux():
    """The fitted image's flux must match the data's, since ``l = 0`` is fitted too.

    A weaker statement than reproducing the image point by point, and honest about the
    rasterisation: the model image is evaluated on a grid, so the integral carries the
    grid's error and the tolerance says so.
    """
    field = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
    fit, stretch, (east, north) = _fit(field, samples=20000)
    star = render.surface_samples(field, EYE, n_samples=20000)
    _, brightness, weights = star

    image = fitted.model_image(fit.coefficients, stretch, half_width=1.5, n=241)
    cell = (3.0 / 240.0) ** 2
    model_flux = float(np.nansum(image)) * cell       # the rim is NaN, like sky_image
    numeric_flux = float(np.sum(weights * brightness))
    print(f"\nflux: numeric {numeric_flux:.6f}, fitted image {model_flux:.6f}")
    assert abs(model_flux - numeric_flux) < 0.02 * abs(numeric_flux)


def test_the_fitted_round_trip():
    """V26: a synthetic image, fitted back to its own coefficients."""
    assert dv.check_fitted_round_trip() < 1e-10


def test_the_disk_gram_against_brute_force():
    """V27: an analytic entry, the block structure, and a midpoint rule."""
    assert dv.check_disk_gram() < 5e-3
