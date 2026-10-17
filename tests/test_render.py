"""The renderer: the surface, its projection, and the direct transform.

The decisive check is the null case.  A zero deformation is the unit sphere with
|grad Phi| = 1 everywhere, so the level solve, the normals, the projection factor and
the quadrature must together reproduce the uniform-disk visibility 2 J_1(x)/x, and
the flux must come out as pi.
"""

import numpy as np
import pytest
from scipy.special import jv

from mostlyharmless.symbolic import render
from mostlyharmless.visylm.deformation import Deformation

EYE = np.array([0.0, 0.0, 1.0])


def samples(field, n_samples=40000):
    east, north = render.sky_basis(EYE)
    points, brightness, weights = render.surface_samples(field, EYE,
                                                         n_samples=n_samples)
    return points, brightness, weights, east, north


def test_the_sphere_renders_as_the_uniform_disk():
    points, brightness, weights, east, north = samples(Deformation())
    flux = render.direct_visibility(points, brightness, weights, east, north, 0.0, 0.0)
    assert flux.real == pytest.approx(np.pi, rel=1e-2)

    for q in (0.5, 1.0, 2.0):
        value = render.direct_visibility(points, brightness, weights, east, north,
                                         q, 0.0)
        x = 2.0 * np.pi * q
        assert (value / flux).real == pytest.approx(2.0 * jv(1, x) / x, abs=1e-2)


def test_the_sphere_is_centrally_symmetric_so_its_transform_is_real():
    """The imaginary part is sampling noise here, which is a check on the projection:
    a real, even brightness distribution has a real transform."""
    points, brightness, weights, east, north = samples(Deformation())
    value = render.direct_visibility(points, brightness, weights, east, north, 1.3, 0.7)
    assert abs(value.imag) < 0.02 * abs(value)


def test_gravity_darkening_makes_the_pole_brighter_than_the_equator():
    """One spin quadrupole: at the pole the quadrupole cancels the l = 0 piece, so g
    is unchanged, while at the equator g is reduced.  Brighter pole, darker equator,
    which is the standard signature."""
    field = Deformation(centrifugal=0.05)
    points, brightness, weights, east, north = samples(field)
    radii = np.linalg.norm(points, axis=1)
    directions = points / radii[:, None]

    polar = directions[:, 2] > 0.99
    equatorial = np.abs(directions[:, 2]) < 0.01
    assert np.any(polar) and np.any(equatorial)
    assert np.mean(brightness[polar]) > np.mean(brightness[equatorial])


def test_an_open_surface_is_reported_rather_than_sampled():
    """Past the validity limit the level surface does not close, and the renderer must
    say so rather than bin nonsense."""
    field = Deformation(tide=0.2, tide_direction=(1.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        render.surface_samples(field, EYE, n_samples=2000)
