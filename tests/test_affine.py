"""Rung 1, the affine model: the stretch, its determinant, and the null case.

The decisive check is that the projected area of the sampled surface equals
``pi |det Lambda|``, the ellipse area the stretch predicts.  That ties the level
radii, the projection onto the sky, the shape matrix and the renderer's weights
together in one number.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import affine, render
from mostlyharmless.visylm.deformation import Deformation

EYE = np.array([0.0, 0.0, 1.0])
TIDE, SPIN = 0.03, 0.01
TRIAXIAL = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                       centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))


def test_the_projected_area_matches_the_stretch_determinant():
    points, _, weights = render.surface_samples(TRIAXIAL, EYE, n_samples=40000)
    stretch = affine.projected_stretch(TRIAXIAL, EYE)
    area = float(np.sum(weights))
    # the deformed surface is an ellipsoid only to O(eps^2), and the quadrature itself
    # carries about half a per cent at this sampling
    assert area == pytest.approx(np.pi * abs(np.linalg.det(stretch)), rel=3e-2)


def test_a_sphere_has_unit_stretch_and_so_affine_agrees_with_the_sphere():
    """With no deformation the axes are all one, Lambda is the identity, and the affine
    visibility is the sphere's own -- so the model is exact where it must be."""
    field = Deformation()
    stretch = affine.projected_stretch(field, EYE)
    assert np.allclose(stretch, np.eye(2), atol=1e-12)

    points, brightness, weights = render.surface_samples(field, EYE, n_samples=40000)
    east, north = render.sky_basis(EYE)

    def sphere_visibility(u, v):
        return render.direct_visibility(points, brightness, weights, east, north, u, v)

    for u, v in ((0.0, 0.0), (0.7, 0.0), (0.3, 0.9)):
        value = affine.affine_visibility(sphere_visibility, stretch, u, v)
        assert value == pytest.approx(sphere_visibility(u, v), rel=1e-12)


def test_a_sphere_projects_to_the_unit_disk_from_any_direction():
    """A sanity check on the shape matrix: a sphere projects to the unit disk from any
    line of sight, so the determinant is one whatever the direction."""
    field = Deformation()
    for eye in ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.3, 0.5, 0.81)):
        stretch = affine.projected_stretch(field, np.array(eye))
        assert abs(np.linalg.det(stretch)) == pytest.approx(1.0, abs=1e-12)
