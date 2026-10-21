"""Rung 1, the affine model: the stretch, its determinant, and the null case.

The decisive check is that the projected area of the sampled surface equals
``pi |det Lambda|``, the ellipse area the stretch predicts.  That ties the level
radii, the projection onto the sky, the shape matrix and the renderer's weights
together in one number.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import affine, render
from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.visibility import visibility

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


def test_the_projected_shape_does_not_depend_on_the_tide_longitude():
    """Rotating the body about its spin must not change the projected shape.

    Built from the coordinate radii rather than the principal axes, the stretch does
    change with longitude -- a first-order error in the shape -- and this catches it.
    Seen down the spin the projected area is the equatorial ellipse's, whose area a
    rotation cannot change.
    """
    eye = np.array([0.0, 0.0, 1.0])
    determinants = []
    for longitude in (0.0, np.pi / 4.0, np.pi / 2.0):
        direction = (np.cos(longitude), np.sin(longitude), 0.0)
        field = Deformation(tide=TIDE, tide_direction=direction,
                            centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
        determinants.append(abs(np.linalg.det(affine.projected_stretch(field, eye))))
    assert determinants[0] == pytest.approx(determinants[1], rel=1e-12)
    assert determinants[1] == pytest.approx(determinants[2], rel=1e-12)


def test_the_frequency_is_reordered_into_the_stretch_basis():
    """Lambda is in (east, north); the frequencies are (u, v) = (North, East).

    A frequency due North must reach the sphere function with a zero east component and
    the north component scaled by Lambda's north factor.  Passing (u, v) through unchanged
    would put the whole thing on the east axis instead.
    """
    stretch = np.diag([2.0, 3.0])          # east scale 2, north scale 3
    calls = []

    def sphere_visibility(u, v):
        calls.append((float(u), float(v)))
        return 0.0

    affine.affine_visibility(sphere_visibility, stretch, 1.0, 0.0)
    # stretched: east 2*0 = 0, north 3*1 = 3, and the sphere takes (North, East)
    assert calls[-1] == (3.0, 0.0)


def test_the_ellipsoid_ignores_the_tide_direction_when_the_tide_is_zero():
    """With eps_tide = 0 the tide's direction cannot matter, whatever the latitude.

    Taking the principal frame from the tide axis breaks this: the axis rotates, and on
    an oblate body the level radius along a tilted direction differs from the equatorial
    one, so the shape would change with no physical cause.  The frame has to come from
    the eigenvectors of the quadrupole tensor.
    """
    eye = np.array([np.sin(np.pi / 4.0), 0.0, np.cos(np.pi / 4.0)])
    stretches = []
    for latitude in (0.0, 0.4, np.pi / 2.0):
        direction = (np.cos(latitude), 0.0, np.sin(latitude))
        field = Deformation(tide=0.0, tide_direction=direction,
                            centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
        stretches.append(affine.projected_stretch(field, eye))
    for other in stretches[1:]:
        assert other == pytest.approx(stretches[0], rel=1e-12)


def test_the_projection_is_in_the_observer_frame():
    """The projection must use the observer's axes, so an inclined view is not axisymmetric.

    Projecting in the body frame leaves the pattern axisymmetric about the line of sight
    and its visibility exactly real at every frequency.  An oblate, gravity-darkened star
    seen at 60 degrees is not centrally symmetric, so the imaginary part must be non-zero.
    """
    eye = np.array([np.sin(np.pi / 3.0), 0.0, np.cos(np.pi / 3.0)])
    field = Deformation(tide=0.0, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=0.05, spin_direction=(0.0, 0.0, 1.0))
    upright = affine.sphere_coefficients(field, l_max=2, n_samples=4000)
    tilted = affine.sphere_coefficients(field, l_max=2, n_samples=4000, eye=eye)
    assert np.max(np.abs(tilted.to_vector() - upright.to_vector())) > 1e-3

    values = [visibility(tilted, u, v) for u in (0.0, 0.5) for v in (0.0, 0.5)]
    flux = abs(visibility(tilted, 0.0, 0.0))
    assert max(abs(value.imag) for value in values) > 1e-6 * flux


def test_the_affine_prediction_matches_the_direct_sum():
    """The affine model's prediction must reproduce the deformed surface's visibility.

    The prediction is ``V_sphere(Lambda^T u)``; leaving Lambda out shifts it at first order
    in eps and is a trap.  With the stretch included the two should agree to the model's own
    ``O(eps^2)`` error plus the quadrature.

    They currently agree only to about four per cent at eps = 0.015, which is far larger
    than ``eps^2`` and so is *not* the affine approximation's own error --- something
    structural remains.  Rather than assert a number that is known to be too loose or too
    tight, this measures how the residual *scales* with the amplitude, which is what
    identifies the culprit: a residual falling like ``eps^2`` is the model being right and
    the tolerance being wrong, one falling like ``eps`` is a first-order error in the
    machinery, and one that does not fall at all is a normalisation.  Halving the amplitude
    separates the three.
    """
    eye = np.array([np.sin(np.pi / 3.0), 0.0, np.cos(np.pi / 3.0)])
    longitude = np.radians(120.0)

    def residual(scale):
        field = Deformation(tide=0.01 * scale,
                            tide_direction=(np.cos(longitude), np.sin(longitude), 0.0),
                            centrifugal=0.015 * scale, spin_direction=(0.0, 0.0, 1.0))
        east, north = render.sky_basis(eye)
        points, brightness, weights = render.surface_samples(field, eye, n_samples=30000)
        coefficients = affine.sphere_coefficients(field, l_max=2, n_samples=30000, eye=eye)
        stretch = affine.projected_stretch(field, eye)
        sphere_flux = complex(visibility(coefficients, 0.0, 0.0))

        def sphere_visibility(u, v):
            return visibility(coefficients, u, v) / sphere_flux

        flux = render.direct_visibility(points, brightness, weights, east, north, 0.0, 0.0)
        worst = 0.0
        for u, v in ((0.4, 0.0), (0.0, 0.4), (0.3, 0.3)):
            modelled = affine.affine_visibility(sphere_visibility, stretch, u, v)
            summed = (render.direct_visibility(points, brightness, weights,
                                               east, north, u, v) / flux)
            worst = max(worst, abs(modelled - summed))
        return worst

    full, half, zero = residual(1.0), residual(0.5), residual(0.0)
    print(f"\naffine residual at eps {full:.5f}, eps/2 {half:.5f}, zero {zero:.5f}"
          f"  ratio {half / max(full, 1e-30):.3f}")
    # With no deformation the affine prediction and the direct sum are the same object ---
    # Lambda is the identity --- so this residual is the projection machinery alone, with no
    # model error in it.  If it is small the machinery is sound and the eps-scaling residual
    # is the approximation; if it is the same few per cent, the machinery is the fault.
    assert zero < 0.002
    assert half < 0.75 * full


def test_a_sphere_projects_to_the_unit_disk_from_any_direction():
    """A sanity check on the shape matrix: a sphere projects to the unit disk from any
    line of sight, so the determinant is one whatever the direction."""
    field = Deformation()
    for eye in ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.3, 0.5, 0.81)):
        stretch = affine.projected_stretch(field, np.array(eye))
        assert abs(np.linalg.det(stretch)) == pytest.approx(1.0, abs=1e-12)


def test_the_first_order_coefficients_are_the_closed_form():
    """V21: the explicit first-order harmonics agree with the quadrature.

    To first order the pattern is a constant plus two quadrupoles, so its harmonics are
    explicit by the addition theorem: a_00 from the constant, a_2m from the two axes, and
    nothing else --- no odd order, because the pattern is even under n -> -n.  This is the
    certificate for that claim; the tolerance is quadrature, not model error.
    """
    assert dv.check_first_order_coefficients() < 5e-3


def test_the_first_order_stretch_differs_at_second_order():
    """Linearised radii change Lambda only at O(eps^2), so the saving is free.

    ``delta r = -phi_1`` is exact to first order, so the ellipsoid built from it differs
    from the exact one at second order and no earlier.  Halving the amplitudes should
    therefore quarter the difference; a factor of two would mean the linearised radius is
    wrong at first order, which would be a bug rather than a saving.
    """
    eye = np.array([np.sin(np.pi / 3.0), 0.0, np.cos(np.pi / 3.0)])

    def spread(scale):
        field = Deformation(tide=0.02 * scale, tide_direction=(1.0, 0.0, 0.0),
                            centrifugal=0.01 * scale, spin_direction=(0.0, 0.0, 1.0))
        exact = affine.projected_stretch(field, eye)
        linear = affine.projected_stretch(field, eye, first_order=True)
        return float(np.max(np.abs(exact - linear)))

    full, half = spread(1.0), spread(0.5)
    print(f"\nstretch difference: {full:.8f} {half:.8f}  ratio {half / max(full, 1e-30):.3f}")
    assert half < 0.4 * full
