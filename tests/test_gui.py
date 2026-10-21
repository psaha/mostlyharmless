"""The GUI's compute path, tested without any widget.

These are the pieces the gate names: the closure bound in closed form, the clamps that
use it, and the two panel builders. The widgets themselves are not exercised.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import gui


def test_the_depth_reproduces_the_equatorial_combination():
    """At latitude zero, along the tide axis, the depth is eps_tide + 1.5 eps_spin."""
    depth = gui.closure_depth(tide=0.06, spin=0.02, latitude=0.0)
    assert depth == pytest.approx(0.06 + 1.5 * 0.02, rel=1e-3)


def test_a_polar_tide_is_limited_only_by_its_own_amplitude():
    """Along the spin the centrifugal quadrupole has no gradient, so it cannot deepen."""
    depth = gui.closure_depth(tide=0.06, spin=0.02, latitude=np.pi / 2.0)
    assert depth == pytest.approx(0.06, rel=1e-3)


def test_the_largest_spin_keeps_the_depth_at_the_limit():
    tide, latitude = 0.04, 0.3
    spin = gui.largest_spin(tide, latitude)
    assert gui.closure_depth(tide, spin, latitude) == pytest.approx(gui.CLOSURE,
                                                                    rel=1e-3)


def test_the_smallest_latitude_keeps_the_depth_at_the_limit():
    """Infeasible at the equator, so the bound bites and the depth lands on the limit."""
    tide, spin = 0.12, 0.03
    assert gui.closure_depth(tide, spin, 0.0) > gui.CLOSURE
    latitude = gui.smallest_latitude(tide, spin)
    assert latitude > 0.0
    assert gui.closure_depth(tide, spin, latitude) == pytest.approx(gui.CLOSURE,
                                                                    rel=1e-3)


def test_clamping_returns_a_feasible_point():
    """An amplitude pair that cannot close is brought inside, by latitude and spin."""
    parameters = gui.Parameters(tide=0.12, spin=0.05, latitude=0.0)
    fixed = gui.clamped(parameters)
    assert gui.closure_depth(fixed.tide, fixed.spin, abs(fixed.latitude)) <= (
        gui.CLOSURE * (1.0 + 1e-3))
    assert abs(fixed.latitude) > 0.0

    feasible = gui.Parameters(tide=0.02, spin=0.005, latitude=0.0)
    assert gui.clamped(feasible) == feasible


def test_a_feasible_default_is_left_alone():
    assert gui.clamped(gui.Parameters()) == gui.Parameters()


def test_compute_image_returns_two_panels():
    panels = gui.compute_image(gui.Parameters(), "numeric", samples=2000, pixels=41)
    assert len(panels) == 2
    assert panels[0].shape == panels[1].shape
    # mean-removed, so the average over the disk is zero
    assert np.nanmean(panels[0]) == pytest.approx(0.0, abs=1e-12)


def test_compute_visibility_vanishes_at_the_origin():
    """The unit-sphere reference is subtracted, so D(0,0) = 0."""
    values = gui.compute_visibility(gui.Parameters(), "numeric", samples=2000,
                                    grid=15, q_max=1.0)
    centre = 15 // 2
    assert abs(values[centre, centre]) == pytest.approx(0.0, abs=2e-3)


def test_all_three_models_produce_finite_panels():
    parameters = gui.Parameters(tide=0.03, spin=0.01)
    for model in ("numeric", "affine", "first-order"):
        values = gui.compute_visibility(parameters, model, samples=3000, grid=9,
                                        q_max=1.0)
        assert np.all(np.isfinite(values))
        panels = gui.compute_image(parameters, model, samples=3000, pixels=31)
        assert len(panels) == 2


def test_the_analytic_and_numeric_panels_are_similar():
    """First order in eps and l <= 2, so the two models should agree closely."""
    parameters = gui.Parameters(tide=0.03, spin=0.01)
    analytic = gui.compute_visibility(parameters, "affine", samples=4000,
                                      grid=11, q_max=1.0)
    numeric = gui.compute_visibility(parameters, "numeric", samples=4000,
                                     grid=11, q_max=1.0)
    assert np.max(np.abs(analytic - numeric)) < 0.15


def test_the_linear_pattern_is_a_second_order_change():
    """Equation (74) differs from the exact pattern at O(eps^2), not at O(eps).

    Both models are first order in the deformation, so the difference between them is
    second order: halving the amplitudes should quarter it.  A factor of two would mean
    the linearisation had kept only one of its two equal first-order contributions --- the
    radial derivative without the surface's displacement.
    """
    def spread(scale):
        parameters = gui.Parameters(tide=0.02 * scale, spin=0.01 * scale)
        exact = gui.compute_visibility(parameters, "affine", samples=4000, grid=7,
                                       q_max=1.0)
        linear = gui.compute_visibility(parameters, "first-order", samples=4000,
                                        grid=7, q_max=1.0)
        return float(np.max(np.abs(exact - linear)))

    full, half = spread(1.0), spread(0.5)
    print(f"\nlinear-pattern spread: {full:.6f} {half:.6f}  ratio {half / max(full, 1e-30):.3f}")
    assert half < 0.4 * full
