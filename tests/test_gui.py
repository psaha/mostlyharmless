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


def test_every_model_label_computes():
    """Each label in the radio produces finite panels and maps."""
    parameters = gui.Parameters(tide=0.03, spin=0.01)
    assert sum(label.startswith("fitted") for label in gui.MODEL_LABELS) == 2
    assert "affine" not in gui.MODEL_LABELS        # dropped from the interface
    for model in gui.MODEL_LABELS:
        values = gui.compute_visibility(parameters, model, samples=3000, grid=9,
                                        q_max=1.0)
        assert np.all(np.isfinite(values))
        panels = gui.compute_image(parameters, model, samples=3000, pixels=31)
        assert len(panels) == 2


def test_the_fitted_row_beats_the_projected_one_in_its_own_norm():
    """The fit is optimal in the image norm, so no projection can do better there.

    Both live in the same span: the affine row's ``l <= 2`` coefficients, padded with
    zeros to ``l <= 4``, against the fit's.  This is the property that makes the row a
    diagnostic rather than another approximation.
    """
    from mostlyharmless.symbolic import affine, fitted, render

    parameters = gui.Parameters(tide=0.05, spin=0.02)
    field = parameters.field()
    eye, (east, north) = gui._view(parameters.inclination)
    stretch = affine.projected_stretch(field, eye)
    star = render.surface_samples(field, eye, n_samples=8000)

    projected = affine.sphere_coefficients(field, l_max=2, n_samples=8000, eye=eye)
    projected_residual, norm = fitted.residual_of(projected, star, stretch, east, north)
    fit = fitted.fit_image(star, stretch, east, north, l_max=4)

    print(f"\nimage residual: projected {projected_residual:.3e}, "
          f"fitted {fit.residual:.3e}, norm {norm:.3e}")
    assert fit.residual <= projected_residual + 1e-9 * norm


def test_the_fitted_report_agrees_with_the_row_it_reports():
    """The headless accessor behind the fitted row's title.

    Two routes to the same residual, and a residual small enough to be truncation rather
    than a modelling failure at the sliders' default amplitudes.
    """
    parameters = gui.Parameters(tide=0.03, spin=0.01)
    report = gui.fitted_report(parameters, samples=6000)
    assert 0.0 < report.residual < 0.2
    assert abs(report.residual - report.residual_direct) < 1e-9 * report.norm
    assert len(report.condition) == 5


def test_every_title_is_mathtext_the_backend_can_parse():
    """A bad symbol in a title crashes the *draw*, not a callback, and `build` needs a display.

    ``\\le`` is LaTeX but not mathtext, and that surfaced as a traceback from inside the
    renderer with nothing pointing at the string.  Parsing the titles here makes the fault
    a test failure instead.
    """
    pytest.importorskip("matplotlib")
    from matplotlib.mathtext import MathTextParser

    parser = MathTextParser("agg")
    for title in list(gui.MODEL_TITLES.values()) + list(gui.VISIBILITY_TITLES):
        parser.parse(title)                      # raises on an unknown symbol


def test_the_imaginary_part_keeps_its_sign_between_the_rows():
    """A conjugate looks plausible and is wrong, so compare Im where it is largest.

    The mode sum and the direct sum are independent routes to the same complex number, so
    their agreement where ``|Im Delta V|`` is a sizeable fraction of the map's range is
    what rules out a sign flip; an O(eps) discrepancy is the model, not the convention.
    """
    parameters = gui.Parameters(tide=0.03, spin=0.01)
    affine = gui.compute_visibility(parameters, "affine", samples=8000, grid=9, q_max=1.0)
    numeric = gui.compute_visibility(parameters, "numeric", samples=8000, grid=9, q_max=1.0)

    fitted_row = gui.compute_visibility(parameters, "fitted, l = 0,2,4",
                                        samples=8000, grid=9, q_max=1.0)
    first_order = gui.compute_visibility(parameters, "first-order", samples=8000,
                                         grid=9, q_max=1.0)

    index = np.unravel_index(np.argmax(np.abs(numeric.imag)), numeric.shape)
    print(f"\nIm at {index}: numeric {numeric.imag[index]:+.4f}, "
          f"affine {affine.imag[index]:+.4f}, "
          f"fitted {fitted_row.imag[index]:+.4f}, "
          f"first-order {first_order.imag[index]:+.4f}")
    for name, values in (("affine", affine), ("fitted", fitted_row),
                         ("first-order", first_order)):
        assert np.sign(values.imag[index]) == np.sign(numeric.imag[index]), name
        # the rows differ by the shape model and the truncation, both O(eps); a convention
        # error is O(1) and would show here long before the tolerance
        assert abs(values.imag[index] - numeric.imag[index]) < 0.02, name
