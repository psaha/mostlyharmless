"""The two-level fork: one brightness channel, no spot, no colour.

The committed map comes from ``make earth-levels``; these tests skip if it has not been
generated, so a fresh checkout can run the suite before anyone extracts anything.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import levels_gui

pytest.importorskip("mostlyharmless.symbolic.earth_levels")


def test_the_two_levels_have_the_right_ratio_and_mean():
    """Sea and land must sit at the committed ratio and average to one over the sphere."""
    from mostlyharmless.symbolic import earth_levels as emitted

    fraction = emitted.FRACTION
    sea, land = emitted.SEA, emitted.LAND
    mean = sea * (1.0 - fraction) + land * fraction
    print(f"\nfraction {fraction:.4f}, sea {sea:.3f}, land {land:.3f}, mean {mean:.4f}")
    assert land / sea == pytest.approx(emitted.RATIO, rel=1e-6)
    # 1.5 is the chosen ratio: visible, but a step small enough that the taper handles it
    assert emitted.RATIO == pytest.approx(1.5, rel=1e-6)
    assert mean == pytest.approx(1.0, abs=1e-4)


def test_the_levels_helper_solves_for_the_mean():
    """The closed form: s = 1/(1 + f(r-1)) makes the cosine-weighted mean exactly one."""
    from mostlyharmless.symbolic.earth import levels

    for fraction in (0.0, 0.292, 0.5, 1.0):
        sea, land = levels(fraction, ratio=1.5)
        mean = sea * (1.0 - fraction) + land * fraction
        assert mean == pytest.approx(1.0, rel=1e-12), f"fraction {fraction}"
        assert sea > 0.0


def test_the_fork_renders_one_channel_without_colour():
    """The panels are 2-D brightnesses --- no third axis, and no RGB compositing anywhere."""
    parameters = levels_gui.Parameters(l_max=4)
    near, far = levels_gui.compute_map(parameters, pixels=31)
    print(f"\npanels {near.shape}, peak {float(np.nanmax(near)):.3f}")
    assert near.shape == (31, 31)
    assert far.shape == (31, 31)
    # a two-level map is non-negative everywhere the disc is finite
    assert float(np.nanmin(near)) > 0.0
    # the faces must differ, or the renderer is showing one hemisphere twice
    assert np.nanmax(np.abs(near - far)) > 0.05


def test_the_fork_visibility_grid_is_finite():
    """The (u,v) view is the ordinary single-channel path, so it must be finite and vanish
    at the origin where Delta V is defined to be zero."""
    values = levels_gui.compute_visibility(levels_gui.Parameters(l_max=4), grid=11)
    assert values.shape == (11, 11)
    assert np.all(np.isfinite(values))
    midpoint = values.shape[0] // 2
    print(f"\nDelta V at the centre: {abs(values[midpoint, midpoint]):.2e}")
    # q = 0 is not exactly a grid point on an odd grid, so this is a small neighbourhood
    assert abs(values[midpoint, midpoint]) < 0.05


def test_every_view_and_resolution_computes():
    """Each radio setting produces finite panels of the documented shape.

    At the lowest degree on purpose: the assertions are about shapes and finiteness, which are
    degree-independent, while a render's cost is the mode count times the pixel count.
    """
    parameters = levels_gui.Parameters(l_max=4)
    for pixels, grid in levels_gui.RESOLUTIONS:
        near, far = levels_gui.compute_map(parameters, pixels=pixels)
        assert near.shape == far.shape == (pixels, pixels)
        values = levels_gui.compute_visibility(parameters, grid=grid)
        assert values.shape == (grid, grid)
        assert np.all(np.isfinite(values))
