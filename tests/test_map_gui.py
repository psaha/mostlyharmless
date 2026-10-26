"""The map GUI: the spot's coefficients, its two hemispheres, and its visibility.

The module has no deformation, no level solve and no fit, so these tests are about the
conventions the manual's figures fix: the map normalised to its mean, a non-negative
brightness, and Delta V vanishing for a uniform disk.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import map_gui
from mostlyharmless.visylm.coefficients import Coeffs


def test_every_l_max_is_band_limited():
    """The coefficients carry nothing above the chosen degree."""
    for l_max in map_gui.L_MAX_LABELS:
        parameters = map_gui.MapParameters(l_max=l_max)
        coeffs = map_gui.coefficients(parameters)
        assert len(coeffs.blocks) == l_max + 1
        assert np.max(np.abs(coeffs.blocks[l_max])) >= 0.0


def test_a_zero_depth_spot_is_a_uniform_disk():
    """No spot, so no departure from the uniform disk: Delta V vanishes identically."""
    parameters = map_gui.MapParameters(depth=0.0, l_max=16)
    values = map_gui.compute_visibility(parameters, grid=11, q_max=2.0)
    print(f"\nuniform map: max |Delta V| = {np.max(np.abs(values)):.2e}")
    assert np.max(np.abs(values)) < 1e-12


def test_the_spot_is_scaled_to_unit_monopole():
    """``a_00 = 1`` fixes the sphere mean, not the flux.

    The flux is the *cosine-weighted* integral, which the visibility uses, so a spotted map's
    ``V(0,0)`` is less than the uniform disk's: darkening the pole removes area that the
    viewing cosine favours.  The coefficient is the testable quantity; the flux is not.
    """
    from mostlyharmless.visylm.visibility import visibility

    parameters = map_gui.MapParameters(depth=0.6, beta=1.05, l_max=16)
    coeffs = map_gui.coefficients(parameters)
    assert abs(coeffs.blocks[0][0] - 1.0) < 1e-12

    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))
    uniform = np.sqrt(np.pi) / 2.0                  # the anchor of the manual's section 4
    print(f"\nflux: spotted {flux:.4f}, uniform {uniform:.4f}")
    assert 0.0 < flux < uniform

    values = map_gui.compute_visibility(parameters, grid=11, q_max=2.0)
    assert np.all(np.isfinite(values))


def test_both_hemispheres_render_and_differ():
    """Near and far are different pictures once the spot is off the spin axis."""
    parameters = map_gui.MapParameters(depth=0.6, beta=1.0, l_max=16)
    near, far = map_gui.compute_map(parameters, pixels=61)
    assert near.shape == far.shape == (61, 61)
    # the rim is NaN by design, as in the renderer's own images
    assert np.isfinite(near).sum() > 0.5 * near.size
    assert np.nanmax(np.abs(near - far)) > 0.1


def test_the_mean_removed_map_is_unnormalised_by_the_square():
    """The map is in units of the mean, so its bright floor sits near 1 - depth."""
    parameters = map_gui.MapParameters(depth=0.0, l_max=16)
    near, _ = map_gui.compute_map(parameters, pixels=61)
    inside = near[np.isfinite(near)]
    print(f"\nuniform map level: mean {inside.mean():.4f}, "
          f"spread {np.ptp(inside):.2e}")
    assert abs(inside.mean() - 1.0) < 1e-2
    assert np.ptp(inside) < 1e-2


def test_every_view_and_resolution_computes():
    """Each radio setting produces finite panels of the documented shape.

    Run at the spot's lowest degree on purpose.  The assertions are about shapes and
    finiteness, which are degree-independent, while the cost of a render is the mode count
    times the pixel count --- at the default degree the 361-pixel setting alone took most of
    the suite's runtime.
    """
    parameters = map_gui.MapParameters(l_max=4)
    for pixels, grid in map_gui.RESOLUTIONS:
        near, far = map_gui.compute_map(parameters, pixels=pixels)
        assert near.shape == (pixels, pixels)
        values = map_gui.compute_visibility(parameters, grid=grid)
        assert values.shape == (grid, grid)
        assert np.all(np.isfinite(values))


def test_the_rotation_phase_is_inert_at_the_pole_and_moves_a_tilted_spot():
    """The degeneracy that made alpha look broken, pinned as a property.

    A rotation about the spin axis is the identity on an axisymmetric map, so with the spot
    at the pole alpha does nothing; off the pole it is the spot's longitude.
    """
    polar = map_gui.MapParameters(colatitude=0.0, alpha=0.0, beta=0.6, l_max=16)
    turned = map_gui.MapParameters(colatitude=0.0, alpha=2.0, beta=0.6, l_max=16)
    near_polar, _ = map_gui.compute_map(polar, pixels=61)
    near_turned, _ = map_gui.compute_map(turned, pixels=61)
    print(f"\npolar:   max |difference| = "
          f"{np.nanmax(np.abs(near_polar - near_turned)):.2e}")

    off_axis = map_gui.MapParameters(colatitude=1.0, alpha=0.0, beta=0.6, l_max=16)
    off_turned = map_gui.MapParameters(colatitude=1.0, alpha=2.0, beta=0.6, l_max=16)
    near_off, _ = map_gui.compute_map(off_axis, pixels=61)
    near_off_turned, _ = map_gui.compute_map(off_turned, pixels=61)
    print(f"tilted:  max |difference| = "
          f"{np.nanmax(np.abs(near_off - near_off_turned)):.2e}")

    assert np.nanmax(np.abs(near_polar - near_turned)) < 1e-9
    assert np.nanmax(np.abs(near_off - near_off_turned)) > 0.1


def test_the_position_angle_rotates_the_picture():
    """PA is no longer a gauge once the spot has an orientation on the sky."""
    upright = map_gui.MapParameters(colatitude=1.0, beta=0.6, pa=0.0, l_max=16)
    turned = map_gui.MapParameters(colatitude=1.0, beta=0.6, pa=0.8, l_max=16)
    near_upright, _ = map_gui.compute_map(upright, pixels=61)
    near_turned, _ = map_gui.compute_map(turned, pixels=61)
    print(f"\nPA: max |difference| = "
          f"{np.nanmax(np.abs(near_upright - near_turned)):.2e}")
    assert np.nanmax(np.abs(near_upright - near_turned)) > 0.05


def test_the_earth_mode_is_three_colour_channels():
    """The Earth's panels are RGB, so the two faces share one colour balance.

    The emitted coefficients are committed, so this runs without the image and without
    pyshtools --- which is the whole point of emitting them.
    """
    # l = 4 and 31 pixels: the properties below are about the RGB structure, not the degree,
    # and six renders at the emitted degree are what made this module slow
    parameters = map_gui.MapParameters(alpha=0.0, beta=np.radians(60.0), pa=0.0, l_max=4)
    near, far = map_gui.compute_map(parameters, pixels=31, model="Earth")
    # the rim is NaN by design, as in the renderer's own images, so every statistic here
    # has to be the nan-aware one
    print(f"\nEarth panels: {near.shape}, peak {float(np.nanmax(near)):.3f}")
    assert near.shape == (31, 31, 3)
    assert far.shape == (31, 31, 3)
    # normalised to the brightest channel, so the colour is a ratio and not an absolute
    assert float(np.nanmax(near)) == pytest.approx(1.0, rel=1e-6)
    assert float(np.nanmin(near)) >= 0.0
    assert np.isfinite(near).sum() > 0.5 * near.size
    # the two faces must differ, or the renderer is showing one hemisphere twice
    assert np.nanmax(np.abs(near - far)) > 0.1


def test_the_earth_is_tapered_so_the_continents_are_visible():
    """Un-tapered, the coastline ringing is larger than the continents themselves."""
    from mostlyharmless.symbolic import earth, earth_coefficients as emitted

    raw = tuple(Coeffs(tuple(np.asarray(block, dtype=complex) for block in getattr(emitted, name)))
                for name in ("RED", "GREEN", "BLUE"))
    top_raw = np.max(np.abs(raw[0].blocks[emitted.L_MAX]))
    top_tapered = np.max(np.abs(map_gui.earth_coefficients()[0].blocks[emitted.L_MAX]))
    print(f"\ntop degree, raw {float(top_raw):.3e} vs tapered {float(top_tapered):.3e}")
    assert top_tapered < 0.3 * top_raw
    assert float(np.max(np.abs(map_gui.earth_coefficients()[0].blocks[0]))) > 0.0
