"""The renderer and the cap: the two pieces the GUIs share.

These survive the retirement of the RGB map GUI, because they test the model layer rather
than a window: the coefficients a cap carries, the two hemispheres, the normalisation to the
mean, and the fast renderer against the slow per-mode sum it replaced.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import maps
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.wigner import rotate_coefficients


def test_a_cap_carries_nothing_above_its_degree():
    """The coefficients are band-limited by construction, at every degree."""
    for l_max in (4, 8, 16, 32):
        coeffs = maps.cap_coefficients(l_max=l_max)
        assert len(coeffs.blocks) == l_max + 1
        assert np.max(np.abs(coeffs.blocks[l_max])) >= 0.0


def test_both_hemispheres_render_and_differ():
    """Near and far are different pictures once the cap is off the spin axis."""
    coeffs = maps.cap_coefficients(l_max=16)
    tilted = rotate_coefficients(coeffs, 0.0, 1.0, 0.0)
    near = maps.sky_map(tilted, n=61)
    far = maps.sky_map(tilted, back=True, n=61)
    assert near.shape == far.shape == (61, 61)
    # the rim is NaN by design, as in the renderer's own images
    assert np.isfinite(near).sum() > 0.5 * near.size
    assert np.nanmax(np.abs(near - far)) > 0.1


def test_a_uniform_map_comes_out_at_one():
    """The map is in units of the sphere mean, so a uniform map is flat at 1.

    The coefficient is sqrt(4 pi), not 1: a_00 multiplies Y_00, which is 1/sqrt(4 pi), so
    a_00 = 1 would give a flat map at 1/sqrt(4 pi) = 0.282 instead of at one.  That is the
    same convention trap as a_00 not being the flux, and it caught this test.
    """
    def disc_mean(monopole):
        blocks = tuple(np.zeros(l + 1, dtype=complex) for l in range(17))
        blocks[0][0] = monopole
        near = maps.sky_map(Coeffs(blocks), n=61)
        inside = near[np.isfinite(near)]
        assert np.ptp(inside) < 1e-12, "a uniform map must be uniform"
        return float(inside.mean())

    # D16: the normalisation is a_00 = 1, which puts the sphere mean at 1/sqrt(4 pi)
    normalised = disc_mean(1.0)
    print(f"\na_00 = 1 gives {normalised:.4f} = 1/sqrt(4 pi) = "
          f"{1 / np.sqrt(4 * np.pi):.4f}")
    assert normalised == pytest.approx(1.0 / np.sqrt(4.0 * np.pi), rel=1e-9)

    # a_00 = sqrt(4 pi) is the same map rescaled to mean one, which is what the viewer
    # multiplies by before display
    assert disc_mean(np.sqrt(4.0 * np.pi)) == pytest.approx(1.0, rel=1e-9)


def test_the_fast_renderer_matches_the_slow_reference():
    """The rewritten sky_map must equal the per-mode sum it replaced.

    The fast path builds the Legendre functions by a single recurrence up each column of the
    table, where the old one called ``lpmv`` per (l, m).  A slip in the phase, the
    normalisation or the sine sign would give a wrong map rather than a slow one, and the
    pictures would still look like a planet, so this compares against the direct sum.
    """
    from mostlyharmless.visylm.spherical_harmonics import real_sph_harm

    l_max, n = 4, 31
    rng = np.random.default_rng(7)
    # the real basis indexes a block by m + l, so each has 2l + 1 entries
    plan = [rng.standard_normal(2 * l + 1) for l in range(l_max + 1)]
    plan[0][0] = 2.0                                           # a definite monopole
    coefficients = Coeffs.from_real(plan)

    axis = np.linspace(-1.0, 1.0, n)
    x_east, x_north = np.meshgrid(axis, axis)
    rho_squared = x_east ** 2 + x_north ** 2
    z = np.sqrt(np.clip(1.0 - rho_squared, 0.0, None))
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    phi = np.arctan2(x_east, x_north)

    expected = np.zeros((n, n))
    real = coefficients.to_real()
    for l, row in enumerate(real):
        expected += row[l] * real_sph_harm(l, 0, theta, phi)
        for m in range(1, l + 1):
            expected += (row[l + m] * real_sph_harm(l, m, theta, phi)
                         + row[l - m] * real_sph_harm(l, -m, theta, phi))

    fast = maps.sky_map(coefficients, n=n)
    inside = rho_squared <= 1.0
    difference = float(np.max(np.abs(fast[inside] - expected[inside])))
    print(f"\nfast vs per-mode sum: {difference:.2e}")
    assert difference < 1e-10
    # and outside the disk it is still NaN, which the comparison above would hide
    assert np.all(np.isnan(fast[~inside]))
