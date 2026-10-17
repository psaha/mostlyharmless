"""The figures render, and the bound behind the uv figure's excess over 1.

Nothing else imports the figure module, so without the rendering test a mistake in
it would only surface when someone runs `make data`.  That is how a tuple-unpacking
typo got through: the certificates were green, and the figure code had never been
executed by the test suite.

The rendering test needs matplotlib, which is optional, and skips without it.  The
bound test does not: cap_coefficients and sky_map use only numpy and scipy.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import figures
from mostlyharmless.visylm.visibility import visibility
from mostlyharmless.visylm.wigner import rotate_coefficients

# Small grids, so this stays fast.  The figures themselves use the defaults.
CASES = [
    ("map.pdf", "figure_map", {"l_max": 4, "n": 61}),
    ("uv.pdf", "figure_uv", {"l_max": 4, "n": 31}),
    ("convergence.pdf", "figure_convergence", {"orders": (2, 4), "n": 61}),
    ("jupiter.pdf", "figure_jupiter", {"l_max": 8, "n": 61}),
    ("triaxial.pdf", "figure_triaxial", {"n_samples": 3000, "n": 41}),
    ("stretched-sphere.pdf", "figure_stretched_sphere",
     {"n_samples": 3000, "n": 41}),
    ("delta-visibility.pdf", "figure_delta_visibility",
     {"n_samples": 1200, "grid": 7, "q_max": 1.0}),
]


@pytest.mark.parametrize("name,function,kwargs", CASES)
def test_the_figure_renders(tmp_path, name, function, kwargs):
    pytest.importorskip("matplotlib")
    target = tmp_path / name
    getattr(figures, function)(target, **kwargs)
    assert target.exists()
    assert target.stat().st_size > 0


def test_a_starspot_has_visibility_inside_the_unit_disc():
    """For a map with I >= 0, |V(u,v)| <= int I = V(0,0), so the normalised
    visibility stays inside the unit disc.

    The starspot map is built to be non-negative: a photosphere of 1 with a spot of
    depth 0.6, whose truncation ringing undershoots by only about 9 per cent of the
    depth and so stays above zero.  The first assertion is what makes the last one
    legitimate --- the bound needs I >= 0 pointwise --- and this replaces the earlier
    absolute-mass test, which needed a map that went negative.  The small slack is
    for the grid not resolving every point of the truncated map.
    """
    n = 601
    weight = (2.0 / (n - 1)) ** 2
    coeffs = rotate_coefficients(figures.cap_coefficients(16),
                                 0.0, figures.BETA, figures.PA)
    image = figures.sky_map(coeffs, n=n)
    flux = weight * float(np.nansum(image))
    analytic = float(np.real(visibility(coeffs, 0.0, 0.0)))

    axis_q = np.linspace(-6.0, 6.0, 121)
    u, v = np.meshgrid(axis_q, axis_q)
    values = visibility(coeffs, u, v) / analytic

    assert float(np.nanmin(image)) > 0.0
    assert flux == pytest.approx(analytic, rel=0.05)
    assert float(np.abs(values).max()) <= 1.0 + 1e-6
