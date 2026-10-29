"""The Earth extraction: the grid, the quadrature, the taper.

These know what the answer should be, which is the point --- the extraction was written
against a real image whose correct coefficients are unknown, and every bug it had was
found by a synthetic case rather than by looking at the picture.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import earth
from mostlyharmless.visylm.spherical_harmonics import real_sph_harm


def test_a_known_colour_plane_comes_back_with_its_own_coefficients():
    """Project a band-limited plane and check the quadrature recovers what was put in.

    A midpoint grid is not exact, so this is a tolerance test rather than an identity ---
    the honest statement about an image, and the reason the low degrees are the ones to
    trust.
    """
    n_lat, n_lon, l_max = 90, 180, 4
    theta = (np.arange(n_lat) + 0.5) * np.pi / n_lat
    phi = (np.arange(n_lon) + 0.5) * 2.0 * np.pi / n_lon - np.pi
    theta_grid = np.broadcast_to(theta[:, None], (n_lat, n_lon))
    phi_grid = np.broadcast_to(phi[None, :], (n_lat, n_lon))

    plane = 1.0 + 0.5 * real_sph_harm(2, 0, theta_grid, phi_grid)
    real = earth.channel(plane, l_max=l_max).to_real()

    print(f"\na00 {real[0][0]:+.4f}, a20 {real[2][2]:+.4f}")
    assert real[0][0] == pytest.approx(1.0, abs=5e-2)
    assert real[2][2] == pytest.approx(0.5, abs=5e-2)
    # a plane with no azimuthal structure must not produce the sine branch
    assert np.max(np.abs(real[2][:2])) < 1e-12


def test_the_taper_leaves_the_monopole_and_kills_the_top_degree():
    """A Lanczos σ is 1 at l = 0 and dies at the truncation, which is the point of it."""
    from mostlyharmless.visylm.coefficients import Coeffs

    blocks = tuple(np.ones(l + 1, dtype=complex) for l in range(5))
    tapered = earth.taper(Coeffs(blocks), l_max=4)

    # the taper multiplies each block by sinc(l / (l_max + 1)), a factor of its degree alone.
    # Monotonicity is a property of that factor, so it is checked on the blocks --- checking
    # it after conversion would fail, because to_real scales the m != 0 entries by sqrt(2)
    # and so puts the l = 1 block above the l = 0 one.
    values = [float(np.max(np.abs(block))) for block in tapered.blocks]
    expected = float(np.sinc(4 / 5.0))
    print(f"\nsigma per degree: {[round(v, 4) for v in values]}")
    assert values[0] == pytest.approx(1.0, abs=1e-12)
    assert values[-1] == pytest.approx(expected, rel=1e-12)
    assert all(a >= b for a, b in zip(values, values[1:]))

    # and the conversion's sqrt(2) lands where it belongs: on the m != 0 entries, with the
    # cosine slots from the real parts and the sine slots from the imaginary ones.  A block
    # of ones is real, so its sine slots are zero --- which is information, not a failure.
    real = tapered.to_real()
    assert real[2][2] == pytest.approx(values[2], rel=1e-12)                 # m = 0
    assert real[2][3] == pytest.approx(values[2] * np.sqrt(2.0), rel=1e-12)  # m = +1, cosine
    assert real[2][1] == pytest.approx(0.0, abs=1e-12)                       # m = -1, sine


def test_the_written_module_survives_a_round_trip(tmp_path):
    """Both parts must be written, because the imaginary part is the sine content.

    Dropping it makes every map even under phi -> -phi: a planet mirror-symmetric about its
    own meridian. That was a real bug, found by looking at the picture; this is the test.
    """
    from mostlyharmless.visylm.coefficients import Coeffs

    def blocks():
        """Blocks with a complex sine content but a real zonal term.

        ``Coeffs`` requires a_l0 to be real --- a real map has no sine term at m = 0 --- so a
        fixture that made every entry complex would be testing against an illegal object.
        """
        out = []
        for l in range(4):
            block = np.full(l + 1, 1.0 - 0.5j, dtype=complex)
            block[0] = 1.0
            out.append(block)
        return tuple(out)

    red, green, blue = Coeffs(blocks()), Coeffs(blocks()), Coeffs(blocks())
    target = earth.write(red, green, blue, l_max=3, target=tmp_path / "emitted.py")

    namespace = {}
    exec(compile(target.read_text(), str(target), "exec"), namespace)
    for name, original in (("RED", red), ("GREEN", green), ("BLUE", blue)):
        written = namespace[name]
        for l, (expected, got) in enumerate(zip(original.blocks, written)):
            assert np.allclose(np.asarray(got).real, expected.real, atol=1e-6)
            assert np.allclose(np.asarray(got).imag, expected.imag, atol=1e-6), \
                f"block {l} lost its imaginary part"
    print(f"\nround trip ok at l_max = {namespace['L_MAX']}")
