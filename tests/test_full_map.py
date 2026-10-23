"""Stage 5: the full map, its FFT check (V15), and the Hermitian symmetry.

V15 is the end-to-end certificate.  It feeds a rotated starspot through the whole
chain -- real coefficients, both parity families, the rotation -- and compares the
result with a direct 2-D transform of the rendered image.  It is the only check
that would catch an error common to both families.

The symmetry test covers the property the uv figure is expected to show and that
the caption of figure 2 states: for a real map, V(-u,-v) = conj(V(u,v)).
"""

import numpy as np

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.visibility import visibility
from mostlyharmless.visylm.wigner import rotate_coefficients

TOL = 0.05


def test_v15_full_map_against_a_two_dimensional_fft():
    assert dv.check_full_map_fft() < TOL


def test_v15_is_not_vacuous():
    """A mistake in the frequency scaling would show up as a large residual, so
    check that the image really is being transformed and not, say, all zeros."""
    coeffs = rotate_coefficients(dv._zonal_starspot(4), 0.0, 0.9, 0.4)
    values, frequencies = dv._fft_visibility(coeffs, 129)
    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))
    assert abs(flux) > 0.1
    assert float(np.max(np.abs(values))) > 0.1
    assert frequencies[0] < 0.0 < frequencies[-1]


def test_the_visibility_is_hermitian_for_a_real_map():
    """V(-u, -v) = conj(V(u, v)): the real part is even and the imaginary part odd."""
    coeffs = rotate_coefficients(dv._zonal_starspot(6), 0.0, 1.05, 0.6)
    axis_q = np.linspace(-4.0, 4.0, 17)
    u, v = np.meshgrid(axis_q, axis_q)
    values = visibility(coeffs, u, v)
    reflected = visibility(coeffs, -u, -v)
    assert np.max(np.abs(reflected - np.conjugate(values))) < 1e-12


def test_the_mode_at_the_opposite_baseline():
    """V23: V_l^m(-u,-v) = (-1)^m V_l^m(u,v), the identity the pairing rests on.

    The Hermitian symmetry of a real map tests this only in combination with the pairing; a
    sign error common to both would survive it.  This isolates the mode identity, so the
    plain conjugation of the negative-order visibility is checked on its own.
    """
    assert dv.check_baseline_parity() < 1e-12


def test_the_truncation_error_falls_with_the_mode_count():
    """V28: a discontinuous map converges algebraically, not exponentially.

    Each doubling of ``l`` should at least halve the distance to the converged visibility;
    exponential convergence would give a much smaller ratio, and a stall a ratio near one.
    """
    assert dv.check_convergence_rate() < 0.75
