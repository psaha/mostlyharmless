"""Rendering a coefficient set as a sky map: hemispheres, the shared basis, the display units.

``sky_map`` renders one hemisphere in the real basis, in the frame the coefficients belong
to; ``sky_map_batch`` does several sets against one shared Legendre table, which is what
makes the Earth's six panels affordable.  Used by ``symbolic.levels_gui``.

There is no matplotlib here.  This module is numerics only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.special import eval_legendre

from mostlyharmless.symbolic import affine, render
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.spherical_harmonics import normalisation
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk
from mostlyharmless.visylm.wigner import rotate_coefficients

__all__ = [
    "legendre_rows",
    "sky_map",
    "sky_map_batch",
]


# The convergence figure wants the largest jump that still stays non-negative.

# The a_00 = 1 convention fixes the sphere mean at 1/sqrt(4 pi), so the map figures
# plot I / <I> = I sqrt(4 pi).  This is a display choice; sky_map still returns model
# units, which the bound test needs.
MEAN_FACTOR = 1.0        # with Y_00 = 1 the panels are already in mean units


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def legendre_rows(order, l_max, x):
    """``P_l^order(x)`` for ``l = order .. l_max``, in one pass up the recurrence.

    With the Condon-Shortley phase, so these are the functions
    ``spherical_harmonics.associated_legendre`` returns, evaluated as a whole column of the
    table rather than one entry at a time.  That distinction is the whole point: calling
    ``lpmv`` per ``(l, m)`` re-walks the recurrence from ``m`` to ``l`` every time, so a
    panel costs ``O(L^3)`` rather than ``O(L^2)`` --- at ``L = 32`` that was fifteen seconds
    for a single panel.
    """
    x = np.asarray(x, dtype=float)
    root = np.sqrt(np.clip(1.0 - x * x, 0.0, None))

    value = np.ones_like(x)                       # P_0^0, then up to P_order^order
    for k in range(order):
        value = -(2 * k + 1) * root * value        # the phase lives here

    rows = [value]
    if l_max > order:
        value_next = (2 * order + 1) * x * value
        rows.append(value_next)
        for l in range(order + 2, l_max + 1):
            value_next, value = (
                ((2 * l - 1) * x * value_next - (l + order - 1) * value) / (l - order),
                value_next)
            rows.append(value_next)
    return rows


def sky_map_batch(coeff_sets, back=False, n=241):
    """The map on one hemisphere for several coefficient sets, sharing one basis.

    The Legendre functions and the azimuthal factors depend on the direction and the degree,
    not on the coefficients, so three colour channels can be rendered in one pass over the
    modes instead of three.  That is the second factor the timing table showed: ``all 6``
    was exactly three times the two faces.

    ``back`` renders the far hemisphere, mirrored so that it reads like a view from the
    other side.  Outside the disk the result is NaN.
    """
    coefficient_sets = tuple(coeff_sets)
    l_max = max(len(s.blocks) - 1 for s in coefficient_sets)

    axis = np.linspace(-1.0, 1.0, n)
    x_east, x_north = np.meshgrid(axis, axis)
    if back:
        x_east = -x_east
    rho_squared = x_east ** 2 + x_north ** 2
    z = (-1.0 if back else 1.0) * np.sqrt(np.clip(1.0 - rho_squared, 0.0, None))
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    # By D24 the body x axis is North, so the map's azimuth is measured from North.
    phi = np.arctan2(x_east, x_north)
    cos_phi = [np.cos(m * phi) for m in range(l_max + 1)]
    sin_phi = [np.sin(m * phi) for m in range(l_max + 1)]

    reals = [s.to_real() for s in coefficient_sets]
    images = [np.zeros_like(x_east) for _ in reals]

    for m in range(l_max + 1):
        rows = legendre_rows(m, l_max, np.cos(theta))
        for l in range(m, l_max + 1):
            p = rows[l - m]
            if m == 0:
                factor = normalisation(l, 0) * p
                for image, real in zip(images, reals):
                    image += real[l][l] * factor
            else:
                # R_{l,m,c} = sqrt(2) Re Y and R_{l,m,s} = sqrt(2) Im Y, and the real block
                # indexes them at l + m and l - m respectively
                scale = np.sqrt(2.0) * normalisation(l, m) * p
                cosine = scale * cos_phi[m]
                sine = scale * sin_phi[m]
                for image, real in zip(images, reals):
                    image += real[l][l + m] * cosine + real[l][l - m] * sine

    inside = rho_squared <= 1.0
    return tuple(np.where(inside, image, np.nan) for image in images)


def sky_map(coeffs, back=False, n=241):
    """The map on one hemisphere, as seen down the frame's z axis.

    One coefficient set.  The several-channel form is :func:`sky_map_batch`, which this
    calls with a single-element tuple so that the sum has one implementation.
    """
    return sky_map_batch((coeffs,), back=back, n=n)[0]


