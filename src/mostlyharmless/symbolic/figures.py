"""Figures, for looking at the results and for the manual.

Run from the repository root:

    python -m mostlyharmless.symbolic.figures

matplotlib is an optional dependency (the ``figures`` extra) and is imported inside
the functions, so the rest of the package does not need it.  Output is vector PDF
into docs/generated/.  The sky maps are built from the real basis with the same
pattern the tests use; there is no public renderer, by decision D16.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.special import eval_legendre

from mostlyharmless.symbolic import affine, render
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.spherical_harmonics import real_sph_harm
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk
from mostlyharmless.visylm.wigner import rotate_coefficients

__all__ = [
    "cap_coefficients",
    "sky_map",
    "figure_map",
    "figure_uv",
    "figure_convergence",
    "write",
]

# A cap this size has a peak brightness about four and a half times its mean, so
# the edge and its ringing are both visible.
HALF_ANGLE = 0.9
BETA = 1.05
PA = 0.6

# A spot at 40 per cent of the photosphere: dark enough to read as a starspot, and
# far from the depth at which the truncation ringing would go negative.
SPOT_DEPTH = 0.6
# The convergence figure wants the largest jump that still stays non-negative.
CONVERGENCE_DEPTH = 0.9

# The a_00 = 1 convention fixes the sphere mean at 1/sqrt(4 pi), so the map figures
# plot I / <I> = I sqrt(4 pi).  This is a display choice; sky_map still returns model
# units, which the bound test needs.
MEAN_FACTOR = float(np.sqrt(4.0 * np.pi))


def _pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def cap_coefficients(l_max=16, half_angle=HALF_ANGLE, depth=SPOT_DEPTH):
    """A polar spot on a uniformly bright surface, as zonal coefficients.

    The surface has brightness 1 outside the cap and 1 - depth inside it, so depth = 0
    is no spot and depth = 1 is black.  The Legendre moments of that map are
    c_l = 2 delta_l0 - depth m_l, with m_l the moments of the cap indicator, and the
    harmonic coefficients follow from a_l0 = c_l sqrt(pi (2l+1)).  Everything is then
    scaled so that a_00 = 1, which fixes the sphere mean at 1/sqrt(4 pi) by D16.

    Keeping depth below about 0.92 keeps the map non-negative: the truncation ringing
    undershoots by roughly 9 per cent of the jump.  That is what makes |V| <= V(0,0)
    hold for this map, which is the bound a starspot should show.

    The blocks are built directly rather than by indexing the flat real vector: there
    each block has length 2l+1 with (l, 0) in the middle, so that index is l*l + l.
    The first attempt used l*l, which put the zonal coefficients in the sin(l*phi)
    slots and drew a banana instead of a cap.
    """
    cos_edge = np.cos(half_angle)
    nodes, weights = np.polynomial.legendre.leggauss(200)
    mu = 0.5 * (nodes + 1.0) * (1.0 - cos_edge) + cos_edge
    weight = 0.5 * (1.0 - cos_edge) * weights
    moments = [float(np.sum(weight * eval_legendre(l, mu))) for l in range(l_max + 1)]
    zeroth = 2.0 - depth * moments[0]
    blocks = []
    for l in range(l_max + 1):
        spherical = (2.0 if l == 0 else 0.0) - depth * moments[l]
        block = np.zeros(l + 1, dtype=complex)
        block[0] = (spherical / zeroth) * np.sqrt(2 * l + 1)
        blocks.append(block)
    return Coeffs(tuple(blocks))


def sky_map(coeffs, back=False, n=241):
    """The map on one hemisphere, as seen down the frame's z axis.

    Each point of the projected disk is a direction on the unit sphere, and the map
    is evaluated there in the real basis, in the frame the coefficients belong to.
    ``back`` renders the far hemisphere, mirrored so that it reads like a view from
    the other side.  Outside the disk it is NaN.
    """
    axis = np.linspace(-1.0, 1.0, n)
    x_east, x_north = np.meshgrid(axis, axis)
    if back:
        x_east = -x_east
    rho_squared = x_east ** 2 + x_north ** 2
    z = (-1.0 if back else 1.0) * np.sqrt(np.clip(1.0 - rho_squared, 0.0, None))
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    # By D24 the body x axis is North, so the map's azimuth is measured from North.
    phi = np.arctan2(x_east, x_north)
    image = np.zeros_like(x_east)
    real = coeffs.to_real()
    for l, row in enumerate(real):
        image = image + row[l] * real_sph_harm(l, 0, theta, phi)
        for m in range(1, l + 1):
            image = image + (row[l + m] * real_sph_harm(l, m, theta, phi)
                             + row[l - m] * real_sph_harm(l, -m, theta, phi))
    return np.where(rho_squared <= 1.0, image, np.nan)


def figure_map(target, l_max=16, n=241):
    """The spot on both hemispheres: about the spin axis, and in the observer's frame.

    Brightness is in units of the mean, I/<I> = I sqrt(4 pi).  The unspotted surface
    is the floor, so the map is non-negative and the ringing at the rim stays above
    zero; the spot reads 1 - depth of the photosphere.  A display choice; the
    library's units and the a_00 = 1 convention are untouched.
    """
    plt = _pyplot()
    body = cap_coefficients(l_max)
    observer = rotate_coefficients(body, 0.0, BETA, PA)
    panels = ((MEAN_FACTOR * sky_map(body, n=n), "spin axis, near side"),
              (MEAN_FACTOR * sky_map(body, back=True, n=n), "spin axis, far side"),
              (MEAN_FACTOR * sky_map(observer, n=n), "observer frame, near side"),
              (MEAN_FACTOR * sky_map(observer, back=True, n=n),
               "observer frame, far side"))
    limits = (float(np.nanmin(panels[0][0])), float(np.nanmax(panels[0][0])))
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 6.6), constrained_layout=True)
    for axis, (data, title) in zip(axes.ravel(), panels):
        # East increases to the left (D24): the extent draws column 0 on the left, so the
        # array is flipped once here rather than built mirrored
        drawn = axis.imshow(data[:, ::-1], origin="lower", extent=(1, -1, -1, 1),
                            cmap="magma", vmin=limits[0], vmax=limits[1])
        axis.set_title(title, fontsize=9)
        axis.set_xlabel(r"$y$", fontsize=9)
    for axis in (axes[0, 0], axes[1, 0]):
        axis.set_ylabel(r"$x$", fontsize=9)
    fig.colorbar(drawn, ax=axes, shrink=0.8, label=r"$I/\langle I\rangle$")
    fig.savefig(target)
    plt.close(fig)


def figure_uv(target, l_max=16, q_max=6.0, n=301):
    """The normalised visibility of the rotated cap, real and imaginary parts.

    V is divided by its value at the origin, which is the total flux, so the real
    part is 1 at the origin by construction and the imaginary part is 0 there: the
    visibility of a real map is real at zero baseline (D15).

    The map is non-negative --- a spot on a bright surface --- so |V| <= V(0,0) and
    both panels stay within the unit range.  The shared symmetric scale gives no other
    way to read off where 1 sits, so each panel is annotated with its own maximum.
    """
    plt = _pyplot()
    coeffs = rotate_coefficients(cap_coefficients(l_max), 0.0, BETA, PA)
    axis_q = np.linspace(-q_max, q_max, n)
    u, v = np.meshgrid(axis_q, axis_q)
    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))
    values = visibility(coeffs, u, v) / flux
    # D24: u is North, so it runs vertically (the array's first index) and v, East,
    # increases leftward.
    values = values.T
    limit = float(max(np.abs(values.real).max(), np.abs(values.imag).max()))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), constrained_layout=True)
    panels = ((values.real, r"$\mathrm{Re}\;V/V(0,0)$"),
              (values.imag, r"$\mathrm{Im}\;V/V(0,0)$"))
    for axis, (data, title) in zip(axes, panels):
        drawn = axis.imshow(data[:, ::-1], origin="lower",
                            extent=(q_max, -q_max, -q_max, q_max),
                            cmap="RdBu_r", vmin=-limit, vmax=limit)
        axis.set_title(title, fontsize=10)
        axis.set_xlabel(r"$v$ (wavelengths, East)", fontsize=9)
        axis.set_ylabel(r"$u$ (wavelengths, North)", fontsize=9)
        axis.text(0.03, 0.97, rf"$\max|\cdot| = {np.abs(data).max():.3f}$",
                  transform=axis.transAxes, verticalalignment="top", fontsize=8,
                  bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"})
        fig.colorbar(drawn, ax=axis, shrink=0.85)
    fig.savefig(target)
    plt.close(fig)


def figure_convergence(target, orders=(4, 8, 16), n=241):
    """The same sharp-edged map at growing l_max, on one colour scale.

    A map that is discontinuous on the sphere converges slowly, so the ripples near
    the edge are the truncation error, not a mistake: this is the ringing that V15
    has to state a scope for.  Brightness is in units of the mean, as in figure_map.
    """
    plt = _pyplot()
    images = [MEAN_FACTOR
              * sky_map(rotate_coefficients(
                  cap_coefficients(l_max, depth=CONVERGENCE_DEPTH), 0.0, BETA, PA),
                  n=n)
              for l_max in orders]
    last = images[-1]
    limits = (float(np.nanmin(last)), float(np.nanmax(last)))
    fig, axes = plt.subplots(1, len(orders), figsize=(9.6, 3.4),
                             constrained_layout=True)
    drawn = None
    for axis, data, l_max in zip(axes, images, orders):
        # East increases to the left (D24): the extent draws column 0 on the left, so the
        # array is flipped once here rather than built mirrored
        drawn = axis.imshow(data[:, ::-1], origin="lower", extent=(1, -1, -1, 1),
                            cmap="magma", vmin=limits[0], vmax=limits[1])
        axis.set_title(rf"$l_{{max}}={l_max}$", fontsize=9)
        axis.set_xticks([])
        axis.set_yticks([])
    fig.colorbar(drawn, ax=axes, shrink=0.85, label=r"$I/\langle I\rangle$")
    fig.savefig(target)
    plt.close(fig)


def write(directory=None):
    """Write the figures, and return the paths written."""
    if directory is None:
        directory = Path(__file__).resolve().parents[3] / "docs" / "generated"
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    # Only the figures the manual embeds are generated for now; Stage 6's figures are
    # regenerated on request while they are still being worked on.
    for name, function in (("map.pdf", figure_map),
                           ("uv.pdf", figure_uv),
                           ("convergence.pdf", figure_convergence)):
        path = directory / name
        function(path)
        written.append(path)
    return written


if __name__ == "__main__":
    for written in write():
        print(written)
