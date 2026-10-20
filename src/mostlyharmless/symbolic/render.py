"""Render a deformed star, and transform it directly.

Stage 6, item 2.  The surface is the level set ``Phi = -1`` of a
:class:`~mostlyharmless.visylm.deformation.Deformation`, and the brightness is set by
the bolometric von Zeipel law ``I ∝ g_eff = |grad Phi|`` (D21; the exponent is a
parameter, so a weaker convective law is a change of number, not of code).

Rather than following rays inward from the observer, the surface is sampled by
direction from the centre and forward-mapped onto the sky plane.  For a convex body
that map is one-to-one, and because the surface is an equipotential the outward
normal is just the gradient, so visibility is the sign of ``N . e_obs`` with no
z-buffer and no nested root find.  Each sample carries the *sky-plane* area it
projects to, which is what our visibility integrates over.

The sample directions are spread by the golden-angle spiral, so each owns the same
solid angle.
"""

from __future__ import annotations

import numpy as np

from mostlyharmless.visylm.deformation import Deformation

__all__ = [
    "sphere_directions",
    "sky_basis",
    "surface_samples",
    "sphere_samples",
    "sky_image",
    "direct_visibility",
]


def sphere_directions(n_samples):
    """n directions spread evenly over the sphere, by the golden-angle spiral."""
    index = np.arange(n_samples) + 0.5
    mu = 1.0 - 2.0 * index / n_samples
    phi = np.pi * (1.0 + 5.0 ** 0.5) * index
    sine = np.sqrt(1.0 - mu ** 2)
    return np.stack([sine * np.cos(phi), sine * np.sin(phi), mu], axis=-1)


def sky_basis(line_of_sight):
    """Two orthonormal vectors spanning the sky plane, east and north."""
    eye = np.asarray(line_of_sight, dtype=float)
    eye = eye / np.linalg.norm(eye)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(eye @ helper)) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    east = np.cross(helper, eye)
    east = east / np.linalg.norm(east)
    north = np.cross(eye, east)
    return east, north


def surface_samples(field, line_of_sight, n_samples=20000, power=1.0):
    """Sample the level surface, keeping the visible points.

    Returns ``(points, brightness, weights)``.  The brightness is ``g^power`` in the
    potential's own units, where the undeformed ``g_eff`` is 1, so views can share a
    colour scale.  The weight of a sample is the area it projects onto the sky plane,
    ``(N . e_obs) r^2 / (n . N)`` times the solid angle it owns, since our visibility
    integrates over the plane rather than over the surface.

    The level radius comes from :meth:`Deformation.level_radii`, a vectorised bisection.
    """
    eye = np.asarray(line_of_sight, dtype=float)
    eye = eye / np.linalg.norm(eye)
    directions = sphere_directions(n_samples)

    radii = field.level_radii(directions)
    points = radii[:, None] * directions

    gradients = field.gradient(points)
    magnitudes = np.linalg.norm(gradients, axis=-1)
    normals = gradients / magnitudes[:, None]
    normal_radial = np.sum(normals * directions, axis=-1)
    normal_eye = np.sum(normals * eye, axis=-1)

    visible = normal_eye > 0.0
    if not np.any(visible):
        raise ValueError("nothing is visible along this line of sight")
    solid = 4.0 * np.pi / n_samples
    weights = (normal_eye * radii ** 2 / np.maximum(normal_radial, 1e-12) * solid)
    weights = np.where(visible, weights, 0.0)
    # g^power in the potential's own units, where the undeformed g_eff is 1.  Not
    # renormalised per view, so two orientations can share one colour scale.
    brightness = np.where(visible, magnitudes ** power, 0.0)
    return points[visible], brightness[visible], weights[visible]


def sphere_samples(field, line_of_sight, n_samples=20000, power=1.0):
    """The star's brightness pattern carried onto the *unit sphere*.

    The affine model's counterpart to :func:`surface_samples`.  The *pattern* is the
    star's own -- ``|grad Phi|`` on the level set, so the contribution of the surface
    having *moved* is included -- while the *geometry* is the sphere's, so visibility,
    projection and the areas are radial.

    Evaluating the gravity at ``r = 1`` instead would be a different object: it keeps
    only the perturbation's own derivative and drops the ``-2 delta r / r`` term, which
    leaves a pattern about half as deep as the star's, since both terms are first order
    and equal for quadrupoles.
    """
    eye = np.asarray(line_of_sight, dtype=float)
    eye = eye / np.linalg.norm(eye)
    directions = sphere_directions(n_samples)

    radii = field.level_radii(directions)
    magnitudes = np.linalg.norm(field.gradient(radii[:, None] * directions), axis=-1)
    facing = directions @ eye

    visible = facing > 0.0
    if not np.any(visible):
        raise ValueError("nothing is visible along this line of sight")
    solid = 4.0 * np.pi / n_samples
    weights = np.where(visible, facing * solid, 0.0)
    brightness = np.where(visible, magnitudes ** power, 0.0)
    return (directions[visible], brightness[visible], weights[visible])


def sphere_brightness(field, power=1.0, n_samples=20000):
    """``|grad Phi|^power`` on the unit sphere, over the *whole* sphere.

    Returns ``(directions, brightness)``.  The directions are the golden-angle spiral,
    so each owns the same solid angle ``4 pi / n_samples`` and a projection onto the
    harmonics is a plain weighted average.  This is the pattern the affine model puts
    on a sphere: taken from the level set, so the contribution of the surface having
    moved is included, but evaluated as a function of direction.
    """
    directions = sphere_directions(n_samples)
    radii = field.level_radii(directions)
    magnitudes = np.linalg.norm(field.gradient(radii[:, None] * directions), axis=-1)
    return directions, magnitudes ** power


def sky_image(points, brightness, weights, east, north, half_width=1.6, n=241):
    """Bin the samples onto the sky plane, as a flux-weighted mean brightness.

    The array is indexed [north, east], the image convention, so it can be handed to
    imshow directly.  Each sample is spread bilinearly over the four pixels around it
    rather than dropped into one: a nearest-pixel deposit speckles whenever the
    samples do not comfortably outnumber the pixels of the disk.
    """
    scale = n / (2.0 * half_width)
    x = (points @ east + half_width) * scale - 0.5
    y = (points @ north + half_width) * scale - 0.5
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    dx, dy = x - x0, y - y0

    inside = (x0 >= 0) & (x0 < n - 1) & (y0 >= 0) & (y0 < n - 1)
    column, row = x0[inside], y0[inside]
    flux = np.zeros((n, n))
    area = np.zeros((n, n))
    for offset_x, offset_y, share in ((0, 0, (1.0 - dx) * (1.0 - dy)),
                                      (1, 0, dx * (1.0 - dy)),
                                      (0, 1, (1.0 - dx) * dy),
                                      (1, 1, dx * dy)):
        share = share[inside]
        np.add.at(flux, (row + offset_y, column + offset_x),
                  (brightness * weights)[inside] * share)
        np.add.at(area, (row + offset_y, column + offset_x),
                  weights[inside] * share)
    return np.where(area > 0.0, flux / np.maximum(area, 1e-30), np.nan)


def direct_visibility(points, brightness, weights, east, north, u, v):
    """V(u,v) by direct summation over the samples, in our planar measure.

    A quadrature, not a periodic transform, so its error is ordinary discretisation
    error that shrinks with the sampling rather than the Gibbs ringing an FFT would
    impose at the limb.

    ``u`` and ``v`` are scalars: the phase runs over the sample array, so a grid of
    frequencies is a loop over this, not one call.  By D24 ``u`` is the North component
    of the baseline and ``v`` the East one.
    """
    x_east = points @ east
    x_north = points @ north
    phase = np.exp(-2j * np.pi * (np.asarray(u, dtype=float) * x_north
                                  + np.asarray(v, dtype=float) * x_east))
    return np.sum(brightness * weights * phase)
