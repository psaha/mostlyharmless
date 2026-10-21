"""Rung 1: the affine model of a deformed star.

Stage 6, item 4, and section 8.  The image of a linearly stretched body is the
stretched image of the original, so stretching the sky coordinates by ``Lambda^-1``
stretches the frequency coordinates by ``Lambda``:

    V_affine(u) = V_sphere(Lambda^T u) .

No flux factor appears: the ``1/|det Lambda|`` that the change of variables puts on
the image cancels against the ``|det Lambda|`` from the Jacobian.  ``V_sphere`` is the
visibility of a *sphere* carrying the same brightness pattern, which is what makes
this a model of the shape alone -- the gravity darkening is carried through, not
discarded along with the deformation.

``Lambda`` is the square root of the projected shape.  For an ellipsoid with
semi-axes ``(a, b, c)`` along the body axes, the sky-plane shape is
``A = P diag(a, b, c)`` with ``P`` the two rows ``(east, north)``, and
``Lambda = sqrt(A A^T)``.  The silhouette of an ellipsoid is an ellipse, so the
projection itself is exact; the approximation is that the deformed body *is* an
ellipsoid, and that the stretch carries the brightness pattern rigidly with it.
"""

from __future__ import annotations

import numpy as np

from mostlyharmless.symbolic import render
from mostlyharmless.symbolic.render import sky_basis
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.spherical_harmonics import sph_harm

__all__ = ["principal_axes", "projected_stretch", "sphere_coefficients",
           "affine_visibility"]


def principal_axes(field, first_order=False):
    """The ellipsoid's principal directions and the radii along them.

    The directions are the **eigenvectors of the quadrupole tensor**

        M = (3/2) (eps_tide t t^T - eps_spin s s^T) ,

    which is what the depth's quadratic form is.  Taking them from the tide axis instead
    is wrong whenever the tide is weak or absent: with eps_tide = 0 the tide's direction
    has no effect on the potential, yet it would still rotate the frame, and on an oblate
    body the level radius along a tilted direction differs from the equatorial one --- so
    the ellipsoid would change for no physical reason.

    The radii are then measured along those directions, and are *not* taken from the
    tensor, so the degeneracy of the equatorial pair costs nothing.

    With ``first_order`` the radii are the linearised ones, ``1 + delta r`` with
    ``delta r = -phi_1``.  That is closed form --- no root-find at all --- and differs from
    the exact radii only at ``O(eps^2)``, below the affine model's own ``O(eps)`` shape
    error, so it costs nothing.  It makes the first-order row of the model fully analytic.

    Returns ``(frame, radii)`` with the directions as the frame's *columns*.
    """
    tide, spin = field.axes()
    tensor = 1.5 * (field.tide * np.outer(tide, tide)
                    - field.centrifugal * np.outer(spin, spin))
    _, frame = np.linalg.eigh(tensor)
    if first_order:
        radii = np.array([1.0 + field.first_order_radii(direction)[0]
                          for direction in frame.T])
    else:
        radii = np.array([field.level_radius(direction) for direction in frame.T])
    return frame, radii


def projected_stretch(field, line_of_sight, first_order=False):
    """Lambda, the 2x2 stretch of the sky plane, in the sky's own ``(east, north)`` order.

    That order is what the image path uses, since ``sky_basis`` returns east and north and
    the samples' sky positions are their projections.  The frequencies are ``(u, v) =
    (North, East)``, so :func:`affine_visibility` reorders them before applying Lambda; one
    matrix cannot serve both orders.

    ``first_order`` takes the ellipsoid's radii from their linearisation instead of the
    level condition; see :func:`principal_axes`.
    """
    east, north = sky_basis(line_of_sight)
    projection = np.stack([east, north])                 # sky order
    frame, radii = principal_axes(field, first_order=first_order)
    shape = projection @ (frame * radii)                 # 2 x 3, columns scaled
    gram = shape @ shape.T                               # 2 x 2
    values, vectors = np.linalg.eigh(gram)
    return vectors @ np.diag(np.sqrt(np.maximum(values, 0.0))) @ vectors.T


def sphere_coefficients(field, l_max=16, power=1.0, n_samples=20000, eye=None,
                        linear=False):
    """The brightness field on the unit sphere, projected onto the harmonics.

    With ``a_lm = (4 pi / N) sum_i I(n_i) conj(Y_lm(n_i))`` and the golden-angle
    points carrying equal solid angles, this is an ordinary quadrature.  The result is
    stored half-complex (D13): ``a_l0`` real, ``a_lm`` complex for ``m > 0``.

    The projection is done in the **observer** frame when ``eye`` is given, so the
    coefficients need no rotation afterwards: ``theta`` is measured from the polar axis and
    ``phi`` from North toward East, which is the azimuth the mode sum uses (D23, D24).  The
    polar axis is not a free choice --- it must be ``north x east`` for that azimuth
    relation to hold --- and ``(north, east)`` come from ``sky_basis``, the same pair the
    direct sum projects onto.  An earlier version rotated body-frame coefficients instead;
    it left the pattern axisymmetric and its visibility real, and the frame error showed up
    as a quarter-turn between the real and imaginary parts.

    With these coefficients ``visibility(coeffs, u, v)`` is the sphere's visibility in
    closed form, which is what makes ``V_sphere(Lambda^T u)`` a short mode sum instead
    of a direct sum over samples --- the whole point of the affine model.

    ``linear`` projects the first-order pattern of :func:`render.first_order_brightness`
    instead of the exact one.  The two differ at ``O(eps^2)``, well below the affine
    model's own ``O(eps)`` shape error, so this is a diagnostic comparison rather than
    an accuracy improvement.
    """
    directions, brightness = render.sphere_brightness(field, power=power,
                                                      n_samples=n_samples,
                                                      linear=linear)
    if eye is None:
        north = np.array([1.0, 0.0, 0.0])
        east = np.array([0.0, 1.0, 0.0])
        polar = np.array([0.0, 0.0, 1.0])
    else:
        east, north = render.sky_basis(np.asarray(eye, dtype=float))
        polar = np.cross(north, east)
        polar = polar / np.linalg.norm(polar)
    theta = np.arccos(np.clip(directions @ polar, -1.0, 1.0))
    phi = np.arctan2(directions @ east, directions @ north)
    solid = 4.0 * np.pi / n_samples
    blocks = []
    for l in range(l_max + 1):
        block = np.zeros(l + 1, dtype=complex)
        for m in range(l + 1):
            block[m] = solid * np.sum(brightness
                                      * np.conjugate(sph_harm(l, m, theta, phi)))
        blocks.append(block)
    return Coeffs(tuple(blocks))


def affine_visibility(sphere_visibility, stretch, u, v):
    """``V_sphere(Lambda^T (u, v))``.

    ``stretch`` is in the sky's ``(east, north)`` order and ``(u, v)`` are ``(North, East)``
    (D24), so the frequency is reordered into the sky's order before Lambda acts on it.
    Passing ``(u, v)`` straight through transposes the stretch and mixes the axes, which
    is invisible for a nearly isotropic Lambda and wrong for an anisotropic one.

    ``u`` and ``v`` are scalars when ``sphere_visibility`` sums a phase over its samples,
    which is the case for the renderer's ``direct_visibility``; they may be arrays only if
    the sphere function itself broadcasts.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    # the frequency in the sky's order, then stretched, then back to (u, v) = (North, East)
    east, north = v, u
    stretched_east = stretch[0, 0] * east + stretch[1, 0] * north
    stretched_north = stretch[0, 1] * east + stretch[1, 1] * north
    return sphere_visibility(stretched_north, stretched_east)
