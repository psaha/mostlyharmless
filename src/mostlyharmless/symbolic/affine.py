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
           "first_order_coefficients", "affine_visibility"]


def principal_axes(field):
    """The ellipsoid's principal directions and the radii along them.

    The directions are the **eigenvectors of the quadrupole tensor**

        M = (3/2) (eps_tide t t^T - eps_spin s s^T) ,

    which is what the depth's quadratic form is.  Taking them from the tide axis instead
    is wrong whenever the tide is weak or absent: with eps_tide = 0 the tide's direction
    has no effect on the potential, yet it would still rotate the frame, and on an oblate
    body the level radius along a tilted direction differs from the equatorial one --- so
    the ellipsoid would change for no physical reason.

    The radii are then measured along those directions, and are *not* taken from the
    tensor, so the degeneracy of the equatorial pair costs nothing.  They come from the level
    condition, solved directly.

    Linearising them instead, ``1 + delta r`` with ``delta r = -phi_1``, removes the last
    root-find and is closed form --- but it is exact only to first order, and the next term is
    not small near the closure limit: at the amplitude ceiling it gives a bulge two fifths of
    the true one.  Since the saving is three bisections per parameter set, some fifty scalar
    evaluations, it buys nothing measurable and is not worth the hazard.  Linearising the
    *brightness* is a different matter and is what the first-order row of the model does; that
    is exact at first order with no such divergence.

    Returns ``(frame, radii)`` with the directions as the frame's *columns*.
    """
    tide, spin = field.axes()
    tensor = 1.5 * (field.tide * np.outer(tide, tide)
                    - field.centrifugal * np.outer(spin, spin))
    _, frame = np.linalg.eigh(tensor)
    radii = np.array([field.level_radius(direction) for direction in frame.T])
    return frame, radii


def projected_stretch(field, line_of_sight):
    """Lambda, the 2x2 stretch of the sky plane, in the sky's own ``(east, north)`` order.

    That order is what the image path uses, since ``sky_basis`` returns east and north and
    the samples' sky positions are their projections.  The frequencies are ``(u, v) =
    (North, East)``, so :func:`affine_visibility` reorders them before applying Lambda; one
    matrix cannot serve both orders.
    """
    east, north = sky_basis(line_of_sight)
    projection = np.stack([east, north])                 # sky order
    frame, radii = principal_axes(field)
    shape = projection @ (frame * radii)                 # 2 x 3, columns scaled
    gram = shape @ shape.T                               # 2 x 2
    values, vectors = np.linalg.eigh(gram)
    return vectors @ np.diag(np.sqrt(np.maximum(values, 0.0))) @ vectors.T


def first_order_coefficients(field, eye, limb=0.0):
    """The first-order brightness pattern's harmonics, in closed form (eq. 82).

    To first order the pattern is a constant plus two quadrupoles, so the addition theorem
    gives its coefficients directly:

        a_00 = (1 - 4 eps_spin) sqrt(4 pi) ,
        a_2m = -(16 pi / 5) [ eps_tide conj(Y_{2m}(t_hat))
                              - eps_spin conj(Y_{2m}(s_hat)) ] ,

    every other coefficient zero --- no odd order at all, since the pattern is even under
    n -> -n.  The two axes are referred to the same observer frame as a sampled direction:
    theta from ``north x east``, phi from North towards East.  That frame rule is not
    optional; it is what makes these coefficients the ones :func:`visibility` expects.

    ``limb`` adds the same law the samplers use, ``1 - u_limb (1 - mu)`` with ``mu`` the cosine
    to the line of sight.  Since ``mu`` is a dipole about the observer's axis, the factor is a
    constant plus an ``l = 1`` term, which costs one coefficient at ``m = 0``.

    That factor *multiplies* the pattern, so strictly the quadrupole times ``mu`` couples into
    ``l = 1`` and ``l = 3`` at ``O(eps limb)``.  Those terms are dropped: the row is a model
    chosen for algebra, not a truncation of a physical law, and dropping them keeps it a sum
    of a constant, a dipole and two quadrupoles.  The leading ``l = 1`` term is what the eye
    reads as limb darkening, and at zero deformation the row is exact.

    No quadrature: ten harmonic evaluations and one constant, against ``n_samples``
    evaluations for the projection in :func:`sphere_coefficients`.  Certificate V21 checks
    the two agree, so this is the fast path and that is the reference.
    """
    east, north = sky_basis(np.asarray(eye, dtype=float))
    polar = np.cross(north, east)
    polar = polar / np.linalg.norm(polar)

    def angles(direction):
        return (float(np.arccos(np.clip(direction @ polar, -1.0, 1.0))),
                float(np.arctan2(direction @ east, direction @ north)))

    tide_axis, spin_axis = field.axes()
    blocks = [np.zeros(l + 1, dtype=complex) for l in range(3)]
    blocks[0][0] = 1.0 - 4.0 * field.centrifugal
    if limb:
        # I = (1 - limb) + limb mu, and mu = sqrt(4 pi / 3) Y_10 about the observer's axis, so
        # the limb is a constant plus a dipole and nothing else to leading order.
        blocks[0][0] *= 1.0 - limb
        blocks[1][0] += limb * np.sqrt(1.0 / 3.0)
    for m in range(3):
        tide_term = sph_harm(2, m, *angles(tide_axis))
        spin_term = sph_harm(2, m, *angles(spin_axis))
        blocks[2][m] = -(4.0 / 5.0) * (
            field.tide * np.conjugate(tide_term)
            - field.centrifugal * np.conjugate(spin_term))
    return Coeffs(tuple(blocks))


def sphere_coefficients(field, l_max=16, power=1.0, n_samples=20000, eye=None,
                        linear=False):
    """The brightness field on the unit sphere, projected onto the harmonics.

    With ``a_lm = (1 / N) sum_i I(n_i) conj(Y_lm(n_i))`` and the golden-angle points
    carrying equal solid angles, this is an ordinary quadrature --- the two 4 pi, from the
    projection and from the solid angle, cancel in this convention.  The result is
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
    # with the average of |Y|^2 over the sphere equal to one, the projection is
    # a_lm = (1/4 pi) int I conj(Y), and the golden-angle samples carry 4 pi / N each, so
    # the two 4 pi cancel and this is a plain average over the samples
    solid = 1.0 / n_samples
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
