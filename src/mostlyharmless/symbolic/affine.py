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

from mostlyharmless.symbolic.render import sky_basis

__all__ = ["body_axes", "projected_stretch", "affine_visibility"]


def body_axes(field):
    """The level radii along the coordinate axes, in the field's own frame."""
    return np.array(field.axis_lengths(), dtype=float)


def projected_stretch(field, line_of_sight):
    """Lambda, the 2x2 stretch of the sky plane, from the axes and the view."""
    east, north = sky_basis(line_of_sight)
    projection = np.stack([east, north])                 # 2 x 3
    shape = projection @ np.diag(body_axes(field))       # 2 x 3
    gram = shape @ shape.T                               # 2 x 2
    values, vectors = np.linalg.eigh(gram)
    return vectors @ np.diag(np.sqrt(np.maximum(values, 0.0))) @ vectors.T


def affine_visibility(sphere_visibility, stretch, u, v):
    """``V_sphere(stretch^T (u, v))``.

    The stretch is applied to the frequency, not the image.  ``u`` and ``v`` are
    scalars when ``sphere_visibility`` sums a phase over its samples, which is the
    case for the renderer's ``direct_visibility``; they may be arrays only if the
    sphere function itself broadcasts.
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    return sphere_visibility(stretch[0, 0] * u + stretch[1, 0] * v,
                             stretch[0, 1] * u + stretch[1, 1] * v)
