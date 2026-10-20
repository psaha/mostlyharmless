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
    "jupiter_coefficients",
    "triaxial_field",
    "sky_map",
    "figure_map",
    "figure_uv",
    "figure_convergence",
    "figure_jupiter",
    "figure_triaxial",
    "figure_stretched_sphere",
    "figure_visibility_maps",
    "figure_delta_visibility",
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


# A Jupiter caricature.  The bands are a function of latitude alone, so they are
# zonal and contribute nothing as the planet turns; the spot is what moves.
JUPITER_BAND_WIDTH = 15.0            # degrees between like bands
JUPITER_POLAR_DARKENING = 0.30
JUPITER_BAND_CONTRAST = 0.35
JUPITER_SPOT_LATITUDE = -22.0        # degrees, planetocentric
JUPITER_SPOT_RADIUS = 0.115          # radians, about the Great Red Spot's size
JUPITER_SPOT_DEPTH = 0.45


def _banded_moments(l_max, n_nodes=400):
    """Legendre moments of the banded latitude profile."""
    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)
    latitude = np.arcsin(nodes)                       # planetocentric, radians
    band = np.cos(2.0 * np.pi * np.degrees(latitude) / JUPITER_BAND_WIDTH)
    profile = ((1.0 - JUPITER_POLAR_DARKENING * np.sin(latitude) ** 2)
               * (1.0 + JUPITER_BAND_CONTRAST * band))
    return [float(np.sum(weights * profile * eval_legendre(l, nodes)))
            for l in range(l_max + 1)]


def jupiter_coefficients(l_max=20, spot_latitude=JUPITER_SPOT_LATITUDE,
                         spot_radius=JUPITER_SPOT_RADIUS,
                         spot_depth=JUPITER_SPOT_DEPTH):
    """A Jupiter caricature: banded zones and belts, plus one oval spot.

    Spherical throughout.  Jupiter's six per cent flattening is deliberately absent,
    non-spherical sources being Stage 6.

    The bands depend on latitude alone, so they are zonal harmonics, converted from
    Legendre moments by a_l0 = c_l sqrt(pi (2l+1)) and scaled to a_00 = 1.  The spot
    is the polar cap of :func:`cap_coefficients`, rotated to its latitude and with
    its uniform part dropped; that rotation is the only source of m != 0 content in
    the map, and so the only thing that changes as the planet turns.
    """
    moments = _banded_moments(l_max)
    zeroth = moments[0] * np.sqrt(np.pi)
    total = []
    for l in range(l_max + 1):
        block = np.zeros(l + 1, dtype=complex)
        block[0] = moments[l] * np.sqrt(np.pi * (2 * l + 1)) / zeroth
        total.append(block)

    # The cap's own l = 0 part is its uniform background, so drop it and keep the
    # localised deficit.  Its rotation lands the spot at body longitude 90 degrees.
    cap = cap_coefficients(l_max, half_angle=spot_radius, depth=spot_depth)
    placed = rotate_coefficients(cap, 0.0, np.radians(90.0 - spot_latitude), 0.0)
    for l, block in enumerate(placed.blocks):
        contribution = np.array(block)
        contribution[0] = 0.0
        total[l] = total[l] + contribution

    scale = total[0][0].real
    return Coeffs(tuple(block / scale for block in total))


def figure_jupiter(target, l_max=20, n=241, phase=0.0):
    """The caricature as a map, and as the observer sees it at two spin phases.

    Top row down the spin axis: the bands are rings about the pole, and the spot sits
    at its own longitude.  Bottom row in the observer's equatorial view, where the
    bands are horizontal and only the spot moves.
    """
    plt = _pyplot()
    body = jupiter_coefficients(l_max)
    first = rotate_coefficients(body, phase, 0.5 * np.pi, 0.0)
    second = rotate_coefficients(body, phase + 0.5 * np.pi, 0.5 * np.pi, 0.0)
    panels = ((MEAN_FACTOR * sky_map(body, n=n), "down the spin axis, near side"),
              (MEAN_FACTOR * sky_map(body, back=True, n=n),
               "down the spin axis, far side"),
              (MEAN_FACTOR * sky_map(first, n=n), "observer's view, phase 0"),
              (MEAN_FACTOR * sky_map(second, n=n),
               "observer's view, phase 90 degrees"))
    limits = (float(np.nanmin(panels[0][0])), float(np.nanmax(panels[0][0])))
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 6.6), constrained_layout=True)
    for axis, (data, title) in zip(axes.ravel(), panels):
        drawn = axis.imshow(data, origin="lower", extent=(1, -1, -1, 1),
                            cmap="magma", vmin=limits[0], vmax=limits[1])
        axis.set_title(title, fontsize=9)
        axis.set_xlabel(r"$y$", fontsize=9)
    for axis in (axes[0, 0], axes[1, 0]):
        axis.set_ylabel(r"$x$", fontsize=9)
    fig.colorbar(drawn, ax=axes, shrink=0.8, label=r"$I/\langle I\rangle$")
    fig.savefig(target)
    plt.close(fig)


# The Stage 6 star: a tidal bulge at azimuth TRIAXIAL_PHASE in the equatorial plane,
# and a spin flattening along z, at a 3:1 ratio.  That ratio makes the perturbation
# vanish at the tide's perpendicular, so that semi-axis sits at exactly the reference
# radius.  At these amplitudes the closure combination is
# eps_tide + 1.5 eps_cent = 0.144 against the limit of 0.148 --- 97 per cent of the
# way --- so the surface is close to pinching along the tide axis and the projection
# weight r^2/(n.N) grows sharply toward the bulge tip.  Sampling artefacts and model
# error have to be told apart there.
# The phase matters: with the tide *in* the plane of the spin and the line of sight,
# reflection through that plane fixes both quadrupoles, so the image keeps the mirror
# symmetry of an axisymmetric star however triaxial the body is.  Only a phase away
# from that plane breaks it, which is what the null test of section 8 looks for.
# Set TRIAXIAL_TIDE = 0.0 for the pure-spin case, which is axisymmetric.
TRIAXIAL_TIDE = 0.096
TRIAXIAL_SPIN = 0.032
TRIAXIAL_PHASE = np.radians(45.0)
TRIAXIAL_INCLINATION = np.radians(60.0)


def triaxial_field(tide=TRIAXIAL_TIDE, spin=TRIAXIAL_SPIN, phase=TRIAXIAL_PHASE):
    """The Stage 6 star: a tidal quadrupole in the equatorial plane, a spin one along z.

    ``phase`` puts the tidal axis at that azimuth from the sky's east direction, and is
    what controls whether the mirror symmetry is broken at all.
    """
    return Deformation(tide=tide,
                       tide_direction=(np.cos(phase), np.sin(phase), 0.0),
                       centrifugal=spin, spin_direction=(0.0, 0.0, 1.0))


def figure_triaxial(target, tide=TRIAXIAL_TIDE, spin=TRIAXIAL_SPIN,
                    phase=TRIAXIAL_PHASE, n_samples=150000, n=241):
    """The deformed star as seen down the spin axis, and inclined, near and far.

    The surface is an equipotential and the brightness the bolometric von Zeipel law
    I ∝ g, in units where the undeformed g is 1, so all four panels share one colour
    scale.  Gravity darkening makes the poles brighter than the equator.  Limb
    darkening is absent by construction until Stage 7, and the caption says so.

    Both hemispheres are drawn through the *same* sky axes for a given orientation, so
    the far panel is the near view rotated by half a turn: the far observer's offset
    (a, b) is the antipodal surface point, not the mirrored one.  A convex body's
    silhouette is the same from either side, so the two panels of a view share an
    outline and differ only in the shading -- which is where the gravity darkening
    lives, and so where the triaxiality shows.
    Down the spin axis the flattening is along the line of sight and only the tidal
    bulge shapes the outline; inclined, both show.
    """
    plt = _pyplot()
    field = triaxial_field(tide, spin, phase)
    views = ((np.array([0.0, 0.0, 1.0]), "down the spin axis"),
             (np.array([np.sin(TRIAXIAL_INCLINATION), 0.0,
                        np.cos(TRIAXIAL_INCLINATION)]),
              "inclined by 60 degrees"))
    panels = []
    for eye, title in views:
        east, north = render.sky_basis(eye)
        for side, sense in (("near side", 1.0), ("far side, from behind", -1.0)):
            direction = sense * eye
            points, brightness, weights = render.surface_samples(
                field, direction, n_samples=n_samples)
            panels.append((render.sky_image(points, brightness, weights,
                                            east, north, half_width=1.3, n=n),
                           f"{title}, {side}"))

    limits = (min(float(np.nanmin(panel[0])) for panel in panels),
              max(float(np.nanmax(panel[0])) for panel in panels))
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 6.8), constrained_layout=True)
    for axis, (data, title) in zip(axes.ravel(), panels):
        drawn = axis.imshow(data, origin="lower", extent=(1.3, -1.3, -1.3, 1.3),
                            cmap="magma", vmin=limits[0], vmax=limits[1])
        axis.set_title(title, fontsize=8)
        axis.set_xlabel(r"$y$", fontsize=9)
        axis.set_aspect("equal")
    for axis in (axes[0, 0], axes[1, 0]):
        axis.set_ylabel(r"$x$", fontsize=9)
    fig.colorbar(drawn, ax=axes, shrink=0.82, label=r"$g_\mathrm{eff}/g_0$")
    fig.savefig(target)
    plt.close(fig)


DELTA_Q_MAX = 2.0
DELTA_GRID = 31
# The truncation of the sphere's brightness pattern, projected onto the harmonics for
# the analytic affine visibility.  Keeping first order in eps makes the pattern exactly
# l = 0 and l = 2: the l = 4 content first appears at eps^2, from P_2^2.  So l_max = 2
# is *exact* at that order rather than a truncation, and the only error left is the
# model's own -- the ellipsoid replacing the true level surface.
SPHERE_L_MAX = 2


def figure_stretched_sphere(target, tide=TRIAXIAL_TIDE, spin=TRIAXIAL_SPIN,
                            phase=TRIAXIAL_PHASE, inclination=TRIAXIAL_INCLINATION,
                            power=1.0, n_samples=150000, n=241):
    """The actual rendering against the affine model's stretched sphere.

    Rows are the two cases and columns the two hemispheres: the affine prediction
    above, the actual star below, matching ``figure_visibility_maps``, and the
    observer's side and the far side.  The prediction is the brightness field of
    the *same* potential evaluated on the unit sphere, then stretched by ``Lambda`` --
    the samples moved to ``Lambda x`` and the brightness divided by ``|det Lambda|``,
    which keeps the total flux while the area grows.

    Both cases are shown as the *fractional* departure from their own mean, which is
    D22 applied in image space: the affine model conserves total flux while the render
    does not (it sets I ∝ g with no renormalisation), so the two would otherwise differ
    by a uniform factor of order the area ratio that swamps the pattern.  *Dividing* by
    the mean removes that factor as well as the level; *subtracting* it would leave the
    model's pattern still short by the same factor, since the model's mean is the
    sphere's divided by |det Lambda|.  The zero contour is the locus g = <g> either way.
    """
    plt = _pyplot()
    field = triaxial_field(tide, spin, phase)
    eye = np.array([np.sin(inclination), 0.0, np.cos(inclination)])
    east, north = render.sky_basis(eye)
    stretch = affine.projected_stretch(field, eye)
    scale = abs(float(np.linalg.det(stretch)))

    def contrast(image):
        """The image's fractional departure from its mean over the disk."""
        return image / np.nanmean(image) - 1.0

    panels = []
    for side, sense in (("observer's side", 1.0), ("far side", -1.0)):
        direction = sense * eye
        truth = render.surface_samples(field, direction, n_samples=n_samples,
                                       power=power)
        sphere = render.sphere_samples(field, direction, n_samples=n_samples,
                                       power=power)
        # The affine prediction: the sphere's image stretched in the sky plane.  The
        # samples move to Lambda x, the area each one represents grows by |det Lambda|,
        # and the brightness falls by the same factor -- so each sample carries its flux
        # and the total is unchanged, while the surface brightness dims.
        sky = np.stack([sphere[0] @ east, sphere[0] @ north], axis=-1)
        moved = sky @ stretch.T
        predicted_points = moved[:, :1] * east + moved[:, 1:] * north
        # the affine row first, matching figure_visibility_maps
        panels.append(((contrast(render.sky_image(predicted_points, sphere[1] / scale,
                                                  sphere[2] * scale, east, north,
                                                  half_width=1.5, n=n)),
                        contrast(render.sky_image(*truth, east, north,
                                                  half_width=1.5, n=n))),
                       side))

    values = np.concatenate([panel.ravel()
                             for pair, _ in panels for panel in pair])
    limit = float(np.nanmax(np.abs(values))) or 1e-30

    fig, axes = plt.subplots(2, 2, figsize=(6.8, 6.9), constrained_layout=True)
    for column, (pair, side) in enumerate(panels):
        for row, (image, label) in enumerate(zip(pair, ("affine, rung 1",
                                                        "the star"))):
            drawn = axes[row, column].imshow(
                image, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit,
                extent=(1.5, -1.5, -1.5, 1.5))
            axes[row, column].set_aspect("equal")
            axes[row, column].set_title(label, fontsize=9)
            if row == 1:
                axes[row, column].set_xlabel(r"$y$", fontsize=9)
            if column == 0:
                axes[row, column].set_ylabel(side + "\n" + r"$x$", fontsize=8)
    fig.colorbar(drawn, ax=axes, shrink=0.85,
                 label=r"$g_\mathrm{eff}/\langle g_\mathrm{eff}\rangle - 1$")
    fig.savefig(target)
    plt.close(fig)


def figure_delta_visibility(target, tide=TRIAXIAL_TIDE, spin=TRIAXIAL_SPIN,
                            phase=TRIAXIAL_PHASE, inclination=TRIAXIAL_INCLINATION,
                            n_samples=20000, q_max=DELTA_Q_MAX, grid=DELTA_GRID):
    """Delta V against the volume-equivalent uniform disk, numerical and affine.

    Visibilities are normalised to 1 at the origin, as the observable is (D22), so the
    difference vanishes there by construction.

    Rows: the numerical result, the affine model at the same points, and their
    difference -- the model error, which is Stage 6's gate.  Columns: Re Delta V,
    Im Delta V, and the *violation* of the mirror symmetry,

        | V(u, v) - conjugate(V(-u, v)) | ,

    which vanishes identically for an axisymmetric star.  Reflection in the plane of
    the spin and the line of sight acts on the visibility as
    V(-u, v) = conjugate(V(u, v)), the projected spin lying along v, so this violation
    *is* the null test of section 8.  It is formed from the normalised visibilities and
    the reference cancels: the uniform disk's visibility is real and depends on q alone.

    The pattern on the spherical counterpart is carried *radially*, by direction from
    the centre.  Carrying it along the stretch instead would make the model exact by
    construction -- the change of variables is then an identity -- and so vacuous.  The
    radial correspondence is a genuine model, and its phase error is part of what the
    error row shows.
    """
    plt = _pyplot()
    field = triaxial_field(tide, spin, phase)
    eye = np.array([np.sin(inclination), 0.0, np.cos(inclination)])
    east, north = render.sky_basis(eye)
    points, brightness, weights = render.surface_samples(field, eye,
                                                         n_samples=n_samples)

    # The spherical counterpart: the same brightness by direction, on the unit sphere.
    directions = points / np.linalg.norm(points, axis=1)[:, None]
    sphere_weights = np.maximum(directions @ eye, 0.0) * (4.0 * np.pi / n_samples)

    equivalent = field.equivalent_radius()
    stretch = affine.projected_stretch(field, eye)

    def deformed_sum(u, v):
        return render.direct_visibility(points, brightness, weights,
                                        east, north, u, v)

    def sphere_sum(u, v):
        return render.direct_visibility(directions, brightness, sphere_weights,
                                        east, north, u, v)

    deformed_flux = deformed_sum(0.0, 0.0)
    spherical_flux = sphere_sum(0.0, 0.0)

    def numerical(u, v):
        return deformed_sum(u, v) / deformed_flux

    def modelled(u, v):
        return affine.affine_visibility(sphere_sum, stretch, u, v) / spherical_flux

    def reference(u, v):
        # already normalised to 1 at the origin, which is the convention of D22
        return visibility_uniform_disk(u, v, equivalent)

    axis = np.linspace(-q_max, q_max, grid)

    def panels(visibility_function):
        """(Re Delta V, Im Delta V, violation) on the grid, rows v and columns u.

        The sums take scalar frequencies, so this is a loop; the reference is real, so
        it contributes to the real part only.
        """
        shape = (grid, grid)
        real = np.empty(shape)
        imaginary = np.empty(shape)
        violation = np.empty(shape)
        for row, b in enumerate(axis):                  # v
            for column, a in enumerate(axis):           # u
                here = visibility_function(a, b)
                mirrored = visibility_function(-a, b)
                real[row, column] = here.real - reference(a, b)
                imaginary[row, column] = here.imag
                violation[row, column] = abs(here - np.conj(mirrored))
        return real, imaginary, violation

    numerical_panels = panels(numerical)
    affine_panels = panels(modelled)
    error_panels = tuple(one - two for one, two in zip(numerical_panels,
                                                       affine_panels))

    rows = ((numerical_panels, "numerical"),
            (affine_panels, "affine, rung 1"),
            (error_panels, "numerical minus affine"))
    headers = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$",
               r"$|V - \overline{V(-u,v)}|$")

    fig, axes = plt.subplots(3, 3, figsize=(8.4, 8.2), constrained_layout=True)
    for column, (header, non_negative) in enumerate(zip(headers,
                                                        (False, False, True))):
        values = np.concatenate([data[column].ravel() for data, _ in rows])
        if non_negative:
            lower, upper = 0.0, float(np.max(values)) or 1e-30
            colour_map = "magma"
        else:
            lower = -float(np.max(np.abs(values))) or -1e-30
            upper = -lower
            colour_map = "RdBu_r"
        drawn = None
        for row, (data, row_label) in enumerate(rows):
            drawn = axes[row, column].imshow(
                data[column], origin="lower", cmap=colour_map, vmin=lower,
                vmax=upper, extent=(-q_max, q_max, -q_max, q_max))
            axes[row, column].set_aspect("equal")
            if row == 0:
                axes[row, column].set_title(header, fontsize=9)
            if row == 2:
                axes[row, column].set_xlabel(r"$u$", fontsize=9)
            if column == 0:
                axes[row, column].set_ylabel(row_label + "\n" + r"$v$",
                                             fontsize=8)
        fig.colorbar(drawn, ax=axes[:, column], shrink=0.85)
    fig.savefig(target)
    plt.close(fig)


def figure_visibility_maps(target, tide=TRIAXIAL_TIDE, spin=TRIAXIAL_SPIN,
                           phase=TRIAXIAL_PHASE, inclination=TRIAXIAL_INCLINATION,
                           power=1.0, n_samples=20000, q_max=2.0, grid=31):
    """Re and Im of the normalised visibility: the affine model against the star.

    Row 1 is the affine prediction from rung 1, ``V_sphere(Lambda^T u)``; row 2 is the
    direct sum over the actual surface.  Columns are the real and imaginary parts, and
    what is drawn is the difference from a *unit sphere* --- the flat disk of unit
    radius, whose normalised visibility is the analytic ``2 J_1(x)/x``:

        Delta V = V / V(0,0) - V_unit-sphere / V_unit-sphere(0,0) ,

    so both rows vanish at the origin and the panels show the deformation's signature
    alone rather than the lobe structure that all three share.  The reference is the
    *unit* sphere, not the volume-equivalent one, so a small scale mismatch remains
    (the star's mean radius is about 1.04); using the volume-equivalent disk would
    isolate the deformation alone.

    The largest difference between the rows, relative to the origin, is written on the
    figure.  It is Stage 6's gate, and it is the one thing these panels cannot show by
    eye: the disagreement is small in absolute terms even where it is clear in shape.
    The reference cancels in that number, so it is unchanged by the subtraction.

    Axes follow D24: u is North and runs vertically, v is East and increases leftward.
    """
    plt = _pyplot()
    field = triaxial_field(tide, spin, phase)
    eye = np.array([np.sin(inclination), 0.0, np.cos(inclination)])
    east, north = render.sky_basis(eye)
    points, brightness, weights = render.surface_samples(field, eye,
                                                         n_samples=n_samples,
                                                         power=power)
    stretch = affine.projected_stretch(field, eye)

    def actual(u, v):
        return render.direct_visibility(points, brightness, weights,
                                        east, north, u, v)

    # The sphere's visibility in closed form: its pattern projected onto the harmonics,
    # so one point costs a short mode sum instead of a sum over samples.
    coefficients = affine.sphere_coefficients(field, l_max=SPHERE_L_MAX,
                                              power=power, n_samples=n_samples,
                                              inclination=inclination)
    sphere_flux = complex(visibility(coefficients, 0.0, 0.0))

    def sphere_sum(u, v):
        return visibility(coefficients, u, v) / sphere_flux

    def modelled(u, v):
        return affine.affine_visibility(sphere_sum, stretch, u, v)

    actual_flux = actual(0.0, 0.0)

    axis = np.linspace(-q_max, q_max, grid)
    # the outer loop is u, the vertical axis, so rows are u and columns v
    u_grid, v_grid = np.meshgrid(axis, axis, indexing="ij")
    reference = visibility_uniform_disk(u_grid, v_grid, 1.0)
    affine_values = np.array([[modelled(u, v) for v in axis]
                              for u in axis]) - reference
    actual_values = np.array([[actual(u, v) / actual_flux for v in axis]
                              for u in axis]) - reference
    worst = float(np.max(np.abs(affine_values - actual_values)))

    rows = ((affine_values, "affine, rung 1"), (actual_values, "the star"))
    headers = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$")
    values = np.concatenate([data.ravel() for data, _ in rows])
    limit = float(np.max(np.abs(values))) or 1e-30

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 7.4), constrained_layout=True)
    for column, header in enumerate(headers):
        part = (lambda data: data.real) if column == 0 else (lambda data: data.imag)
        for row, (data, label) in enumerate(rows):
            drawn = axes[row, column].imshow(
                part(data), origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit,
                extent=(q_max, -q_max, -q_max, q_max))
            axes[row, column].set_aspect("equal")
            if row == 0:
                axes[row, column].set_title(header, fontsize=10)
            if row == 1:
                axes[row, column].set_xlabel(r"$v$ (East)", fontsize=9)
            if column == 0:
                axes[row, column].set_ylabel(label + "\n" + r"$u$ (North)",
                                             fontsize=8)
        fig.colorbar(drawn, ax=axes[:, column], shrink=0.85)
    fig.text(0.5, -0.01, rf"$\max|V_\mathrm{{affine}} - V_\mathrm{{star}}| = "
             rf"{worst:.4f}$", ha="center", fontsize=9)
    fig.savefig(target)
    plt.close(fig)


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
        drawn = axis.imshow(data, origin="lower", extent=(1, -1, -1, 1),
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
        drawn = axis.imshow(data, origin="lower",
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
        drawn = axis.imshow(data, origin="lower", extent=(1, -1, -1, 1),
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
                           ("convergence.pdf", figure_convergence),
                           ("stretched-sphere.pdf", figure_stretched_sphere),
                           ("visibility-maps.pdf", figure_visibility_maps)):
        path = directory / name
        function(path)
        written.append(path)
    return written


if __name__ == "__main__":
    for written in write():
        print(written)
