"""The interactive front end: sliders and toggles over the Stage 6 model.

Stage 8, D25.  Matplotlib widgets only --- no Jupyter, and nothing added to the library.

Two radio buttons: the *view* is a rendering or a pair of ``(u,v)`` maps, and the *model* has
four settings, so that the shape approximation, the brightness approximation and the choice of
coefficients can be varied separately.  ``numeric`` is the exact level surface; ``fitted`` keeps
the ellipsoid's stretch but fits the coefficients at ``l <= 4`` instead of projecting the
pattern (section 8.8 of the manual); ``affine`` carries the exact brightness pattern at
``l <= 2``; ``first-order`` also replaces the pattern by its linearised form, equation (74).  The sliders cover
the five parameters.  The brightness law itself is fixed at von Zeipel, ``I ∝ g`` (D21), and
the position angle is pinned to zero as a gauge, so North is vertical.

The closure bound is explicit.  Increasing the tide's latitude *loosens* it --- at the
equator the centrifugal quadrupole deepens the potential by 1.5 eps_spin, while along the
spin it contributes nothing --- so the infeasible region is near latitude zero and the
latitude slider carries a *lower* limit.  The eps sliders get upper limits from the same
function.  Clamping is applied after every change rather than by restricting ranges, since
the feasible set is a polygon in the amplitudes rather than a box.

``compute_image`` and ``compute_visibility`` hold the whole numerical content and are
testable without a widget; ``build`` only wires sliders to them.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq

from mostlyharmless.symbolic import affine, fitted, gui_core, render
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk

__all__ = ["Parameters", "closure_depth", "largest_spin", "smallest_latitude",
           "clamped", "compute_image", "compute_visibility", "fitted_report", "build",
           "main"]

# D25: the analytic row keeps l <= 2, the eps^2 truncation of the l >= 4 content with the
# l = 2 coefficient resummed.  Not a choice to make.
SPHERE_L_MAX = 2
# The fitted row's truncation (manual section 8.8).  Wider than the analytic row's two,
# because the point of it is to measure what the truncation costs.
FIT_L_MAX = 4
# The model radio, in order.  The fitted row sits last, below first-order, because it is the
# one that replaces the coefficients rather than the shape or the brightness.
MODEL_LABELS = ("numeric", "first-order", "fitted, l = 0,2,4",
                "fitted, l = 0,2,4,6,8")
# The truncation each fitted row uses.  Two rows rather than a slider: the point is to see
# what the extra order buys, and every change costs a full redraw.
FIT_ORDERS = {"fitted, l = 0,2,4": 4, "fitted, l = 0,2,4,6,8": 8}

# Named parameter sets, labelled by what they are for rather than by their angles.
#
#   reference  --- below the closure ceiling, tide and spin both on, moderate tilt: the
#                  configuration to quote when comparing rows.  eps_tide = 0.03 and
#                  eps_spin = 0.01 keep E = eps_tide + 1.5 eps_spin at about 0.045, a third
#                  of the 0.148 ceiling, and i = 60 deg shows both the tide's bulge and the
#                  spin's flattening.
#   tilt check --- the deformations at the ceiling, with the tide axis off both the sky's
#                  axes (lambda = -45, varphi = 45) so the bulge lies on a diagonal.  A
#                  mirrored or mis-rotated image is invisible for an axis-aligned bulge and
#                  obvious for this one, which is why it is the geometry check.
#   ceiling    --- the same amplitudes with the axis in the equatorial plane, so the
#                  deformation is as large as the model allows: the stress test for the
#                  analytic rows, where the affine error reaches 0.11.
PRESETS = {
    "reference": dict(tide=0.03, spin=0.01, inclination=np.radians(60.0),
                      latitude=0.0, longitude=0.0),
    "tilt check": dict(tide=0.138, spin=0.0, inclination=np.radians(60.0),
                       latitude=np.radians(-45.0), longitude=np.radians(45.0)),
    "ceiling": dict(tide=0.138, spin=0.0, inclination=np.radians(60.0),
                    latitude=0.0, longitude=0.0),
}
# The panel titles.  Module level rather than inline because mathtext is parsed at draw
# time, so an unsupported symbol crashes the draw rather than a callback --- and draw-time
# code cannot be covered by a headless test unless the strings are reachable.
MODEL_TITLES = {"numeric": "the star",
                "affine": "affine",
                "first-order": "affine, linear $g$",
                "fitted, l = 0,2,4": "fitted, l = 0,2,4",
                "fitted, l = 0,2,4,6,8": "fitted, l = 0,2,4,6,8"}
VISIBILITY_TITLES = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$")
CLOSURE = 0.148
# Numerical resolution: discrete levels, chosen by radio buttons, because every change
# costs a full redraw.  Each sets the surface samples, the image pixels and the (u, v) grid
# together, keeping the sample density per pixel roughly constant.  Only the samples are
# physics --- the pixels are display and the grid is the (u, v) sampling --- so the control
# trades speed against noise, which is what makes it the way to tell a sampling artefact
# from a model error.  The finest level is chosen to match the Stage 6 figures: 150 000
# samples for a rendering and a 31-point grid for a visibility map.
RESOLUTIONS = ((10000, 140, 17), (30000, 220, 23), (150000, 360, 31))
RESOLUTION_LABELS = ("coarse", "medium", "fine")
SAMPLES, PIXELS, GRID = RESOLUTIONS[0]
Q_MAX = 2.0


@dataclass(frozen=True)
class Parameters:
    """The five free parameters, in radians where they are angles.

    The GUI shows degrees and converts at the widget boundary; this class and everything
    outside the GUI stay in radians.  The spin phase is deliberately absent: this front
    end is for a spin and tide alone, and on a smooth map the phase is degenerate with
    the tide's longitude --- only their difference is observable.
    """

    tide: float = 0.05
    spin: float = 0.01
    inclination: float = np.radians(60.0)
    latitude: float = 0.0
    longitude: float = np.radians(45.0)
    # limb darkening, I = I(1 - limb rho^2) with rho the projected radius.  A viewing effect
    # rather than a property of the surface, so every row has to be told about it: the
    # numeric one samples it, the fitted ones absorb it by fitting those samples, and the
    # first-order one carries it as a quadrupole about the observer's axis.  One law for all
    # three, so the rows differ in the shape and the brightness pattern and not in the limb.
    limb: float = 0.0

    def field(self):
        """The Deformation: the tidal axis tilted by latitude, swung by longitude."""
        axis = (np.cos(self.latitude) * np.cos(self.longitude),
                np.cos(self.latitude) * np.sin(self.longitude),
                np.sin(self.latitude))
        return Deformation(tide=self.tide, tide_direction=axis,
                           centrifugal=self.spin, spin_direction=(0.0, 0.0, 1.0))


def closure_depth(tide, spin, latitude):
    """The largest depth over directions; the surface closes while this is below CLOSURE.

    Exact, not sampled.  With the tide's axis at latitude ``latitude`` in the ``x-z``
    plane and the spin along ``z``, the depth at ``r = 1`` is

        E(n) = (3/2) [ tide (n.t)^2 - spin (n.s)^2 ] + (3 spin - tide) / 2 ,

    a quadratic form whose maximum over unit ``n`` is its largest eigenvalue.  It lies in
    the span of the two axes, so with ``cos(theta) = t.s = sin(latitude)``

        E_max = (3/4)[(tide - spin) + sqrt((tide - spin)^2 + 4 tide spin cos^2 theta)]
                + (3 spin - tide) / 2 .

    The special cases follow: at latitude zero this is ``tide + 1.5 spin``, and at a
    polar tide (``cos theta = 0``) it is ``tide`` alone when the tide dominates, because
    the centrifugal term then has no gradient along the spin.
    """
    t, s = float(tide), float(spin)
    cosine_squared = np.cos(latitude) ** 2
    largest = 0.75 * ((t - s) + np.sqrt((t - s) ** 2 + 4.0 * t * s * cosine_squared))
    return float(largest + 0.5 * (3.0 * s - t))


def largest_spin(tide, latitude, limit=CLOSURE):
    """The largest ``eps_spin`` that still closes, by root solve.

    Monotone in the spin, but not linear --- the depth carries the spin inside a square
    root --- so this is a root solve rather than a division.
    """
    if closure_depth(tide, 0.0, latitude) > limit:
        return 0.0
    if closure_depth(tide, 1.0, latitude) <= limit:
        return 1.0
    return brentq(lambda spin: closure_depth(tide, spin, latitude) - limit,
                  0.0, 1.0, xtol=1e-12)


def smallest_latitude(tide, spin, limit=CLOSURE):
    """The *smallest* latitude that closes.  Latitude loosens the bound."""
    if closure_depth(tide, spin, np.pi / 2.0) > limit:
        return np.pi / 2.0                  # never closes
    if closure_depth(tide, spin, 0.0) <= limit:
        return 0.0                          # closes even at the equator
    return brentq(lambda lat: closure_depth(tide, spin, lat) - limit,
                  0.0, np.pi / 2.0, xtol=1e-4)


def clamped(parameters):
    """The nearest feasible parameters: spin down, then latitude up, until it closes."""
    low = smallest_latitude(parameters.tide, parameters.spin)
    if low >= np.pi / 2.0:
        return replace(parameters, tide=max(0.0, min(parameters.tide, 0.14)),
                       spin=0.0, latitude=parameters.latitude)
    latitude = parameters.latitude
    if abs(latitude) < low:
        # the bound depends on |latitude| only, since the depth is even in it
        latitude = low if latitude >= 0.0 else -low
    spin = min(parameters.spin, largest_spin(parameters.tide, latitude))
    return replace(parameters, latitude=latitude, spin=spin)


def _view(inclination):
    """The eye and the sky axes, with the spin phase absent: this model has none."""
    eye = np.array([np.sin(inclination), 0.0, np.cos(inclination)])
    return eye, render.sky_basis(eye)


def compute_image(parameters, model="numeric", samples=SAMPLES, pixels=PIXELS,
                  report=None):
    """The near and far hemispheres, as brightness in units of their own mean.

    ``report``, if given, is filled with the near side's fit when the row is the fitted
    one.  That is what the panel titles print: recomputing it would render the samples a
    second time, and would also report a slightly different fit from the one displayed.

    The affine image is the sphere's image stretched in the sky plane --- exact as a
    coordinate stretch, no per-point sum --- and the numeric one samples the level
    surface.  Both are mean-normalised so the flux is out of the way (D22).
    """
    field = parameters.field()
    eye, (east, north) = _view(parameters.inclination)
    linear = model == "first-order"
    stretch = affine.projected_stretch(field, eye)
    scale = abs(float(np.linalg.det(stretch)))

    def brightness(image):
        # in units of the mean, so the colourbar's black is zero brightness.  The panels used
        # to show the departure from the mean, which puts the floor at -1 and needs a
        # diverging map to read; a brightness is what the star actually has.
        return image / np.nanmean(image)

    panels = []
    for sense in (1.0, -1.0):
        direction = sense * eye
        if model.startswith("fitted"):
            # Lambda still comes from the level surface; only the coefficients are chosen.
            star = render.surface_samples(field, direction, n_samples=samples,
                                          limb=parameters.limb)
            fit = fitted.fit_image(star, stretch, east, north,
                                   l_max=FIT_ORDERS.get(model, FIT_L_MAX))
            if report is not None and sense > 0.0:
                report["fit"] = fit                     # the near side is the reported one
            image = fitted.model_image(fit.coefficients, stretch, half_width=1.5,
                                       n=pixels)
        elif model == "numeric":
            star = render.surface_samples(field, direction, n_samples=samples,
                                          limb=parameters.limb)
            image = render.sky_image(*star, east, north, half_width=1.5, n=pixels)
        else:
            sphere = render.sphere_samples(field, direction, n_samples=samples,
                                           linear=linear, limb=parameters.limb)
            sky = np.stack([sphere[0] @ east, sphere[0] @ north], axis=-1)
            moved = sky @ stretch.T
            points = moved[:, :1] * east + moved[:, 1:] * north
            image = render.sky_image(points, sphere[1] / scale, sphere[2] * scale,
                                     east, north, half_width=1.5, n=pixels)
        panels.append(brightness(image))
    return panels


def compute_visibility(parameters, model="numeric", samples=SAMPLES,
                       grid=GRID, q_max=Q_MAX):
    """Re and Im of Delta V against the unit sphere, for one model.

    Zero at the origin by construction, so the lobes the models share are gone and only
    the deformation's signature remains.  ``u`` is North, so it is the first index.
    """
    field = parameters.field()
    eye, (east, north) = _view(parameters.inclination)

    if model == "numeric":
        star = render.surface_samples(field, eye, n_samples=samples,
                                      limb=parameters.limb)
        flux = render.direct_visibility(*star, east, north, 0.0, 0.0)

        def model_visibility(u, v):
            return (render.direct_visibility(*star, east, north, u, v) / flux)
    else:
        if model.startswith("fitted"):
            star = render.surface_samples(field, eye, n_samples=samples,
                                          limb=parameters.limb)
            fitted_fit = fitted.fit_image(star, affine.projected_stretch(field, eye),
                                          east, north,
                                          l_max=FIT_ORDERS.get(model, FIT_L_MAX))
            coefficients = fitted_fit.coefficients
        elif model == "first-order":
            # the closed form of (83): no quadrature, ten harmonic evaluations.  The
            # ellipse is still the exact one, so only the brightness is linearised.
            coefficients = affine.first_order_coefficients(field, eye,
                                                           limb=parameters.limb)
        else:
            coefficients = affine.sphere_coefficients(
                field, l_max=SPHERE_L_MAX, n_samples=samples, eye=eye)
        sphere_flux = complex(visibility(coefficients, 0.0, 0.0))
        stretch = affine.projected_stretch(field, eye)

        def sphere_sum(u, v):
            return visibility(coefficients, u, v) / sphere_flux

        def model_visibility(u, v):
            return affine.affine_visibility(sphere_sum, stretch, u, v)

    axis = np.linspace(-q_max, q_max, grid)
    u_grid, v_grid = np.meshgrid(axis, axis, indexing="ij")
    reference = visibility_uniform_disk(u_grid, v_grid, 1.0)
    return np.array([[model_visibility(u, v) for v in axis] for u in axis]) - reference


def fitted_report(parameters, samples=SAMPLES, model="fitted, l = 0,2,4"):
    """The fitted row's diagnostics for the near side (manual section 8.8).

    The residual is the point of the row: with ``Lambda`` fixed it is the truncation floor
    given the ellipsoidal shape, so comparing it with the affine row's error separates the
    truncation from the shape model's own.  ``residual`` and ``residual_direct`` are the
    same quantity by two routes and must agree; ``condition`` is a property of the basis
    alone and does not depend on the parameters.
    """
    field = parameters.field()
    eye, (east, north) = _view(parameters.inclination)
    star = render.surface_samples(field, eye, n_samples=samples)
    return fitted.fit_image(star, affine.projected_stretch(field, eye), east, north,
                            l_max=FIT_ORDERS.get(model, FIT_L_MAX))


def build(parameters=None):  # pragma: no cover - needs a display
    """The widget window; returns (figure, state) with the live parameters in ``state``."""
    from matplotlib import pyplot as plt
    from matplotlib.widgets import RadioButtons

    state = {"parameters": clamped(parameters or Parameters()),
             "model": "numeric", "visibility": False, "refreshing": False,
             "resolution": 0, "limit_key": None, "limit": None, "colorbar": None}

    figure = plt.figure(figsize=(9.2, 6.6))
    axes, colorbar_axis = gui_core.panel_axes(figure)
    brightness_cmap = gui_core.image_colormap()

    def redraw(*_):
        """Mark the view stale and schedule one redraw: a drag fires on every motion."""
        coalescer.request()

    def _draw():
        current = state["parameters"]
        samples, pixels, grid = RESOLUTIONS[state["resolution"]]
        report = {}                     # filled by compute_image for the fitted row
        if state["visibility"]:
            values = compute_visibility(current, state["model"],
                                        samples=samples, grid=grid)
            panels = (values.real, values.imag)
            titles = VISIBILITY_TITLES
            extent = (Q_MAX, -Q_MAX, -Q_MAX, Q_MAX)
            labels = (r"$v$ (East)", r"$u$ (North)")
            bar_label = r"$\Delta V$"
        else:
            panels = compute_image(current, state["model"],
                                   samples=samples, pixels=pixels, report=report)
            # the analytic rows differ only in the brightness pattern: exact versus the
            # first-order form of equation (88)
            name = MODEL_TITLES[state["model"]]
            titles = (f"{name}, near", f"{name}, far")
            extent = (1.5, -1.5, -1.5, 1.5)
            labels = (r"$y$ (East)", r"$x$ (North)")
            bar_label = r"$g/\langle g\rangle$"

        # the flip and the shared scale are the scaffolding's, not this model's
        panels = tuple(gui_core.flip(panel) for panel in panels)
        limit = gui_core.colour_limit(
            (current, state["resolution"], state["visibility"]), state, panels)

        if "fit" in report:
            # the residual is the point of the row: at fixed Lambda it is the truncation
            # floor, so it is not a defect of the fit but the price of the degree
            titles = tuple(f"{title}  (residual {report['fit'].residual:.1e})"
                           for title in titles)

        # a brightness wants a linear grey ramp; a signed Delta V wants a diverging map
        brightness = not state["visibility"]
        cmap = brightness_cmap if brightness else "RdBu_r"

        drawn = None
        for axis_handle, data, title in zip(axes, panels, titles):
            axis_handle.clear()
            if brightness:
                gui_core.darken(axis_handle)     # after clear(), which resets the face
            drawn = axis_handle.imshow(
                data, origin="lower", cmap=cmap,
                vmin=0.0 if brightness else -limit, vmax=limit, extent=extent)
            axis_handle.set_title(title, fontsize=9)
            axis_handle.set_xlabel(labels[0], fontsize=8)
            axis_handle.set_ylabel(labels[1], fontsize=8)
        # one colourbar for both panels: they share a scale by construction.  Updated in
        # place rather than removed --- Colorbar.remove() takes its axes out of the
        # figure, so recreating on those axes gives a bar nothing draws to.
        if state["colorbar"] is None:
            state["colorbar"] = figure.colorbar(drawn, cax=colorbar_axis)
        else:
            state["colorbar"].update_normal(drawn)
        state["colorbar"].set_label(bar_label, fontsize=8)
        figure.canvas.draw_idle()

    coalescer = gui_core.Coalescer(figure, _draw)

    def write_back():
        """Show the clamped parameters in the widgets, without re-firing their callbacks.

        ``clamped`` may move a different parameter from the one being dragged --- raising
        eps_tide can push lambda up, raising eps_spin with a large tide can move it too ---
        so every slider is refreshed, not just the active one.  ``Slider.set_val`` fires
        ``on_changed``, hence the guard: without it the write-back would recurse.
        """
        state["refreshing"] = True
        try:
            for slider, name, in_degrees in slider_info:
                value = getattr(state["parameters"], name)
                shown = np.degrees(value) if in_degrees else value
                if abs(slider.val - shown) > 1e-9:
                    slider.set_val(shown)
        finally:
            state["refreshing"] = False

    def apply_preset(updates):
        """Set every parameter at once, clamp, and refresh the widgets."""
        state["parameters"] = clamped(replace(state["parameters"], **updates))
        write_back()
        redraw()

    def slider_callback(name, in_degrees):
        def callback(value):
            if state["refreshing"]:
                return
            # The widgets show degrees; Parameters holds radians.  The conversion lives
            # here so nothing outside the GUI has to know about it.
            updated = np.radians(value) if in_degrees else value
            state["parameters"] = clamped(replace(state["parameters"], **{name: updated}))
            preset.manual()             # these are the user's numbers now
            write_back()
            redraw()
        return callback

    def model_callback(label):
        state["model"] = label
        redraw()

    def view_callback(label):
        state["visibility"] = label.startswith("u,v")
        redraw()

    def resolution_callback(label):
        state["resolution"] = RESOLUTION_LABELS.index(label)
        redraw()

    sliders = []
    slider_info = []
    for index, (name, low, high, label, in_degrees, number_format) in enumerate((
            ("tide", 0.0, 0.14, r"$\epsilon_\mathrm{tide}$", False, "%.3f"),
            ("spin", 0.0, 0.09, r"$\epsilon_\mathrm{spin}$", False, "%.3f"),
            ("inclination", 0.0, 90.0, r"$i$  [deg]", True, "%.0f"),
            ("latitude", -90.0, 90.0, r"$\lambda_\mathrm{tide}$  [deg]", True, "%.0f"),
            ("longitude", 0.0, 360.0, r"$\varphi_\mathrm{tide}$  [deg]", True, "%.0f"),
            ("limb", 0.0, 1.0, r"limb darkening", False, "%.2f"))):
        current = getattr(state["parameters"], name)
        sliders.append(gui_core.add_slider(
            figure, index, label, low, high,
            np.degrees(current) if in_degrees else current,
            number_format, slider_callback(name, in_degrees)))
        slider_info.append((sliders[-1], name, in_degrees))

    model_radio = RadioButtons(figure.add_axes([0.06, 0.02, 0.20, 0.21]), MODEL_LABELS)
    model_radio.on_clicked(model_callback)
    view_radio = RadioButtons(figure.add_axes([0.32, 0.06, 0.18, 0.13]),
                              ("image", "u,v"))
    view_radio.on_clicked(view_callback)
    resolution_radio = RadioButtons(figure.add_axes([0.56, 0.02, 0.18, 0.21]),
                                    RESOLUTION_LABELS, active=state["resolution"])
    resolution_radio.on_clicked(resolution_callback)

    preset = gui_core.PresetRadio(figure, [0.76, 0.02, 0.18, 0.21], PRESETS, apply_preset)

    gui_core.retain_widgets(
        figure, (sliders, model_radio, view_radio, resolution_radio, preset.radio))

    for label, order in FIT_ORDERS.items():
        print(f"{label}: Gram conditioning by order --- "
              + ", ".join(f"m={m}: {c:.1f}"
                          for m, c in enumerate(fitted.block_condition(order))))

    coalescer.now()
    return figure, state


def main():  # pragma: no cover - interactive
    from matplotlib import pyplot as plt
    build()
    plt.show()


if __name__ == "__main__":  # pragma: no cover
    main()
