"""The interactive front end: sliders and toggles over the Stage 6 model.

Stage 8, D25.  Matplotlib widgets only --- no Jupyter, and nothing added to the library.

Two radio buttons, matching the four rows of the Stage 6 figures: the *view* is a rendering
or a pair of ``(u,v)`` maps, and the *model* is the affine (analytic) one or the actual
(numeric) surface.  The sliders cover the six parameters.  The brightness is fixed at
von Zeipel, ``I ∝ g`` (D21), and the position angle is pinned to zero as a gauge, so North
is vertical.

The closure bound is explicit.  Increasing the tide's latitude *loosens* it --- at the
equator the centrifugal quadrupole deepens the potential by 1.5 eps_cent, while along the
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

from mostlyharmless.symbolic import affine, render
from mostlyharmless.visylm.deformation import Deformation
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk

__all__ = ["Parameters", "closure_depth", "largest_spin", "smallest_latitude",
           "clamped", "compute_image", "compute_visibility", "build", "main"]

# D25: the analytic row keeps l <= 2, the eps^2 truncation of the l >= 4 content with the
# l = 2 coefficient resummed.  Not a choice to make.
SPHERE_L_MAX = 2
CLOSURE = 0.148
SAMPLES = 8000
PIXELS = 120
GRID = 15
Q_MAX = 2.0


@dataclass(frozen=True)
class Parameters:
    """The five free parameters.

    The spin phase is deliberately absent: this front end is for a spin and tide alone, and
    with no starspots the phase is degenerate with the tide's longitude --- only their
    difference is observable.
    """

    tide: float = 0.0
    spin: float = 0.05
    inclination: float = np.pi / 4.0
    latitude: float = 0.0
    longitude: float = 0.0

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
    """The largest ``eps_cent`` that still closes, by root solve.

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


def compute_image(parameters, analytic, samples=SAMPLES, pixels=PIXELS):
    """The near and far hemispheres, as fractional departures from their own mean.

    The affine image is the sphere's image stretched in the sky plane --- exact as a
    coordinate stretch, no per-point sum --- and the numeric one samples the level
    surface.  Both are mean-normalised so the flux is out of the way (D22).
    """
    field = parameters.field()
    eye, (east, north) = _view(parameters.inclination)
    stretch = affine.projected_stretch(field, eye)
    scale = abs(float(np.linalg.det(stretch)))

    def contrast(image):
        return image / np.nanmean(image) - 1.0

    panels = []
    for sense in (1.0, -1.0):
        direction = sense * eye
        if analytic:
            sphere = render.sphere_samples(field, direction, n_samples=samples)
            sky = np.stack([sphere[0] @ east, sphere[0] @ north], axis=-1)
            moved = sky @ stretch.T
            points = moved[:, :1] * east + moved[:, 1:] * north
            image = render.sky_image(points, sphere[1] / scale, sphere[2] * scale,
                                     east, north, half_width=1.5, n=pixels)
        else:
            star = render.surface_samples(field, direction, n_samples=samples)
            image = render.sky_image(*star, east, north, half_width=1.5, n=pixels)
        panels.append(contrast(image))
    return panels


def compute_visibility(parameters, analytic, samples=SAMPLES, grid=GRID, q_max=Q_MAX):
    """Re and Im of Delta V against the unit sphere, for one model.

    Zero at the origin by construction, so the lobes the models share are gone and only
    the deformation's signature remains.  ``u`` is North, so it is the first index.
    """
    field = parameters.field()
    eye, (east, north) = _view(parameters.inclination)

    if analytic:
        coefficients = affine.sphere_coefficients(field, l_max=SPHERE_L_MAX,
                                                 n_samples=samples,
                                                 inclination=parameters.inclination)
        sphere_flux = complex(visibility(coefficients, 0.0, 0.0))
        stretch = affine.projected_stretch(field, eye)

        def sphere_sum(u, v):
            return visibility(coefficients, u, v) / sphere_flux

        def model(u, v):
            return affine.affine_visibility(sphere_sum, stretch, u, v)
    else:
        star = render.surface_samples(field, eye, n_samples=samples)
        flux = render.direct_visibility(*star, east, north, 0.0, 0.0)

        def model(u, v):
            return (render.direct_visibility(*star, east, north, u, v) / flux)

    axis = np.linspace(-q_max, q_max, grid)
    u_grid, v_grid = np.meshgrid(axis, axis, indexing="ij")
    reference = visibility_uniform_disk(u_grid, v_grid, 1.0)
    return np.array([[model(u, v) for v in axis] for u in axis]) - reference


def build(parameters=None):  # pragma: no cover - needs a display
    """The widget window; returns (figure, state) with the live parameters in ``state``."""
    from matplotlib import pyplot as plt
    from matplotlib.widgets import RadioButtons, Slider

    state = {"parameters": clamped(parameters or Parameters()),
             "analytic": False, "visibility": False}

    figure = plt.figure(figsize=(9.2, 6.6))
    axes = (figure.add_axes([0.06, 0.10, 0.36, 0.36]),
            figure.add_axes([0.50, 0.10, 0.36, 0.36]))
    # a fixed slot for the shared colourbar, clear of the radio buttons
    colorbar_axis = figure.add_axes([0.90, 0.10, 0.015, 0.36])
    state["colorbar"] = None

    def redraw(*_):
        # A widget callback that raises would otherwise fail silently.
        try:
            _redraw()
        except Exception:                                   # pragma: no cover
            import traceback
            traceback.print_exc()

    def _redraw():
        current = state["parameters"]
        if state["visibility"]:
            values = compute_visibility(current, state["analytic"])
            panels = (values.real, values.imag)
            titles = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$")
            limit = float(np.max(np.abs(values))) or 1e-30
            extent = (Q_MAX, -Q_MAX, -Q_MAX, Q_MAX)
            labels = (r"$v$ (East)", r"$u$ (North)")
            bar_label = r"$\Delta V$"
        else:
            panels = compute_image(current, state["analytic"])
            titles = (("affine, near", "affine, far") if state["analytic"]
                      else ("the star, near", "the star, far"))
            limit = float(np.nanmax(np.abs(np.concatenate(
                [p.ravel() for p in panels])))) or 1e-30
            extent = (1.5, -1.5, -1.5, 1.5)
            labels = (r"$y$ (East)", r"$x$ (North)")
            bar_label = r"$g/\langle g\rangle - 1$"
        for axis_handle, data, title in zip(axes, panels, titles):
            axis_handle.clear()
            drawn = axis_handle.imshow(data, origin="lower", cmap="RdBu_r",
                                       vmin=-limit, vmax=limit, extent=extent)
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

    def slider_callback(name):
        def callback(value):
            state["parameters"] = clamped(replace(state["parameters"], **{name: value}))
            redraw()
        return callback

    def model_callback(label):
        state["analytic"] = label.startswith("analytic")
        redraw()

    def view_callback(label):
        state["visibility"] = label.startswith("u,v")
        redraw()

    # Matplotlib widgets must stay in scope: a widget with no reference is collected and
    # its callbacks stop firing, which looks like a frozen control.  Keep them on the
    # figure, where they live as long as the window does.
    sliders = []
    for index, (name, low, high, label) in enumerate((
            ("tide", 0.0, 0.14, r"$\epsilon_\mathrm{tide}$"),
            ("spin", 0.0, 0.09, r"$\epsilon_\mathrm{cent}$"),
            ("inclination", 0.0, np.pi / 2.0, r"$i$"),
            ("latitude", -np.pi / 2.0, np.pi / 2.0, r"$\lambda$"),
            ("longitude", 0.0, 2.0 * np.pi, r"$\varphi$"))):
        axis_handle = figure.add_axes([0.06, 0.93 - 0.05 * index, 0.80, 0.026])
        slider = Slider(axis_handle, label, low, high,
                        valinit=getattr(state["parameters"], name))
        slider.on_changed(slider_callback(name))
        sliders.append(slider)

    model_radio = RadioButtons(figure.add_axes([0.90, 0.78, 0.09, 0.10]),
                               ("numeric", "analytic"))
    model_radio.on_clicked(model_callback)
    view_radio = RadioButtons(figure.add_axes([0.90, 0.60, 0.09, 0.10]),
                              ("image", "u,v"))
    view_radio.on_clicked(view_callback)
    figure._widgets = (sliders, model_radio, view_radio)

    redraw()
    return figure, state


def main():  # pragma: no cover - interactive
    from matplotlib import pyplot as plt
    build()
    plt.show()


if __name__ == "__main__":  # pragma: no cover
    main()
