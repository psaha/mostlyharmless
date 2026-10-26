"""A fork of the map GUI for the two-level Earth: one brightness, no colour.

The parent module shows a polar spot in grey, or the Earth in three colour channels.  This one
shows a single map: sea and land at two fixed brightnesses, chosen so that the sphere mean is
one (``make earth-levels``).  There is no spot, no colour compositing, no colourbar special
case and no presets --- the panels are an ordinary brightness image throughout, which is what
makes it the simpler thing to read.

Sea sits at 0.87 and land at 1.31 for a land fraction of 0.29: a fifty per cent jump, which
the eye resolves easily while keeping the step small enough that the coastline ringing stays
modest under the taper the coefficients already carry.
"""

from __future__ import annotations

import numpy as np

from mostlyharmless.symbolic import gui_core, maps
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk
from mostlyharmless.visylm.wigner import rotate_coefficients

__all__ = ["Parameters", "DEGREES", "RESOLUTIONS", "Q_MAX", "VIEW_LABELS",
           "levels_coefficients", "compute_map", "compute_visibility", "build"]

DEGREES = (4, 8, 16, 32)
RESOLUTIONS = ((121, 23), (241, 31), (361, 41))
Q_MAX = 4.0
VIEW_LABELS = ("map", "u,v")
RESOLUTION_LABELS = ("coarse", "medium", "fine")
DEFAULT_L_MAX = 16


class Parameters:
    """The rotation and the cache key: the map itself is fixed."""

    def __init__(self, alpha=0.0, beta=np.radians(60.0), pa=0.0, l_max=DEFAULT_L_MAX):
        self.alpha = alpha
        self.beta = beta
        self.pa = pa
        self.l_max = l_max

    def as_tuple(self):
        return (self.alpha, self.beta, self.pa, self.l_max)


_LEVELS = {}


def levels_coefficients():
    """The committed two-level map, already tapered, in the body frame."""
    if not _LEVELS:
        from mostlyharmless.symbolic import earth_levels as emitted

        blocks = tuple(np.asarray(block, dtype=complex) for block in emitted.LEVELS)
        _LEVELS["coefficients"] = Coeffs(blocks)
        _LEVELS["fraction"] = emitted.FRACTION
        _LEVELS["sea"], _LEVELS["land"] = emitted.SEA, emitted.LAND
        _LEVELS["l_max"] = emitted.L_MAX
    return _LEVELS["coefficients"]


def _faced(parameters):
    """The coefficients as the observer sees them, truncated to the requested degree."""
    from mostlyharmless.symbolic.earth import truncate

    coefficients = truncate(levels_coefficients(), parameters.l_max)
    return rotate_coefficients(coefficients, parameters.alpha, parameters.beta, parameters.pa)


def compute_map(parameters, pixels=241):
    """Near and far hemispheres in units of the sphere mean, East to the left.

    One coefficient set, so one ``sky_map_batch`` call per face --- the ordinary brightness
    path, with no compositing anywhere.
    """
    faced = _faced(parameters)
    near = maps.sky_map(faced, n=pixels)
    far = maps.sky_map(faced, back=True, n=pixels)
    return near, far


def compute_visibility(parameters, grid=23, q_max=Q_MAX):
    """Re and Im of Delta V against the uniform disk, on a grid of baselines."""
    coefficients = _faced(parameters)
    flux = float(np.real(visibility(coefficients, 0.0, 0.0)))
    axis = np.linspace(-q_max, q_max, grid)
    u_grid, v_grid = np.meshgrid(axis, axis, indexing="ij")
    reference = visibility_uniform_disk(u_grid, v_grid, 1.0)
    return visibility(coefficients, u_grid, v_grid) / flux - reference


def build(parameters=None):  # pragma: no cover - needs a display
    """The widget window; returns ``(figure, state)``."""
    from matplotlib import pyplot as plt
    from matplotlib.widgets import RadioButtons

    state = {"parameters": parameters or Parameters(),
             "view": "map", "resolution": 1, "refreshing": False}

    figure = plt.figure(figsize=(9.2, 6.6))

    def _draw():
        pixels, grid = RESOLUTIONS[state["resolution"]]
        current = state["parameters"]
        if state["view"] == "u,v":
            values = compute_visibility(current, grid=grid)
            panels = (values.real, values.imag)
            titles = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$")
            extent = (Q_MAX, -Q_MAX, -Q_MAX, Q_MAX)
            labels = (r"$v$ (East)", r"$u$ (North)")
            bar_label = r"$\Delta V$"
        else:
            panels = compute_map(current, pixels=pixels)
            titles = ("two-level Earth, near", "two-level Earth, far")
            extent = (1.0, -1.0, -1.0, 1.0)
            labels = (r"$y$ (East)", r"$x$ (North)")
            bar_label = r"$I/\langle I\rangle$"

        panels = tuple(gui_core.flip(panel) for panel in panels)
        limit = gui_core.colour_limit(
            (current.as_tuple(), state["view"], state["resolution"]), state, panels)

        brightness = state["view"] == "map"
        cmap = gui_core.image_colormap() if brightness else "RdBu_r"

        drawn = None
        for axis_handle, data, title in zip(axes, panels, titles):
            axis_handle.clear()
            if brightness:
                gui_core.darken(axis_handle)
            drawn = axis_handle.imshow(
                data, origin="lower", cmap=cmap,
                vmin=0.0 if brightness else -limit, vmax=limit, extent=extent)
            axis_handle.set_title(title, fontsize=9)
            axis_handle.set_xlabel(labels[0], fontsize=8)
            axis_handle.set_ylabel(labels[1], fontsize=8)

        if state.get("colorbar") is None:
            state["colorbar"] = figure.colorbar(drawn, cax=colorbar_axis)
        else:
            state["colorbar"].update_normal(drawn)
        state["colorbar"].set_label(bar_label, fontsize=8)
        figure.canvas.draw_idle()

    coalescer = gui_core.Coalescer(figure, _draw)

    def write_back():
        """Show the parameters in the widgets without re-firing their callbacks."""
        state["refreshing"] = True
        current = state["parameters"]
        for slider, value in zip(sliders, (current.alpha, current.beta, current.pa)):
            slider.set_val(value)
        degree_radio.set_active(DEGREES.index(current.l_max))
        view_radio.set_active(VIEW_LABELS.index(state["view"]))
        resolution_radio.set_active(state["resolution"])
        state["refreshing"] = False

    def redraw():
        write_back()
        coalescer.request()

    def changed(name):
        def callback(value):
            if state["refreshing"]:
                return
            setattr(state["parameters"], name, value)
            coalescer.request()
        return callback

    axes, colorbar_axis = gui_core.panel_axes(figure)
    sliders = [
        gui_core.add_slider(figure, 2, "alpha", -np.pi, np.pi,
                            state["parameters"].alpha, "%.2f", changed("alpha")),
        gui_core.add_slider(figure, 1, "beta", 0.0, np.pi,
                            state["parameters"].beta, "%.2f", changed("beta")),
        gui_core.add_slider(figure, 0, "pa", -np.pi, np.pi,
                            state["parameters"].pa, "%.2f", changed("pa")),
    ]

    def degree_callback(label):
        if state["refreshing"]:
            return
        state["parameters"].l_max = int(label.split("=")[-1].strip())
        coalescer.request()

    degree_radio = RadioButtons(
        figure.add_axes([0.22, 0.02, 0.16, 0.21]),
        tuple(f"l = {value}" for value in DEGREES),
        active=DEGREES.index(state["parameters"].l_max))
    degree_radio.on_clicked(degree_callback)

    def view_callback(label):
        if state["refreshing"]:
            return
        state["view"] = label
        coalescer.request()

    view_radio = RadioButtons(figure.add_axes([0.44, 0.06, 0.16, 0.13]), VIEW_LABELS)
    view_radio.on_clicked(view_callback)

    def resolution_callback(label):
        if state["refreshing"]:
            return
        state["resolution"] = RESOLUTION_LABELS.index(label)
        coalescer.request()

    resolution_radio = RadioButtons(
        figure.add_axes([0.62, 0.02, 0.16, 0.21]), RESOLUTION_LABELS, active=1)
    resolution_radio.on_clicked(resolution_callback)

    gui_core.retain_widgets(
        figure, (sliders, degree_radio, view_radio, resolution_radio))
    coalescer.request()
    return figure, state


def main(argv=None):  # pragma: no cover - needs a display
    from matplotlib import pyplot as plt

    build()
    plt.show()


if __name__ == "__main__":  # pragma: no cover
    main()
