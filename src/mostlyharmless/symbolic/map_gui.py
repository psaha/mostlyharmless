"""The map GUI: a starspot in spherical harmonics, its hemispheres and its visibility.

T5.  The manual's picture figures made live --- the same coefficients, the same two views,
the same flips and the same colour-scale rule, driven by sliders instead of hard-coded.
A second module rather than a mode of the deformation GUI, because there is no deformation
here: the model is a plain ``Coeffs``, with no level solve, no projection and no fit, so a
redraw costs a fraction of the other one's.

The map is a polar cap on a uniformly bright surface, carried in the body frame, rotated
into the observer frame by the Wigner matrices.  The spin's position angle is a gauge for a
single spot --- it is degenerate with the rotation phase --- so it is pinned at zero and the
rotation phase is the slider.

The scaffolding deliberately mirrors ``gui``: one colour scale per parameter set, sampled
on a coalescing timer so that a drag renders once rather than once per mouse motion, and
panels flipped so that East increases to the left (D24).  Separating that scaffolding from
this module is the remainder of T5, so that the two GUIs share it rather than drift.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from mostlyharmless.symbolic import gui_core, maps
from mostlyharmless.visylm.visibility import visibility, visibility_uniform_disk
from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.wigner import rotate_coefficients

__all__ = ["MapParameters", "L_MAX_LABELS", "VIEW_LABELS", "RESOLUTIONS",
           "coefficients", "compute_map", "compute_visibility", "build", "main"]

# The truncations offered.  A radio rather than a slider because changing it recomputes the
# coefficients, and because the comparison is between named orders.
L_MAX_LABELS = (4, 6, 8, 16)
VIEW_LABELS = ("map", "u,v")
# The Earth is a different kind of model: three colour channels rather than one brightness,
# so its panels are composited rather than shaded.
MODEL_LABELS = ("spot", "Earth")
RESOLUTION_LABELS = ("coarse", "medium", "fine")
# (pixels, (u,v) grid).  There are no samples here: the map is analytic, so the only
# numerical choices are display resolution and the size of the frequency grid.
RESOLUTIONS = ((121, 23), (241, 31), (361, 41))

Q_MAX = 4.0

# Named parameter sets for the geometry checks, labelled by what you should see rather than
# by their angles.  The first four are the ones that pin the display: a spot whose position
# on the page is stated can only appear there if East is left, North is up and the far
# hemisphere is mirrored correctly.  See the closed entry for GUI verification in the
# roadmap for why these are worth having.
PRESETS = {
    "spot left": dict(half_angle=np.radians(20.0), depth=0.92,
                       colatitude=np.radians(60.0), alpha=np.radians(90.0),
                       beta=0.0, pa=0.0),
    "spot right": dict(half_angle=np.radians(20.0), depth=0.92,
                        colatitude=np.radians(60.0), alpha=np.radians(270.0),
                        beta=0.0, pa=0.0),
    "spot above": dict(half_angle=np.radians(20.0), depth=0.92,
                        colatitude=np.radians(60.0), alpha=0.0, beta=0.0, pa=0.0),
    "spot behind": dict(half_angle=np.radians(20.0), depth=0.92,
                         colatitude=np.radians(60.0), alpha=0.0,
                         beta=np.radians(90.0), pa=0.0),
}


@dataclass(frozen=True)
class MapParameters:
    """The spot and the orientation: six sliders and a radio.

    ``colatitude`` is what makes ``alpha`` mean something.  The cap is built zonal --- at
    the body's north pole --- and an axisymmetric map is invariant under a rotation about
    its own axis, so with the spot at the pole the rotation phase is inert.  Off the pole
    it is the spot's longitude, and the position angle of the sky becomes a genuinely
    independent rotation rather than a gauge.
    """

    half_angle: float = maps.HALF_ANGLE        # radians
    depth: float = maps.SPOT_DEPTH
    colatitude: float = 0.0                       # spot's angle from the spin axis
    alpha: float = 0.0                            # spot's longitude, radians
    beta: float = maps.BETA                    # spin axis to line of sight, radians
    pa: float = 0.0                               # position angle of the sky, radians
    l_max: int = 16


def coefficients(parameters):
    """The spot's coefficients in the observer frame.

    The body-frame cap is zonal, so the rotation is what gives it a longitude: ``alpha``
    turns the spot about the spin, ``beta`` tilts the spin to the line of sight, and the
    position angle stays zero as the gauge it is.
    """
    cap = maps.cap_coefficients(parameters.l_max, parameters.half_angle,
                                    parameters.depth)
    placed = rotate_coefficients(cap, 0.0, parameters.colatitude, 0.0)
    return rotate_coefficients(placed, parameters.alpha, parameters.beta,
                               parameters.pa)


_EARTH = {}


def earth_coefficients():
    """The Earth's three channels, Lanczos-tapered, in the body frame.

    Cached because the taper is a pass over 3 x 1681 coefficients and the renderer asks for
    them on every redraw.  The taper is not optional here: a colour image of the Earth is
    close to a jump at every coastline, and without it the panels show the ringing and not
    the continents.
    """
    if not _EARTH:
        from mostlyharmless.symbolic import earth, earth_coefficients as emitted

        for name in ("RED", "GREEN", "BLUE"):
            blocks = tuple(np.asarray(block, dtype=complex) for block in getattr(emitted, name))
            _EARTH[name] = earth.taper(Coeffs(blocks), l_max=emitted.L_MAX)
    return _EARTH["RED"], _EARTH["GREEN"], _EARTH["BLUE"]


def _earth_panels(parameters, pixels, l_max=None):
    """The Earth's near and far faces as RGB arrays, one colour per channel.

    Each channel is rendered separately and then scaled by its own maximum over both faces,
    so the two panels share one colour balance and the disc is the only thing that differs.
    ``l_max`` truncates the emitted set: three channels times two faces is six renders per
    redraw, and the cost is the mode count.
    """
    from mostlyharmless.symbolic.earth import truncate

    near, far = [], []
    for coeffs in earth_coefficients():
        if l_max is not None:
            coeffs = truncate(coeffs, l_max)
        faced = rotate_coefficients(coeffs, parameters.alpha, parameters.beta, parameters.pa)
        near.append(maps.sky_map(faced, n=pixels))
        far.append(maps.sky_map(faced, back=True, n=pixels))
    panels = []
    for channels in (near, far):
        stack = np.stack(channels, axis=-1)
        peak = float(np.nanmax(stack)) or 1.0
        panels.append(np.clip(stack / peak, 0.0, 1.0))
    return panels


def compute_map(parameters, pixels=241, model="spot"):
    """Near and far hemispheres, in units of the mean, East increasing leftward.

    ``sky_map`` already mirrors the far hemisphere so that it reads as a view from the
    other side.  The display flip is the scaffolding's, applied once by ``build``.

    For ``model == "Earth"`` the panels are RGB arrays rather than brightnesses, so
    ``build`` must draw them with ``imshow`` and no colour map at all.
    """
    if model == "Earth":
        return _earth_panels(parameters, pixels, l_max=parameters.l_max)
    coeffs = coefficients(parameters)
    """Near and far hemispheres, in units of the mean, East increasing leftward.

    ``sky_map`` already mirrors the far hemisphere so that it reads as a view from the
    other side.  The display flip is the scaffolding's, applied once by ``build``.
    """
    coeffs = coefficients(parameters)
    factor = maps.MEAN_FACTOR
    return (factor * maps.sky_map(coeffs, n=pixels),
            factor * maps.sky_map(coeffs, back=True, n=pixels))


def compute_visibility(parameters, grid=23, q_max=Q_MAX, model="spot"):
    """Re and Im of Delta V against the uniform disk, normalised at the origin.

    Subtracting the uniform disk is what makes the deformation's signature visible; the
    map is normalised so that ``V(0,0)`` is the flux, and the map being non-negative makes
    ``|V| <= V(0,0)`` hold, which is the bound the manual's figure uses.
    """
    if model == "Earth":
        # V is linear in the map, so the three channels add before the transform: the
        # visibility of the colour image is the sum of the channel visibilities.
        from mostlyharmless.symbolic.earth import truncate

        channels = [truncate(channel, parameters.l_max)
                    if parameters.l_max else channel
                    for channel in earth_coefficients()]
        total = Coeffs(tuple(sum(block for block in blocks)
                             for blocks in zip(*(channel.blocks for channel in channels))))
        coeffs = rotate_coefficients(total, parameters.alpha, parameters.beta, parameters.pa)
    else:
        coeffs = coefficients(parameters)
    flux = float(np.real(visibility(coeffs, 0.0, 0.0)))

    axis = np.linspace(-q_max, q_max, grid)
    u_grid, v_grid = np.meshgrid(axis, axis, indexing="ij")
    reference = visibility_uniform_disk(u_grid, v_grid, 1.0)
    values = np.array([[visibility(coeffs, u, v) for v in axis] for u in axis])
    return values / flux - reference


def build(parameters=None):  # pragma: no cover - needs a display
    """The widget window; returns ``(figure, state)`` with the live parameters in it."""
    from matplotlib import pyplot as plt
    from matplotlib.widgets import RadioButtons

    state = {"parameters": parameters or MapParameters(),
             "model": "spot", "view": "map", "resolution": 1, "limit_key": None,
             "limit": None, "colorbar": None, "refreshing": False}

    figure = plt.figure(figsize=(9.2, 6.6))
    axes, colorbar_axis = gui_core.panel_axes(figure)

    def redraw(*_):
        coalescer.request()

    def _draw():
        pixels, grid = RESOLUTIONS[state["resolution"]]
        current = state["parameters"]
        model = state["model"]
        is_earth = model == "Earth" and state["view"] == "map"

        if state["view"] == "u,v":
            values = compute_visibility(current, grid=grid, model=model)
            panels = (values.real, values.imag)
            titles = (r"$\mathrm{Re}\,\Delta V$", r"$\mathrm{Im}\,\Delta V$")
            extent = (Q_MAX, -Q_MAX, -Q_MAX, Q_MAX)
            labels = (r"$v$ (East)", r"$u$ (North)")
            bar_label = r"$\Delta V$"
        elif is_earth:
            # RGB: three rendered channels, no colour map and no colourbar
            panels = compute_map(current, pixels=pixels, model=model)
            titles = ("Earth, near side", "Earth, far side")
            extent = (1.0, -1.0, -1.0, 1.0)
            labels = (r"$y$ (East)", r"$x$ (North)")
            bar_label = None
        else:
            panels = compute_map(current, pixels=pixels, model=model)
            titles = (r"$I/\langle I\rangle$, near", r"$I/\langle I\rangle$, far")
            extent = (1.0, -1.0, -1.0, 1.0)
            labels = (r"$y$ (East)", r"$x$ (North)")
            bar_label = r"$I/\langle I\rangle$"

        # the flip is the scaffolding's, not this model's
        panels = tuple(gui_core.flip(panel) for panel in panels)

        if is_earth:
            for axis_handle, data, title in zip(axes, panels, titles):
                axis_handle.clear()
                gui_core.darken(axis_handle)
                axis_handle.imshow(np.clip(data, 0.0, 1.0), origin="lower",
                                   extent=extent, interpolation="nearest")
                axis_handle.set_title(title, fontsize=9)
                axis_handle.set_xlabel(labels[0], fontsize=8)
                axis_handle.set_ylabel(labels[1], fontsize=8)
            if state["colorbar"] is not None:
                state["colorbar"].ax.set_visible(False)     # a colour image needs no bar
            figure.canvas.draw_idle()
            return

        limit = gui_core.colour_limit(
            (current, state["view"], state["resolution"]), state, panels)

        # a brightness wants a linear grey ramp; a signed Delta V wants a diverging map
        brightness = state["view"] == "map"
        cmap = gui_core.image_colormap() if brightness else "RdBu_r"

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

        if state["colorbar"] is None:
            state["colorbar"] = figure.colorbar(drawn, cax=colorbar_axis)
        else:
            state["colorbar"].update_normal(drawn)
        state["colorbar"].ax.set_visible(True)
        state["colorbar"].set_label(bar_label, fontsize=8)
        figure.canvas.draw_idle()

    coalescer = gui_core.Coalescer(figure, _draw)

    def write_back():
        """Show the current parameters in the widgets, without re-firing their callbacks.

        Needed because a preset sets them all at once, and because ``Slider.set_val`` fires
        ``on_changed``: without the guard the write-back would recurse through the slider
        callback and then through ``preset.manual``.
        """
        state["refreshing"] = True
        try:
            for slider, name, in_degrees in slider_info:
                shown = getattr(state["parameters"], name)
                shown = np.degrees(shown) if in_degrees else shown
                if abs(slider.val - shown) > 1e-9:
                    slider.set_val(shown)
        finally:
            state["refreshing"] = False

    def slider_callback(name, in_degrees):
        def callback(value):
            if state["refreshing"]:
                return
            converted = np.radians(value) if in_degrees else value
            state["parameters"] = replace(state["parameters"], **{name: converted})
            preset.manual()             # these are the user's numbers now
            redraw()
        return callback

    def apply_preset(updates):
        state["parameters"] = replace(state["parameters"], **updates)
        write_back()
        redraw()

    sliders = []
    slider_info = []
    for index, (name, low, high, label, in_degrees, number_format) in enumerate((
            ("half_angle", 1.0, 90.0, r"half-angle  [deg]", True, "%.0f"),
            ("depth", 0.0, 0.92, r"depth", False, "%.2f"),
            ("colatitude", 0.0, 90.0, r"$\theta_\mathrm{spot}$  [deg]", True, "%.0f"),
            ("alpha", 0.0, 360.0, r"$\alpha_\mathrm{spot}$  [deg]", True, "%.0f"),
            ("beta", 0.0, 90.0, r"$\beta$  [deg]", True, "%.0f"),
            ("pa", 0.0, 90.0, r"PA  [deg]", True, "%.0f"))):
        current = getattr(state["parameters"], name)
        sliders.append(gui_core.add_slider(
            figure, index, label, low, high,
            np.degrees(current) if in_degrees else current,
            number_format, slider_callback(name, in_degrees)))
        slider_info.append((sliders[-1], name, in_degrees))

    def l_max_callback(label):
        state["parameters"] = replace(state["parameters"], l_max=int(label.strip("l =")))
        redraw()

    def view_callback(label):
        state["view"] = label
        redraw()

    def resolution_callback(label):
        state["resolution"] = RESOLUTION_LABELS.index(label)
        redraw()

    # the radio must start on the degree the parameters actually carry, or it shows l = 4
    # while the map is built at 16 --- a mismatch nothing in the picture reveals
    l_max_radio = RadioButtons(
        figure.add_axes([0.22, 0.02, 0.16, 0.21]),
        tuple(f"l = {value}" for value in L_MAX_LABELS),
        active=L_MAX_LABELS.index(state["parameters"].l_max))
    l_max_radio.on_clicked(l_max_callback)
    view_radio = RadioButtons(figure.add_axes([0.40, 0.06, 0.16, 0.13]), VIEW_LABELS)
    view_radio.on_clicked(view_callback)
    resolution_radio = RadioButtons(figure.add_axes([0.58, 0.02, 0.16, 0.21]),
                                    RESOLUTION_LABELS, active=state["resolution"])
    resolution_radio.on_clicked(resolution_callback)

    def model_callback(label):
        state["model"] = label
        redraw()

    model_radio = RadioButtons(figure.add_axes([0.04, 0.02, 0.16, 0.13]), MODEL_LABELS)
    model_radio.on_clicked(model_callback)

    preset = gui_core.PresetRadio(figure, [0.76, 0.02, 0.18, 0.21], PRESETS, apply_preset)

    gui_core.retain_widgets(
        figure, (sliders, model_radio, l_max_radio, view_radio, resolution_radio,
                 preset.radio))

    coalescer.now()
    return figure, state


def main():  # pragma: no cover - interactive
    from matplotlib import pyplot as plt
    build()
    plt.show()


if __name__ == "__main__":  # pragma: no cover
    main()
