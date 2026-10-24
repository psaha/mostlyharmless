"""Shared scaffolding for the two GUIs: coalescing, colour scale, layout, flips.

Both ``gui`` and ``map_gui`` need the same five things around their own model code, and
each was written once in each file.  That is the duplication T5 exists to remove: two
copies of a debounce timer, a colour-scale rule and a column flip are exactly how the
forward-reference and lost-widget bugs happened the first time.

Only the pure parts are tested --- ``flip`` and ``colour_limit`` --- because ``build()``
needs a display.  Keeping those two pure, and the matplotlib wiring down to a handful of
one-line calls here, is the reason this module exists rather than a copy in each GUI.
"""

from __future__ import annotations

import numpy as np

__all__ = ["flip", "colour_limit", "Coalescer", "PresetRadio", "PANEL_RECTS",
           "COLOURBAR_RECT", "panel_axes", "add_slider", "retain_widgets"]

# The panels and the colourbar sit in one band, with the radios below and the slider stack
# above.  A radio row placed at the panels' top overlaps their titles, which are set just
# above the axes and are not part of the rectangle.
PANEL_RECTS = ([0.06, 0.26, 0.36, 0.36], [0.50, 0.26, 0.36, 0.36])
COLOURBAR_RECT = (0.90, 0.26, 0.015, 0.36)
SLIDER_RECT = (0.20, 0.94, 0.64, 0.026)          # x, top, width, height
SLIDER_SPACING = 0.043


def flip(panel):
    """East increases to the left (D24); ``imshow`` draws column 0 on the left.

    Flipping once, here, rather than building the arrays mirrored, is what keeps the model
    code free of a display convention.
    """
    return np.asarray(panel)[:, ::-1]


def colour_limit(key, state, panels):
    """One colour scale per parameter set, shared by every row.

    Recomputing it from whichever row is showing makes the same star change colour when
    the model changes, which reads as a difference in the model rather than in the display.
    ``key`` must be everything *except* the row: parameters, view and resolution.
    """
    if state.get("limit_key") != key:
        stacked = np.concatenate([np.asarray(panel).ravel() for panel in panels])
        state["limit_key"] = key
        state["limit"] = float(np.nanmax(np.abs(stacked))) or 1e-30
    return state["limit"]


class Coalescer:
    """One render per pause in the motion, not one per event.

    A slider fires its callback on every mouse motion while it is dragged, and each change
    costs a full render.  Rendering inside the callback enqueues work faster than the event
    loop can drain it, and the window stops responding --- which looks like a hang but is a
    backlog.  So: mark stale, restart a short timer, render once when the motion stops.
    Without a timer (a non-interactive canvas) it renders at once.

    The timer must be created after the draw callable exists, since naming a local before
    its ``def`` runs is an ``UnboundLocalError`` and not a forward reference.  That is why
    this is a class: construction takes the callable, so the ordering cannot be got wrong.
    """

    def __init__(self, figure, draw, interval=200):
        self.draw = draw
        self.stale = False
        self.timer = None
        try:
            timer = figure.canvas.new_timer(interval=interval)
        except Exception:                       # a backend without timers
            timer = None
        if timer is not None:
            timer.add_callback(self.fire)
            self.timer = timer

    def fire(self, *_):
        """Render if stale; a callback that raises must not vanish."""
        if not self.stale:
            return
        self.stale = False
        try:
            self.draw()
        except Exception:                                   # pragma: no cover
            import traceback
            traceback.print_exc()

    def request(self, *_):
        """Mark the view stale and schedule one render."""
        self.stale = True
        if self.timer is None:
            self.fire()
        else:
            self.timer.stop()
            self.timer.start()

    def now(self):
        """Paint once immediately, without waiting for the timer."""
        self.stale = True
        self.fire()


def panel_axes(figure):
    """The two panel axes and the shared colourbar axis, in the standard layout."""
    axes = (figure.add_axes(PANEL_RECTS[0]), figure.add_axes(PANEL_RECTS[1]))
    return axes, figure.add_axes(COLOURBAR_RECT)


def add_slider(figure, index, label, low, high, value, number_format, callback,
               spacing=SLIDER_SPACING):
    """One slider in the stack above the panels, positioned by ``index``."""
    from matplotlib.widgets import Slider

    x, top, width, height = SLIDER_RECT
    handle = figure.add_axes([x, top - spacing * index, width, height])
    slider = Slider(handle, label, low, high, valinit=value, valfmt=number_format)
    slider.label.set_fontsize(9)
    slider.valtext.set_fontsize(8)
    slider.on_changed(callback)
    return slider


def retain_widgets(figure, widgets):
    """A widget with no reference is collected and its callbacks stop firing.

    That looks like a frozen control, and it is the second bug this scaffolding exists to
    prevent, so the reference lives on the figure rather than on a local.
    """
    figure._widgets = widgets


class PresetRadio:
    """A radio group of named parameter sets, for reproducing a geometry check.

    The checks that pin the display --- East on the left, North up, the far hemisphere ---
    are not expressed in any test, because they are statements about a picture.  Making them
    a radio entry rather than a paragraph of instructions means the next reader can see the
    geometry without being told which sliders to move, and cannot get the expected answer
    wrong: the labels say what you should see, not what the angles are.

    A radio group has no "off", so there is a ``manual`` entry that means "these are the
    user's own numbers".  It is selected automatically whenever a slider moves, and
    selecting it does nothing --- which is why the callback is safe to fire on a programmatic
    ``set_active``.
    """

    MANUAL = "manual"

    def __init__(self, figure, rect, presets, apply):
        from matplotlib.widgets import RadioButtons

        self.presets = dict(presets)
        self.apply = apply
        self.selected = self.MANUAL
        self.radio = RadioButtons(figure.add_axes(rect), (self.MANUAL,) + tuple(self.presets))
        self.radio.on_clicked(self.clicked)

    def clicked(self, label):
        """Apply a preset; ``manual`` is a no-op."""
        if label == self.MANUAL:
            return
        self.selected = label
        self.apply(self.presets[label])

    def manual(self):
        """The user moved a slider: the parameters are their own again."""
        if self.selected != self.MANUAL:
            self.selected = self.MANUAL
            self.radio.set_active(0)
