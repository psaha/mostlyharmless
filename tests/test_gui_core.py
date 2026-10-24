"""The shared GUI scaffolding: the flip, the colour scale, and the debounce.

Only the parts that do not need a display are tested here.  That is deliberate: the
duplication this module removes produced two bugs that only appeared in a running window,
so the pieces worth testing are the pure ones, and the matplotlib wiring is kept down to a
handful of one-line calls in ``gui_core``.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import gui_core


def test_the_flip_puts_east_on_the_left():
    """D24: East increases leftward, and imshow draws column 0 on the left."""
    panel = np.arange(6.0).reshape(2, 3)
    flipped = gui_core.flip(panel)
    assert flipped.shape == panel.shape
    assert np.array_equal(flipped[0], panel[0][::-1])


def test_the_colour_limit_is_reused_until_the_key_changes():
    """Switching rows must not rescale the same star; changing a slider must."""
    state = {}
    first = gui_core.colour_limit(("a", 0), state, [np.array([[1.0, -3.0]])])
    assert first == pytest.approx(3.0)
    # a different panel at the same key: the scale does not move
    second = gui_core.colour_limit(("a", 0), state, [np.array([[0.1]])])
    assert second == pytest.approx(3.0)
    # a new key recomputes
    third = gui_core.colour_limit(("a", 1), state, [np.array([[0.25]])])
    assert third == pytest.approx(0.25)


def test_the_colour_limit_ignores_nan_and_survives_an_empty_panel():
    """The rim of a rendered image is NaN by design, and a zero panel must not divide by it."""
    state = {}
    assert gui_core.colour_limit(("k",), state, [np.array([[np.nan, 2.0]])]) == pytest.approx(2.0)
    empty = {}
    assert gui_core.colour_limit(("k",), empty, [np.zeros((3, 3))]) > 0.0


class _Timer:
    def __init__(self):
        self.starts = 0
        self.stops = 0
        self.callback = None

    def add_callback(self, callback):
        self.callback = callback

    def start(self):
        self.starts += 1

    def stop(self):
        self.stops += 1


class _Canvas:
    def __init__(self, timer=None):
        self._timer = timer

    def new_timer(self, interval=200):
        if self._timer is None:
            raise AttributeError("this backend has no timers")
        return self._timer


class _Figure:
    def __init__(self, timer=None):
        self.canvas = _Canvas(timer)


def test_the_coalescer_renders_once_per_burst():
    """Every event marks stale; one timer fire produces one render."""
    calls = []
    timer = _Timer()
    coalescer = gui_core.Coalescer(_Figure(timer), lambda: calls.append(1))

    for _ in range(5):                      # a drag: five events
        coalescer.request()
    assert calls == []                      # nothing rendered yet
    assert timer.starts == 5 and timer.stops == 5

    timer.callback()                        # the motion stops
    assert len(calls) == 1
    timer.callback()                        # a spurious second fire
    assert len(calls) == 1
    print(f"\ncoalesced: {5} events -> {len(calls)} render")


def test_a_backend_without_timers_renders_inline():
    """A non-interactive canvas must still draw, just without the pause."""
    calls = []
    coalescer = gui_core.Coalescer(_Figure(), lambda: calls.append(1))
    coalescer.request()
    assert len(calls) == 1


def test_a_raising_draw_does_not_leave_the_view_stale():
    """The callback must swallow the error, not the state update."""
    def boom():
        raise RuntimeError("render failed")
    coalescer = gui_core.Coalescer(_Figure(), boom)
    coalescer.request()                     # prints a traceback, does not raise
    assert coalescer.stale is False


def test_the_brightness_colormap_runs_black_to_white_and_hides_the_rim():
    """A brightness wants a linear ramp; an unrendered rim must not be a pale halo."""
    cmap = gui_core.image_colormap()
    low = np.asarray(cmap(0.0))[:3]
    high = np.asarray(cmap(1.0))[:3]
    rim = np.asarray(cmap(np.nan))[:3]
    print(f"\ngrey ramp: {low.round(2)} -> {high.round(2)}, rim {rim.round(2)}")
    assert low.max() < 0.1                  # black floor
    assert high.min() > 0.9                 # white peak
    assert rim.max() < 0.1                  # the rim is black, not pale
