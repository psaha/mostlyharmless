"""
Spot test: where does a body-fixed spot appear in the observer's view?

A throwaway exploration script for pinning down the rotation conventions of the
frame change from the star's body frame to the observer's line-of-sight frame.

Observer frame (u, v, w) -- the standard radio-interferometry convention
(Measurement Set 2.0 / CASA; Thompson, Moran & Swenson):
    u-hat points East, v-hat points North, and w-hat is the line of sight,
    pointing from the observer towards the star, so the near hemisphere is
    w > 0 and the visible disk is the unit circle.  The direction cosines l and
    m are parallel to u and v, so l increases with Right Ascension (East) and m
    with Declination (North).  Consequently astronomical images put East to the
    LEFT: the u axis of this plot is inverted so the figure shows the sky as
    seen rather than a map.

Body frame (x, y, z):
    z-hat is the spin axis.  The spot is fixed in this frame.

Rotation
--------
Active, right-handed rotations, composed as

    v_observer = R_z(90 deg - PA) @ R_y(beta) @ R_z(alpha) @ v_body

so that

    alpha  spin phase: rotation of the star about its own axis,
           (the projected spin axis does not depend on alpha)
    beta   inclination: angle between the spin axis and the line of sight,
           0 = pole-on, 90 = equator-on,
    PA     position angle of the sky-projected spin axis, measured from North
           through East (the astronomical convention), in degrees.

The spin axis expressed in observer coordinates is then

    s = (sin(beta) sin(PA), sin(beta) cos(PA), cos(beta))

i.e. (East, North, line of sight).  Its sky projection therefore has length
sin(beta) -- it vanishes when the star is seen pole-on, as it must -- and points
at position angle PA, which on this plot runs counter-clockwise from North.

Because the observer is in the far field the sky is flat, and a surface point at
observer-frame direction v appears on the unit disk at

    East = v[0],   North = v[1]

-- an orthographic projection -- and is on the near hemisphere iff v[2] > 0.

Everything here is exact geometry: no spherical harmonics, no band limit.

Run it, then drag the sliders.  The dashed curve is the locus of the spot over
one full turn; it depends on beta and PA but not on alpha.  Watch for:

  * the spot running around the dashed locus as alpha advances,
  * disappearing (hollow marker, grey) when it goes round the far side,
  * the red arrow holding position angle PA, and shrinking to nothing pole-on.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider

# ----------------------------------------------------------------- fixed inputs
SPOT_THETA_S_DEG = 70.0  # spot colatitude in the body frame, from the spin axis
SPOT_PHI_S_DEG = 150.0  # spot azimuth in the body frame
TRACK_SAMPLES = 361  # points used for the locus of the spot over one turn


# ------------------------------------------------------------------- geometry
def rot_z(angle):
    """Active right-handed rotation of a vector by `angle` about +z."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_y(angle):
    """Active right-handed rotation of a vector by `angle` about +y."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def body_to_observer(alpha, beta, pa):
    """Active rotation taking body-frame components to observer-frame components.

    The observer axes are (East, North, line of sight); `beta` is the inclination
    from the line of sight and `pa` the position angle of the projected spin axis
    measured from North through East.
    """
    return rot_z(0.5 * np.pi - pa) @ rot_y(beta) @ rot_z(alpha)


def unit_vector(theta, phi):
    """Unit vector with polar angle theta and azimuth phi."""
    return np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta),
        ]
    )


def spot_observer_direction(alpha, beta, pa, theta_s, phi_s):
    """Observer-frame direction of the body-fixed spot, as (East, North, w)."""
    return body_to_observer(alpha, beta, pa) @ unit_vector(theta_s, phi_s)


def spin_axis_observer_direction(beta, pa):
    """Observer-frame direction of the spin axis (independent of alpha).

    Written out rather than obtained from the rotation matrix, so it doubles as an
    independent check of `body_to_observer`.
    """
    return np.array(
        [
            np.sin(beta) * np.sin(pa),
            np.sin(beta) * np.cos(pa),
            np.cos(beta),
        ]
    )


def observer_angles(v):
    """(theta, phi) of a unit vector in the observer frame, in radians.

    `phi` is measured from East towards North, the usual mathematical azimuth for
    these axes; use `position_angle` for the astronomical quantity.
    """
    theta = np.arccos(np.clip(v[2], -1.0, 1.0))
    phi = np.arctan2(v[1], v[0])
    return theta, phi


def position_angle(east, north):
    """Position angle in radians, measured from North through East."""
    return np.arctan2(east, north)


def spot_track(beta, pa, theta_s, phi_s, n=TRACK_SAMPLES):
    """Projected (East, North) locus of the spot over one full turn in alpha."""
    alphas = np.linspace(0.0, 2.0 * np.pi, n)
    return np.array(
        [spot_observer_direction(a, beta, pa, theta_s, phi_s)[:2] for a in alphas]
    )


# ---------------------------------------------------------------------- figure
def draw(ax, alpha, beta, pa, theta_s, phi_s):
    """Draw the observer's view for one set of angles.  Angles in radians."""
    ax.clear()
    ax.set_aspect("equal")
    ax.set_xlim(-1.18, 1.18)
    ax.set_ylim(-1.18, 1.18)
    ax.invert_xaxis()  # u is East, and East is to the left on the sky
    ax.set_xlabel("u   (East)          East is left, North is up")
    ax.set_ylabel("v   (North)", rotation=0, labelpad=10)
    ax.set_title("Spot test -- observer's view (far field)", fontsize=11)

    # Visible hemisphere: the unit disk.  Outside the circle is sky.
    ax.add_patch(
        Circle((0.0, 0.0), 1.0, facecolor="0.94", edgecolor="0.35", lw=1.2, zorder=0)
    )
    ax.axhline(0.0, color="0.85", lw=0.6, zorder=1)
    ax.axvline(0.0, color="0.85", lw=0.6, zorder=1)

    # Compass letters, so the sky orientation is unmistakable.  They sit just
    # inside the limb, clear of the axis labels.  Note that "E" is on the left
    # because the u axis is inverted.
    for label, (cu, cv) in (
        ("E", (0.93, 0.0)),
        ("W", (-0.93, 0.0)),
        ("N", (0.0, 0.93)),
        ("S", (0.0, -0.93)),
    ):
        ax.text(cu, cv, label, ha="center", va="center", fontsize=9, color="0.3")

    # Locus of the spot over one full turn (depends on beta and PA only).
    track = spot_track(beta, pa, theta_s, phi_s)
    ax.plot(
        track[:, 0],
        track[:, 1],
        ls="--",
        lw=0.9,
        color="C0",
        alpha=0.85,
        zorder=2,
        label="locus over one turn",
    )

    # Sky-projected spin axis: position angle PA, length sin(beta).  The legend
    # entry is made by an empty line, so the label cannot collide with the arrow.
    axis_len = np.sin(beta)
    axis_xy = axis_len * np.array([np.sin(pa), np.cos(pa)])
    if axis_len > 1e-9:  # degenerate when pole-on: the projection genuinely vanishes
        ax.annotate(
            "",
            xy=axis_xy,
            xytext=(0.0, 0.0),
            arrowprops=dict(arrowstyle="-|>", color="C3", lw=1.8),
            zorder=3,
        )
    ax.plot([], [], color="C3", lw=1.8, label="projected spin axis (length = sin beta)")

    # The spot itself.
    v = spot_observer_direction(alpha, beta, pa, theta_s, phi_s)
    theta_o, phi_o = observer_angles(v)
    visible = v[2] > 0.0
    if visible:
        ax.plot([v[0]], [v[1]], marker="o", ms=11, mfc="C1", mec="k", mew=1.0, zorder=5)
    else:
        ax.plot(
            [v[0]],
            [v[1]],
            marker="o",
            ms=11,
            mfc="none",
            mec="0.45",
            mew=1.8,
            zorder=5,
        )

    # Live invariant: the spot keeps a fixed angular distance from the spin axis.
    axis_v = spin_axis_observer_direction(beta, pa)
    angle_to_axis = np.degrees(np.arccos(np.clip(np.dot(v, axis_v), -1.0, 1.0)))
    pa_spot = np.degrees(position_angle(v[0], v[1]))

    lines = [
        f"spot in body frame   theta_S = {np.degrees(theta_s):7.2f}   phi_S = {np.degrees(phi_s):7.2f}",
        f"spot in obs. frame   theta_O = {np.degrees(theta_o):7.2f}   phi_O = {np.degrees(phi_o):7.2f}",
        f"sky offset           u = {v[0]:+7.4f} East   v = {v[1]:+7.4f} North",
        f"spot position angle  PA_spot = {pa_spot:7.2f}   rho = {np.hypot(v[0], v[1]):6.4f}",
        f"mu = cos(theta_O)    {v[2]:+7.4f}   ->  {'near side (visible)' if visible else 'far side (hidden)'}",
        f"angle to spin axis   {angle_to_axis:7.2f}   (constant = theta_S)",
    ]
    ax.text(
        0.015,
        0.985,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        ha="left",
        family="monospace",
        fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", alpha=0.92),
        zorder=6,
    )
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)


def main():
    theta_s = np.radians(SPOT_THETA_S_DEG)
    phi_s = np.radians(SPOT_PHI_S_DEG)

    fig = plt.figure(figsize=(7.4, 8.4))
    ax = fig.add_axes([0.09, 0.30, 0.84, 0.64])

    ax_alpha = fig.add_axes([0.26, 0.185, 0.60, 0.028])
    ax_beta = fig.add_axes([0.26, 0.125, 0.60, 0.028])
    ax_pa = fig.add_axes([0.26, 0.065, 0.60, 0.028])

    s_alpha = Slider(ax_alpha, r"$\alpha$  spin phase", 0.0, 360.0, valinit=0.0)
    s_beta = Slider(ax_beta, r"$\beta$  inclination", 0.0, 180.0, valinit=60.0)
    s_pa = Slider(ax_pa, r"PA   position angle (N$\rightarrow$E)", 0.0, 360.0, valinit=90.0)

    def update(_value):
        draw(
            ax,
            np.radians(s_alpha.val),
            np.radians(s_beta.val),
            np.radians(s_pa.val),
            theta_s,
            phi_s,
        )
        fig.canvas.draw_idle()

    s_alpha.on_changed(update)
    s_beta.on_changed(update)
    s_pa.on_changed(update)

    update(None)
    plt.show()


if __name__ == "__main__":
    main()
