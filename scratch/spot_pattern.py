"""
Spherical-harmonic spot: build a band-limited spot at the body North pole,
rotate it off the pole with a fixed initial rotation, then rotate the
coefficients into the observer's frame with Wigner D matrices and display the
pattern on the visible disk.

Throwaway exploration code, but the Wigner machinery and the coefficient
rotation are written as separable functions so they can graduate into the
package later.

Conventions (same as scratch/spot_test.py)
------------------------------------------
Observer frame (u, v, w) = (East, North, line of sight), w pointing from the
observer towards the star, so the visible disk is the unit circle and the near
hemisphere is w > 0.  East is to the LEFT: the u axis is inverted on the plot.

Body frame (x, y, z): z is the spin axis.

The body-to-observer rotation is

    R(alpha, beta, PA) = R_z(pi/2 - PA) @ R_y(beta) @ R_z(alpha)

with alpha the spin phase, beta the inclination from the line of sight, and PA
the position angle of the projected spin axis measured from North through East.

Coefficient transformation (verified against a direct geometric rotation of the
function on the sphere, to 1e-15):

    b_{l m'} = sum_m D^l_{m' m}(R) a_{l m}
    D^l_{m' m}(R) = exp(-i m' (pi/2 - PA)) d^l_{m' m}(beta) exp(-i m alpha)

coefficients indexed by m + l, m ascending, and NO complex conjugation -- R here
is body -> observer.

Pipeline
--------
1. spot at the body pole:  a0_{l0} = sqrt((2l+1)/4pi), a0_{lm} = 0 otherwise.
   That is a band-limited delta, i.e. the spherical Dirichlet kernel: a broad
   main lobe with negative rings.
2. fixed initial rotation R_y(theta0) moves the spot to body colatitude theta0.
   This matters: with the spot ON the spin axis, alpha would do nothing at all,
   because R_z(alpha) leaves the axis invariant.
3. slider rotation R(alpha, beta, PA) takes it to the observer frame.
4. synthesise on the observer's disk and show both hemispheres.

Because the spot sits at angular distance theta0 from the spin axis, alpha
carries it around a cone of half-angle theta0: the blob orbits the projected
spin axis, running out through the limb and back in on the far side.
"""

from __future__ import annotations

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider
from scipy.linalg import expm
from scipy.special import lpmv

# --------------------------------------------------------------------- inputs
LMAX = 10  # band limit
SPOT_BODY_COLATITUDE_DEG = 45.0  # fixed initial rotation: spot off the pole
GRID_N = 129  # samples across the disk, per panel
COLOUR_LIMIT = 1.0  # symmetric colour scale; lower it (e.g. 0.3) to show the rings


# =====================================================================
# Wigner d matrices
# =====================================================================
def jy_matrix(l):
    """Angular-momentum operator J_y in the |l m> basis, m ascending."""
    n = 2 * l + 1
    J = np.zeros((n, n), dtype=complex)
    for m in range(-l, l + 1):
        if m + 1 <= l:
            J[m + 1 + l, m + l] = np.sqrt((l - m) * (l + m + 1)) / 2j
        if m - 1 >= -l:
            J[m - 1 + l, m + l] = -np.sqrt((l + m) * (l - m + 1)) / 2j
    return J


def wigner_d_expm(l, beta):
    r"""d^l_{m'm}(beta) as exp(-i beta J_y); the definition, no formula involved.

    Returned as an array indexed [m' + l, m + l].
    """
    return expm(-1j * beta * jy_matrix(l)).real


def wigner_d_explicit(l, beta):
    r"""d^l_{m'm}(beta) from the explicit finite sum, indexed [m' + l, m + l].

    .. math::

       d^l_{m'm}(\beta) = \sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}\;
         \sum_s \frac{(-1)^{m'-m+s}}{s!\,(l+m-s)!\,(l-m'-s)!\,(m'-m+s)!}\;
         \left(\cos\tfrac\beta2\right)^{2l+m-m'-2s}
         \left(\sin\tfrac\beta2\right)^{m'-m+2s}

    with s running over all values for which every factorial argument is
    non-negative.  Stable to machine precision for l <~ 20; beyond that the
    alternating sum loses digits and a recurrence is needed.
    """
    n = 2 * l + 1
    out = np.zeros((n, n))
    cb, sb = math.cos(beta / 2.0), math.sin(beta / 2.0)
    for mp in range(-l, l + 1):
        for m in range(-l, l + 1):
            pref = math.sqrt(
                math.factorial(l + m)
                * math.factorial(l - m)
                * math.factorial(l + mp)
                * math.factorial(l - mp)
            )
            s_min = max(0, m - mp)
            s_max = min(l + m, l - mp)
            total = 0.0
            for s in range(s_min, s_max + 1):
                den = (
                    math.factorial(s)
                    * math.factorial(l + m - s)
                    * math.factorial(l - mp - s)
                    * math.factorial(mp - m + s)
                )
                total += (
                    (-1.0) ** (mp - m + s)
                    / den
                    * cb ** (2 * l + m - mp - 2 * s)
                    * sb ** (mp - m + 2 * s)
                )
            out[mp + l, m + l] = pref * total
    return out


def check_wigner(lmax=LMAX, n_beta=7):
    """Largest disagreement between the explicit sum and exp(-i beta J_y)."""
    worst = 0.0
    for l in range(lmax + 1):
        for beta in np.linspace(0.0, np.pi, n_beta):
            worst = max(
                worst,
                float(np.abs(wigner_d_explicit(l, beta) - wigner_d_expm(l, beta)).max()),
            )
    return worst


# =====================================================================
# Rotations of the coefficient vector
# =====================================================================
def rot_z(angle):
    """Active right-handed rotation by `angle` about +z."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rot_y(angle):
    """Active right-handed rotation by `angle` about +y."""
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def rotation_matrix(alpha, beta, pa):
    """R(alpha, beta, PA) taking body-frame components to observer-frame ones."""
    return rot_z(0.5 * np.pi - pa) @ rot_y(beta) @ rot_z(alpha)


def pole_spot_coefficients(lmax=LMAX):
    """Band-limited delta at the body North pole, normalised to unit peak.

    a_{l0} = sqrt((2l+1)/4pi), all other m zero; the synthesis peaks at
    sum_l (2l+1)/4pi = (lmax+1)^2/4pi, which is divided out.
    """
    peak = (lmax + 1) ** 2 / (4.0 * np.pi)
    out = []
    for l in range(lmax + 1):
        c = np.zeros(2 * l + 1, dtype=complex)
        c[l] = math.sqrt((2 * l + 1) / (4.0 * np.pi)) / peak
        out.append(c)
    return out


def rotate_coefficients(a, alpha, beta, pa, d_matrix=wigner_d_explicit):
    """Apply b_l = D^l(R)^T a_l with R = R(alpha, beta, pa).

    `a` is a list of length lmax+1 of complex arrays indexed by m + l.
    """
    out = []
    for l, al in enumerate(a):
        ms = np.arange(-l, l + 1)
        d = d_matrix(l, beta)
        phase_observer = np.exp(-1j * ms * (0.5 * np.pi - pa))  # on m'
        phase_body = np.exp(-1j * ms * alpha)  # on m
        operator = phase_observer[:, None] * d * phase_body[None, :]
        out.append(operator @ al)
    return out


def flatten_coefficients(a):
    """Pack the ragged coefficient list into one vector, offset(l) = l^2."""
    out = np.empty((len(a)) ** 2, dtype=complex)
    for l, al in enumerate(a):
        out[l * l : l * l + 2 * l + 1] = al
    return out


# =====================================================================
# Spherical harmonics and synthesis on the observer's disk
# =====================================================================
def spherical_harmonic(l, m, theta, phi):
    r"""Y_l^m(theta, phi): orthonormal, Condon-Shortley phase, theta polar.

    Written from the definition rather than calling a version-specific SciPy
    routine, so it works across SciPy releases.
    """
    mp = abs(int(m))
    norm = math.sqrt(
        (2 * l + 1)
        / (4.0 * np.pi)
        * math.factorial(l - m)
        / math.factorial(l + m)
    )
    p = lpmv(mp, l, np.clip(np.cos(theta), -1.0, 1.0))
    if m < 0:
        p = (-1.0) ** mp * (math.factorial(l - mp) / math.factorial(l + mp)) * p
    return norm * p * np.exp(1j * m * phi)


def build_bases(lmax=LMAX, n=GRID_N):
    """Precompute everything that does not change as the sliders move.

    Returns a dict with the near and far hemisphere basis matrices (one row per
    grid point, one column per coefficient), the Cartesian grid, and the disk
    mask.  Building these once is what keeps slider dragging responsive.
    """
    x = np.linspace(-1.0, 1.0, n)
    xx, yy = np.meshgrid(x, x)  # xx[i,j] = x[j] (East), yy[i,j] = y[i] (North)
    rho = np.hypot(xx, yy)
    inside = rho <= 1.0
    rho_c = np.clip(rho, 0.0, 1.0)
    theta_near = np.arcsin(rho_c)  # near hemisphere: w > 0
    theta_far = np.pi - theta_near  # far hemisphere:  w < 0
    phi = np.arctan2(yy, xx)

    bases = {}
    for name, theta in (("near", theta_near), ("far", theta_far)):
        basis = np.empty((theta.size, (lmax + 1) ** 2), dtype=complex)
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                basis[:, l * l + (m + l)] = spherical_harmonic(
                    l, m, theta.ravel(), phi.ravel()
                )
        bases[name] = basis
    return {
        "bases": bases,
        "shape": (n, n),
        "xx": xx,
        "yy": yy,
        "inside": inside,
    }


def synthesise(basis, coefficients):
    """Evaluate sum b_lm Y_lm on the grid; the map is real by construction."""
    return (basis @ coefficients).real


# =====================================================================
# Drawing
# =====================================================================
def compass(ax):
    """E/W/N/S letters just inside the limb.  E sits on the left because the
    u axis is inverted, matching the sky."""
    for label, (cu, cv) in (
        ("E", (0.92, 0.0)),
        ("W", (-0.92, 0.0)),
        ("N", (0.0, 0.92)),
        ("S", (0.0, -0.92)),
    ):
        ax.text(
            cu,
            cv,
            label,
            ha="center",
            va="center",
            fontsize=9,
            color="0.15",
            bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none", alpha=0.55),
            zorder=6,
        )


def draw(
    ax_near,
    ax_far,
    readout,
    a0,
    grid,
    d_matrix,
    alpha,
    beta,
    pa,
    theta0,
):
    """Rotate, synthesise, and redraw both panels for one set of angles."""
    # 1. coefficients: pole spot -> fixed initial rotation -> slider rotation.
    #    R_y(theta0) is R(alpha=0, beta=theta0, PA=pi/2), so the same code path.
    a = rotate_coefficients(a0, 0.0, theta0, 0.5 * np.pi, d_matrix=d_matrix)
    b = rotate_coefficients(a, alpha, beta, pa, d_matrix=d_matrix)
    bvec = flatten_coefficients(b)

    # 2. geometry: where the spin axis and the spot point.
    R = rotation_matrix(alpha, beta, pa)
    s_pole = R @ np.array([0.0, 0.0, 1.0])
    s_spot = R @ (rot_y(theta0) @ np.array([0.0, 0.0, 1.0]))

    # 3. synthesis on both hemispheres.
    images = {}
    for name, ax in (("near", ax_near), ("far", ax_far)):
        values = synthesise(grid["bases"][name], bvec).reshape(grid["shape"])
        images[name] = values
        masked = np.ma.masked_where(~grid["inside"], values)
        ax.clear()
        ax.imshow(
            masked,
            extent=(-1.0, 1.0, -1.0, 1.0),
            origin="lower",
            cmap="RdBu_r",
            vmin=-COLOUR_LIMIT,
            vmax=COLOUR_LIMIT,
            interpolation="bilinear",
            zorder=0,
        )
        ax.add_patch(
            Circle((0.0, 0.0), 1.0, facecolor="none", edgecolor="0.25", lw=1.2, zorder=3)
        )
        ax.axhline(0.0, color="0.6", lw=0.5, zorder=2)
        ax.axvline(0.0, color="0.6", lw=0.5, zorder=2)
        compass(ax)
        ax.set_aspect("equal")
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.15, 1.15)
        ax.invert_xaxis()  # u is East, and East is to the left
        ax.set_xlabel("u   (East)")
        ax.set_ylabel("v   (North)", rotation=0, labelpad=10)

        # Projected spin axis: a sky-plane direction, so it appears in both
        # panels.  Solid where the pole actually is, dotted otherwise.
        axis_xy = math.sin(beta) * np.array([math.sin(pa), math.cos(pa)])
        pole_here = (s_pole[2] > 0.0) == (name == "near")
        if math.sin(beta) > 1e-9:
            ax.annotate(
                "",
                xy=axis_xy,
                xytext=(0.0, 0.0),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color="crimson",
                    lw=1.8,
                    linestyle="-" if pole_here else ":",
                ),
                zorder=4,
            )

        # The spot itself: a point on the sphere, so it lives in one panel only.
        if (s_spot[2] > 0.0) == (name == "near"):
            ax.plot(
                [s_spot[0]],
                [s_spot[1]],
                marker="o",
                ms=7,
                mfc="none",
                mec="black",
                mew=1.4,
                zorder=5,
            )

    ax_near.set_title("near hemisphere (what we see)", fontsize=10)
    ax_far.set_title("far hemisphere (same sky axes)", fontsize=10)

    # 4. quantitative check: does the image peak sit where geometry says?
    peak = {}
    for name in ("near", "far"):
        values = np.where(grid["inside"], images[name], -np.inf)
        idx = int(np.argmax(values))
        ii, jj = np.unravel_index(idx, grid["shape"])
        peak[name] = (float(grid["xx"][ii, jj]), float(grid["yy"][ii, jj]), float(values[ii, jj]))
    best = "near" if peak["near"][2] >= peak["far"][2] else "far"
    px, py, pval = peak[best]
    pz = math.sqrt(max(0.0, 1.0 - px * px - py * py))
    if best == "far":
        pz = -pz
    peak_dir = np.array([px, py, pz])
    offset_deg = math.degrees(
        math.acos(float(np.clip(np.dot(peak_dir, s_spot), -1.0, 1.0)))
    )

    def where(v):
        return "near" if v[2] > 0.0 else "far"

    def angles(v):
        return math.degrees(math.acos(float(np.clip(v[2], -1.0, 1.0)))), math.degrees(
            math.atan2(v[1], v[0])
        )

    th_p, ph_p = angles(s_pole)
    th_s, ph_s = angles(s_spot)
    lines = [
        f"lmax = {LMAX}    spot body colatitude theta0 = {math.degrees(theta0):6.2f} deg",
        f"alpha = {math.degrees(alpha) % 360:7.2f}   beta = {math.degrees(beta):7.2f}   PA = {math.degrees(pa) % 360:7.2f}   (degrees)",
        "",
        f"spot      theta_O = {th_s:7.2f}  phi_O = {ph_s:8.2f}   u = {s_spot[0]:+7.4f}  v = {s_spot[1]:+7.4f}",
        f"          rho = {math.hypot(s_spot[0], s_spot[1]):6.4f}  mu = {s_spot[2]:+7.4f}  ->  {where(s_spot)} hemisphere",
        f"spin axis theta_O = {th_p:7.2f}  phi_O = {ph_p:8.2f}   ->  {where(s_pole)} hemisphere",
        "",
        f"image peak in the {best} panel at u = {px:+7.4f}  v = {py:+7.4f}"
        f"   offset from the geometric spot: {offset_deg:5.2f} deg",
        f"angular distance from spot to spin axis: "
        f"{math.degrees(math.acos(float(np.clip(np.dot(s_spot, s_pole), -1.0, 1.0)))):6.2f} deg"
        f"   (constant = theta0)",
    ]
    readout.set_text("\n".join(lines))


# =====================================================================
def main():
    # Self-check before anything is drawn: the closed form against the definition.
    err = check_wigner(LMAX)
    d_matrix = wigner_d_explicit
    print(
        f"Wigner d self-check (explicit sum vs exp(-i beta J_y)), l <= {LMAX}: "
        f"max |difference| = {err:.2e}"
    )
    if err < 1e-10:
        print("  -> using wigner_d_explicit")
    else:
        d_matrix = wigner_d_expm
        print("  -> DISAGREEMENT: falling back to wigner_d_expm (please report this)")

    theta0 = math.radians(SPOT_BODY_COLATITUDE_DEG)
    a0 = pole_spot_coefficients(LMAX)
    grid = build_bases(LMAX, GRID_N)
    print(f"built {GRID_N}x{GRID_N} bases for l <= {LMAX}; starting the window")

    fig = plt.figure(figsize=(12.5, 9.6))
    ax_near = fig.add_axes([0.05, 0.42, 0.41, 0.52])
    ax_far = fig.add_axes([0.54, 0.42, 0.41, 0.52])
    readout = fig.text(0.05, 0.385, "", va="top", ha="left", family="monospace", fontsize=8.5)

    ax_alpha = fig.add_axes([0.30, 0.215, 0.55, 0.022])
    ax_beta = fig.add_axes([0.30, 0.145, 0.55, 0.022])
    ax_pa = fig.add_axes([0.30, 0.075, 0.55, 0.022])

    # alpha = 120 deg rather than 0: with theta0 = 70 and beta = 60 the spot
    # starts behind the limb at alpha = 0, which is a poor first view.
    s_alpha = Slider(ax_alpha, r"$\alpha$  spin phase", 0.0, 360.0, valinit=120.0)
    s_beta = Slider(ax_beta, r"$\beta$  inclination", 0.0, 180.0, valinit=60.0)
    s_pa = Slider(ax_pa, r"PA   position angle (N$\rightarrow$E)", 0.0, 360.0, valinit=90.0)

    def update(_value):
        draw(
            ax_near,
            ax_far,
            readout,
            a0,
            grid,
            d_matrix,
            math.radians(s_alpha.val),
            math.radians(s_beta.val),
            math.radians(s_pa.val),
            theta0,
        )
        fig.canvas.draw_idle()

    s_alpha.on_changed(update)
    s_beta.on_changed(update)
    s_pa.on_changed(update)

    update(None)
    plt.show()


if __name__ == "__main__":
    main()
