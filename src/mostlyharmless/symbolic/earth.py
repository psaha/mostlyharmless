"""The Earth, as three coefficient sets, extracted once from a colour image.

The map GUI can show a recognisable planet instead of a circular spot, and the
coefficients cannot be synthesised: the whole point is that a user recognises the
continents.  So they are extracted here, from an equirectangular image, and then committed
as numbers by ``make earth``.  Nothing at run time reads an image, and nothing at run time
depends on pyshtools.

An equirectangular image *is* a lat-lon grid, so the projection is the identity: pixel
$(i,j)$ is $(\theta_i, \varphi_j)$ with uniform spacing.  The analysis is a quadrature,

    a_lm = sum_ij sin(theta_i) dtheta dphi f(i,j) real_sph_harm(l, m, theta_i, phi_j),

with ``real_sph_harm`` being our own real basis, so each channel feeds ``Coeffs.from_real``
directly and every downstream consumer --- the renderer, the visibility, the rotation ---
takes it unchanged.

Three channels, because the image is colour: red, green and blue are carried as three
independent brightness maps and recombined only at the end.  Land, ice and ocean then differ
by colour rather than by shape alone, which is what makes a small image readable.

A midpoint grid is not exact for a band-limited image --- it is what a picture gives --- so
the low degrees are the trustworthy ones.  A degree of about 40 resolves roughly 4.5 degrees,
which is enough for the continents to be recognisable and is the value used here.
"""

from __future__ import annotations

import numpy as np

from mostlyharmless.visylm.coefficients import Coeffs
from mostlyharmless.visylm.spherical_harmonics import real_sph_harm

__all__ = ["read_grid", "channel", "channels", "taper", "truncate", "emit",
           "timing", "main"]

DEFAULT_L_MAX = 40


def read_grid(path, n_lat=180, n_lon=360):
    """The image as a ``(n_lat, n_lon, 3)`` array of colour in ``[0, 1]``.

    Downsampled by block-averaging, since at ``l ~ 40`` the finest feature is about four
    degrees and the full-resolution image is only being averaged anyway.  ``matplotlib``
    reads row 0 from the top of the image, which for a map is the north pole, so there is
    no flip.
    """
    from matplotlib import pyplot as plt

    image = np.asarray(plt.imread(path), dtype=float)
    if image.ndim == 2:                                 # a greyscale image
        image = np.stack([image] * 3, axis=-1)
    image = image[..., :3]
    if image.max() > 1.5:                               # 8-bit
        image = image / 255.0
    rows, columns = image.shape[:2]
    lat_edges = np.linspace(0, rows, n_lat + 1).astype(int)
    lon_edges = np.linspace(0, columns, n_lon + 1).astype(int)
    grid = np.empty((n_lat, n_lon, 3))
    for i in range(n_lat):
        block = image[lat_edges[i]:lat_edges[i + 1]]
        for j in range(n_lon):
            cell = block[:, lon_edges[j]:lon_edges[j + 1], :]
            grid[i, j] = cell.reshape(-1, 3).mean(axis=0)
    return grid


def channel(values, l_max=DEFAULT_L_MAX):
    """Project one colour plane onto the real basis.

    ``theta`` is the colatitude, so the area element is ``sin(theta) dtheta dphi``; a
    ``cos`` weight would be the latitude weight and vanishes on a symmetric grid.
    """
    values = np.asarray(values, dtype=float)
    n_lat, n_lon = values.shape
    theta = (np.arange(n_lat) + 0.5) * np.pi / n_lat
    phi = (np.arange(n_lon) + 0.5) * 2.0 * np.pi / n_lon - np.pi
    weight = np.sin(theta)[:, None] * (np.pi / n_lat) * (2.0 * np.pi / n_lon)

    weighted = values * weight
    theta_grid = np.broadcast_to(theta[:, None], (n_lat, n_lon))
    phi_grid = np.broadcast_to(phi[None, :], (n_lat, n_lon))

    plan = [np.zeros(2 * l + 1) for l in range(l_max + 1)]
    for l in range(l_max + 1):
        for m in range(-l, l + 1):
            plan[l][m + l] = np.sum(weighted * real_sph_harm(l, m, theta_grid, phi_grid))
    return Coeffs.from_real(plan)


def channels(path, l_max=DEFAULT_L_MAX, n_lat=180, n_lon=360):
    """``(red, green, blue)`` as three coefficient sets, in the body frame."""
    grid = read_grid(path, n_lat=n_lat, n_lon=n_lon)
    return tuple(channel(grid[..., k], l_max=l_max) for k in range(3))


def taper(coeffs, l_max=DEFAULT_L_MAX):
    """Lanczos σ-factors, applied per degree.

    A colour image of the Earth is close to a jump at every coastline, so its coefficients
    fall slowly and the ringing at the truncation is severe --- without this the continents
    are not visible at all, only the ripples.  The taper trades the ripples for a soft coast.
    """
    blocks = []
    for l, block in enumerate(coeffs.blocks):
        sigma = 1.0 if l == 0 else float(np.sinc(l / (l_max + 1.0)))
        blocks.append(block * sigma)
    return Coeffs(tuple(blocks))


def truncate(coeffs, l_max):
    """The same map at a lower degree, for a cheaper render.

    The emitted set is as high as the extraction went; the GUI may want less, and at l = 40
    six renders per redraw (three channels, two faces) is slow enough to feel stuck.
    """
    blocks = tuple(block for block in coeffs.blocks[:l_max + 1])
    return Coeffs(blocks)


def write(red, green, blue, l_max=DEFAULT_L_MAX, target=None):
    """Write three coefficient sets as a Python module of literals.

    Committed, so that the package keeps its no-data-file property and never depends on
    pyshtools: the extraction is a one-off, the numbers are the product.

    Split from :func:`emit` so that the writer can be round-tripped on a set whose
    coefficients are known --- which is how the missing imaginary parts were going to be
    caught the second time.
    """
    from pathlib import Path

    target = Path(target or "src/mostlyharmless/symbolic/earth_coefficients.py")
    lines = [
        '"""The Earth as three coefficient sets, emitted by ``make earth``.',
        "",
        "Do not edit: regenerate instead.  The extraction is in ``earth.py``; this file",
        "carries the numbers so that the package needs neither the image nor pyshtools.",
        f"Extracted at l_max = {l_max} from a Blue Marble image; NASA, public domain.",
        '"""',
        "",
        "L_MAX = " + repr(int(l_max)),
        "",
    ]
    for name, coeffs in (("RED", red), ("GREEN", green), ("BLUE", blue)):
        lines.append(f"{name} = (")
        for block in coeffs.blocks:
            # complex, not just the real part: in the half-complex form the imaginary part
            # carries the sine content, and dropping it makes the map even under
            # phi -> -phi --- a planet symmetric about its own meridian, which is wrong
            values = ", ".join(f"complex({float(value.real):.6e}, {float(value.imag):.6e})"
                               for value in block)
            lines.append(f"    ({values}),")
        lines.append(")")
        lines.append("")
    target.write_text("\n".join(lines))
    return target


def emit(path, l_max=DEFAULT_L_MAX, target=None):
    """Extract three colour channels from an image and write them as a module."""
    return write(*channels(path, l_max=l_max), l_max=l_max, target=target)


def timing(l_values=(4, 8, 16, 32), pixels_values=(121, 241, 361)):  # pragma: no cover
    """Cost of the render path, per degree and per resolution.

    The Earth needs six panels per redraw --- three colour channels by two faces --- so the
    renderer's cost is six times one panel's.  This breaks the redraw into its parts, so that
    a slow setting can be attributed rather than guessed at: the rotation is Wigner-D, the
    panel is one ``sky_map`` over every mode, and the pair is what a channel costs.
    """
    import time

    from mostlyharmless.symbolic import map_gui, maps
    from mostlyharmless.visylm.wigner import rotate_coefficients

    red = map_gui.earth_coefficients()[0]        # one channel is representative
    print(f"{'l':>3} {'px':>4} {'rotate':>8} {'warm':>6} {'1 panel':>8} {'2 faces':>8} {'all 6':>8}")
    for l_max in l_values:
        for pixels in pixels_values:
            channel = truncate(red, l_max)
            start = time.perf_counter()
            faced = rotate_coefficients(channel, 0.3, 0.9, 0.2)
            rotate_time = time.perf_counter() - start

            # a second rotation at the same inclination: the difference between the two is
            # the d-matrix build, and what remains is the per-mode work
            start = time.perf_counter()
            rotate_coefficients(channel, 0.5, 0.9, 0.2)
            warm_time = time.perf_counter() - start

            start = time.perf_counter()
            maps.sky_map(faced, n=pixels)
            one_panel = time.perf_counter() - start

            start = time.perf_counter()
            maps.sky_map(faced, n=pixels)
            maps.sky_map(faced, back=True, n=pixels)
            both_faces = time.perf_counter() - start

            start = time.perf_counter()
            map_gui.compute_map(map_gui.MapParameters(l_max=l_max), pixels=pixels,
                                model="Earth")
            all_six = time.perf_counter() - start

            print(f"{l_max:>3} {pixels:>4} {rotate_time:>8.3f} {warm_time:>6.3f} "
                  f"{one_panel:>8.3f} {both_faces:>8.3f} {all_six:>8.3f}")


def main(argv=None):  # pragma: no cover - command line
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image", nargs="?",
                        help="an equirectangular global colour image, north at top")
    parser.add_argument("--time", action="store_true",
                        help="time the render path instead of extracting")
    parser.add_argument("--lmax", type=int, default=DEFAULT_L_MAX)
    parser.add_argument("--target", default=None)
    arguments = parser.parse_args(argv)

    if arguments.time:
        timing()
        return
    if not arguments.image:
        parser.error("an image is needed unless --time is given")
    written = emit(arguments.image, arguments.lmax, arguments.target)
    print(f"wrote {written} at l_max = {arguments.lmax}")


if __name__ == "__main__":  # pragma: no cover
    main()
