"""Emit the LaTeX fragments that the manual inputs.

Run from the repository root:

    python -m mostlyharmless.symbolic.emit

This is the only module that knows about LaTeX.  Everything else in the package is
numerics, and the coefficients it prints are computed just in time rather than read
from a file.  The output is committed; a test regenerates it in memory and fails on
a mismatch.
"""

from __future__ import annotations

from pathlib import Path

import sympy as sp

from mostlyharmless.visylm.zernike import zernike_coefficients

__all__ = ["coefficients_table", "write"]

L_MAX = 3


def coefficients_table(l_max=L_MAX):
    """The even-family coefficients as a LaTeX table, to l_max."""
    lines = [
        r"\begin{center}",
        r"\begin{tabular}{llll}",
        r"\hline",
        r"$l$ & $m$ & $k$ & $c_k$ \\",
        r"\hline",
    ]
    for l in range(l_max + 1):
        for m in range(l + 1):
            if (l + m) % 2 == 1:
                continue
            for k, value in enumerate(zernike_coefficients(l, m)):
                lines.append(f"${l}$ & ${m}$ & ${k}$ & ${sp.latex(value)}$ \\\\")
    lines += [r"\hline", r"\end{tabular}", r"\end{center}"]
    return "\n".join(lines) + "\n"


def write(directory=None):
    """Write the fragments, and return the paths written."""
    if directory is None:
        root = Path(__file__).resolve().parents[3]
        directory = root / "docs" / "generated"
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "coefficients.tex"
    path.write_text(coefficients_table())
    return [path]


if __name__ == "__main__":
    for written in write():
        print(written)
