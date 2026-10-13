"""Spherical-harmonic modes and their interferometric visibilities.

Modules planned (ROADMAP.md section 2):

``spherical_harmonics``  definitions, evaluation, orthonormality      (F1-F3)
``coefficients``         the Coeffs container for a real map          (D13)
``wigner``               d and D matrices, coefficient rotation       (F4-F8)
``zernike``              radial polynomials and the exact even-family
                         coefficients, computed just in time           (F14-F15)
``visibility``           the master reduction, both parity families,
                         and the classical anchors                     (F9-F18)

There is no stored coefficient table and no data file (D19): the even-family
coefficients are exact rationals computed on demand and memoised.

Conventions are fixed in ``docs/sections/01-conventions.tex`` and are not
restated here.

Public entry points take ``(u, v)`` in wavelengths.  ``q`` and ``psi`` are
derived internally and do not appear in a caller's code.
"""
