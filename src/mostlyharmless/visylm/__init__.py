"""Spherical-harmonic modes and their interferometric visibilities.

Modules planned (ROADMAP.md section 2):

``spherical_harmonics``  definitions, evaluation, orthonormality      (F1-F3)
``wigner``               d and D matrices, coefficient rotation       (F4-F8)
``zernike``              radial polynomials built from Jacobi         (F14)
``visibility``          the master reduction and both parity families (F9-F18)
``tables``               the exact rational coefficient table

Conventions are fixed in ``docs/sections/01-conventions.tex`` and are not
restated here.

Public entry points take ``(u, v)`` in wavelengths.  ``q`` and ``psi`` are
derived internally and do not appear in a caller's code.
"""
