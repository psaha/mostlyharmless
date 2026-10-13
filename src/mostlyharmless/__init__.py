"""mostlyharmless -- analytic interferometric visibilities of surface maps.

The public surface lives in two subpackages:

``mostlyharmless.visylm``
    Spherical harmonics, Wigner rotation matrices, Zernike polynomials, the
    visibility chain, and the classical anchors used as end-to-end checks.
``mostlyharmless.symbolic``
    SymPy derivations that back the machine-checked lemmas, and the emitters for
    the coefficient table and the LaTeX fragments.

Bessel and Legendre functions come straight from scipy and sympy.  There is no
wrapper module for them, because there is no convention of ours for it to
protect.

Submodules are not imported eagerly; the caller imports what it needs.
"""

__version__ = "0.0.0"
