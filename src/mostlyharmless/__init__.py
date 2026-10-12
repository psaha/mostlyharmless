"""mostlyharmless -- analytic interferometric visibilities of surface maps.

The public surface lives in three subpackages:

``mostlyharmless.visylm``
    Spherical harmonics, Wigner rotation matrices, Zernike polynomials and the
    visibility chain -- the physics.
``mostlyharmless.special``
    Thin wrappers over the special functions, with conventions documented, plus
    the classical anchors used as end-to-end checks.
``mostlyharmless.symbolic``
    SymPy derivations that back the machine-checked lemmas, and the emitters for
    the coefficient table and the LaTeX fragments.

Submodules are not imported eagerly; the caller imports what it needs.
"""

__version__ = "0.0.0"
