"""Special-function wrappers and the classical anchors.

Everything here exists to keep convention decisions in one place.  The two
anchors -- the uniform disk, ``2 J1(x)/x``, and the linearly limb-darkened disk
-- are the end-to-end checks that the whole convention stack is right.

Naming collision to avoid: the limb-darkening literature calls its coefficients
``u``, ``u1``, ``u2``, whereas ``u`` here is a baseline coordinate.  In this
package limb-darkening coefficients are named ``ld_coeffs``.
"""
