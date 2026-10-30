"""The ``Coeffs`` container: spherical-harmonic coefficients of a real map.

A real map is expanded as

    I(theta, phi) = sum_{l=0}^{L} sum_{m=-l}^{l} a_lm Y_{lm}(theta, phi) ,

but only the non-negative orders are stored.  The negative orders are derived from
the reality condition

    a_{l,-m} = conj(a_lm) ,

so the map is real *by construction*.  There is no stored negative order that can
disagree with its positive partner, no roundoff defect to police, and no
unphysical direction for an optimiser to wander into.

Storage
-------
One complex array per degree, of length l + 1:

    blocks[l][m] = a_lm ,   m = 0, ..., l

with a_l0 real (the reality condition forces it) and the rest complex.  The number
of independent real numbers is

    sum_{l=0}^{L} (2l + 1) = (L + 1)^2 ,

which is the count any real parameterisation of the same map must have, and is
what ``to_vector`` returns.

No magnitude |m| is ever used: the order is non-negative wherever it is stored,
and where a negative order is needed it is written as -m.

Real basis
----------
``from_real`` and ``to_real`` convert to and from the real cos/sin basis of
section 02, in which the coefficients are indexed by signed m running from -l to
l: cosine for m > 0, sine for m < 0, and a single real number at m = 0.  For a
real map,

    a_l0 = r_l0 ,
    a_lm = (r_lm - i r_{l,-m}) / sqrt(2) ,   m > 0 ,

with inverse

    r_lm = sqrt(2) Re a_lm ,   r_{l,-m} = -sqrt(2) Im a_lm .

That sign is the part nobody guesses.  It is the conversion documented by
pyshtools for the coefficients of a real function, which is a useful independent
confirmation.

A limitation worth knowing
--------------------------
``Coeffs`` describes a *real* map.  A single complex mode, such as one Y_{lm} on
its own, is not a real map and will be rejected; use the mode-level visibility
functions for that, which is what the derivations are written in terms of.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ["Coeffs", "negative_order"]

#: Relative tolerance for accepting an a_l0 that ought to be real, and for the
#: reality check that ``from_unfolded`` performs.  A correct rotation leaves a
#: defect at the roundoff level, some 1e-15, so this has room without being slack.
_TOL = 1e-12


def negative_order(block, m):
    """``a_{l,-m}`` for ``m > 0``, from the stored ``block[m]``.

    This is the whole of the negative-order convention, in one place: the map is real,
    so the order $-m$ is the conjugate of $+m$ and carries no phase.  Everything that
    needs a negative order --- the accessor, the unfold and the fold --- calls this, so
    a change of convention is one edit rather than five.
    """
    return np.conjugate(block[m])


@dataclass(frozen=True, eq=False)
class Coeffs:
    """Spherical-harmonic coefficients of a real map, stored for m >= 0 only.

    Instances are immutable, and so are the arrays they hold.
    """

    blocks: tuple

    # ---------------------------------------------------------------- setup

    def __post_init__(self):
        blocks = []
        for l, block in enumerate(self.blocks):
            arr = np.asarray(block, dtype=complex).reshape(-1)
            if arr.size != l + 1:
                raise ValueError(
                    f"block {l} has {arr.size} entries, expected {l + 1}")
            if abs(arr[0].imag) > _TOL * max(1.0, abs(arr[0])):
                raise ValueError(
                    f"block {l}: a_l0 must be real for a real map, "
                    f"got {arr[0]!r}")
            arr = arr.copy()
            arr[0] = arr[0].real
            arr.setflags(write=False)
            blocks.append(arr)
        object.__setattr__(self, "blocks", tuple(blocks))

    @classmethod
    def zeros(cls, lmax):
        """A zero map of the given maximum degree."""
        if lmax < 0:
            raise ValueError(f"lmax must be non-negative, got {lmax}")
        return cls(tuple(np.zeros(l + 1, dtype=complex) for l in range(lmax + 1)))

    @property
    def lmax(self):
        return len(self.blocks) - 1

    # ------------------------------------------------------------- access

    def __getitem__(self, key):
        """a_lm for any order -l <= m <= l, deriving the negative ones."""
        l, m = key
        if m > l or m < -l:
            raise ValueError(f"need -l <= m <= l, got l={l}, m={m}")
        if m >= 0:
            return complex(self.blocks[l][m])
        return complex(negative_order(self.blocks[l], -m))

    def __eq__(self, other):
        if not isinstance(other, Coeffs):
            return NotImplemented
        if self.lmax != other.lmax:
            return False
        return all(np.array_equal(a, b)
                   for a, b in zip(self.blocks, other.blocks))

    # ------------------------------------------- real basis and its vector

    def to_real(self):
        """Real cos/sin coefficients, one array per degree, index m + l."""
        out = []
        for l, block in enumerate(self.blocks):
            row = np.zeros(2 * l + 1)
            row[l] = block[0].real
            for m in range(1, l + 1):
                row[l + m] = math.sqrt(2.0) * block[m].real
                row[l - m] = -math.sqrt(2.0) * block[m].imag
            out.append(row)
        return tuple(out)

    @classmethod
    def from_real(cls, blocks):
        """Build from real cos/sin coefficients, index m + l.

        Reality holds exactly, because the coefficients are constructed from real
        inputs rather than projected onto them.
        """
        out = []
        for l, block in enumerate(blocks):
            row = np.asarray(block, dtype=float).reshape(-1)
            if row.size != 2 * l + 1:
                raise ValueError(
                    f"block {l} has {row.size} entries, expected {2 * l + 1}")
            a = np.zeros(l + 1, dtype=complex)
            a[0] = row[l]
            for m in range(1, l + 1):
                a[m] = (row[l + m] - 1j * row[l - m]) / math.sqrt(2.0)
            out.append(a)
        return cls(tuple(out))

    def to_vector(self):
        """The flat real parameter vector, of length (lmax + 1)^2.

        In starry's order: increasing degree, and within a degree m from -l to l.
        This is the form to hand to an optimiser.
        """
        return np.concatenate(self.to_real())

    @classmethod
    def from_vector(cls, vector):
        """Inverse of :meth:`to_vector`."""
        vector = np.asarray(vector, dtype=float).reshape(-1)
        lmax = int(round(math.sqrt(vector.size))) - 1
        if (lmax + 1) ** 2 != vector.size:
            raise ValueError(
                "a flat coefficient vector must have (lmax+1)^2 entries, "
                f"that is 1, 4, 9, 16, ...; got {vector.size}")
        blocks, i = [], 0
        for l in range(lmax + 1):
            blocks.append(vector[i:i + 2 * l + 1])
            i += 2 * l + 1
        return cls.from_real(tuple(blocks))

    # ------------------------------------------------ full-range complex form

    def unfolded(self):
        """Full-range complex blocks, index m + l, m from -l to l.

        What the rotation and the mode sum consume.  Derived on demand from the
        stored half, so the two can never disagree.
        """
        out = []
        for l, block in enumerate(self.blocks):
            row = np.zeros(2 * l + 1, dtype=complex)
            row[l:] = block
            for m in range(1, l + 1):
                row[l - m] = negative_order(block, m)
            out.append(row)
        return tuple(out)

    @classmethod
    def from_unfolded(cls, blocks, tol=_TOL):
        """Fold full-range complex blocks back onto the stored half.

        The non-negative orders are kept; the negative ones are checked against
        the reality condition and discarded.  Since the real subspace is invariant
        under rotation, a correct rotation passes and a broken one is caught here
        rather than three stages downstream.
        """
        out = []
        for l, block in enumerate(blocks):
            row = np.asarray(block, dtype=complex).reshape(-1)
            if row.size != 2 * l + 1:
                raise ValueError(
                    f"block {l} has {row.size} entries, expected {2 * l + 1}")
            scale = max(1.0, float(np.max(np.abs(row))))
            for m in range(1, l + 1):
                want = negative_order(row, l + m)
                defect = abs(row[l - m] - want)
                if defect > tol * scale:
                    raise ValueError(
                        f"block {l}, order -{m}: violates the reality condition "
                        f"by {defect:.3e} (scale {scale:.3e})")
            out.append(row[l:].copy())
        return cls(tuple(out))
