"""The deformation field: one potential, its level surface, and the gravity there.

Stage 6, section 8.  The shape and the surface gravity of a deformed star are not
independent.  On a level surface the normal is along the gradient, so

    g_eff = |grad Phi| ,   delta r / a = (E - <E>) a^3 ,   delta g / g_0 = -2 delta r / a - 2 E a^3 ,

the last two to first order in the perturbation.  The potential is a monopole plus
one quadrupole per deformation, in units where GM = 1 and the unperturbed radius is 1:

    Phi(x) = -1/r - eps_tide  r^2 P_2(n . e_tide)
                  + eps_spin  r^2 P_2(n . e_spin)
                  - eps_spin  r^2 ,

both amplitudes positive for a tidal bulge and a centrifugal flattening.  The last
term is the l = 0 part of the centrifugal potential: that term is not a pure
quadrupole, and its l = 0 piece renormalises the effective monopole.

Amplitudes are quoted as A r_0^3 / GM, so eps = 0.2 is a twenty per cent deformation
of the potential at the reference radius, whatever axis ratio that produces.  The
surface is a level set of Phi, and the level value is a *choice*: it is the star's size.
``level_value`` fixes it by requiring the angular mean radius to be exactly 1, which makes
the first-order displacement delta r = E - <E> a pure quadrupole --- the l = 0 part of the
centrifugal term then shows up in the brightness alone, not in the shape.  The
volume-equivalent sphere is a separate quantity, since the deformation changes the volume.

**Range of validity.** The quadrupole potential grows as r^2 while the monopole falls
as 1/r, so past a modest amplitude the level set becomes unreachable along
some directions and the surface opens, as a Roche lobe does.  The binding direction is
the one where the perturbation is deepest, which for positive amplitudes is the tide
axis, where phi_1 = -E r^2 with

    E = eps_tide + 1.5 eps_spin .

The level condition along it, 1/r + E r^2 = 1, has its minimum at r^3 = 1/(2E), where
it equals (2^(1/3) + 2^(-2/3)) E^(1/3) = 1.89 E^(1/3).  The surface closes only while
that is below 1:

    eps_tide + 1.5 eps_spin  <  0.148 .

So about 0.15 for the tide alone, 0.10 for the centrifugal term alone (whose
coefficient is 1.5 eps_spin), and 0.059 for two equal amplitudes.  ``level_radius``
raises rather than returning nonsense when the surface has opened, which is the
intended way to discover it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

__all__ = ["Deformation"]

def _level_value(spin):
    """The potential at the surface, in the gauge the manual fixes.

    The level value is ours to choose, and it is the star's size.  Choosing it so that
    the angular mean radius is 1 gives a = 1 - eps_spin + O(eps^2), i.e. a level of
    -1/a here.  The gauge matters: with the level pinned at -1 the l = 0 part of the
    centrifugal potential shows up as a uniform change of radius, and every first-order
    quantity then carries a term that the ellipsoid and Lambda immediately discard.
    """
    return -1.0 / (1.0 - spin)


def _p2(cosine):
    """The second Legendre polynomial, 3c^2/2 - 1/2."""
    return 1.5 * cosine ** 2 - 0.5


def _p2_derivative(cosine):
    return 3.0 * cosine


@dataclass(frozen=True)
class Deformation:
    """A monopole plus quadrupoles, in units GM = 1 and reference radius 1.

    ``tide`` is the amplitude of a quadrupole along ``tide_direction``, which a
    companion would raise; ``centrifugal`` is that of the spin quadrupole along
    ``spin_direction``, with its l = 0 piece carried separately.
    """

    tide: float = 0.0
    tide_direction: tuple = (1.0, 0.0, 0.0)
    centrifugal: float = 0.0
    spin_direction: tuple = (0.0, 0.0, 1.0)

    def level_value(self):
        """The potential at the surface, from the gauge of section 8.2."""
        return _level_value(self.centrifugal)

    def axes(self):
        """The two quadrupole directions, normalised."""
        tide = np.asarray(self.tide_direction, dtype=float)
        spin = np.asarray(self.spin_direction, dtype=float)
        return tide / np.linalg.norm(tide), spin / np.linalg.norm(spin)

    def potential(self, x):
        """Phi at x, dimensionless.  x is (..., 3); the last axis is the position."""
        x = np.asarray(x, dtype=float)
        radius = np.linalg.norm(x, axis=-1)
        safe = np.where(radius == 0.0, 1.0, radius)
        direction = x / safe[..., None]
        tide_axis, spin_axis = self.axes()

        total = -1.0 / safe
        if self.tide:
            total = total - self.tide * safe ** 2 * _p2(direction @ tide_axis)
        if self.centrifugal:
            total = total + self.centrifugal * safe ** 2 * (
                _p2(direction @ spin_axis) - 1.0)
        return total

    def gradient(self, x):
        """grad Phi at x, dimensionless.

        Each quadrupole term contributes 2 r P_2 n + r P_2'(c) (e - c n), which is
        the radial part plus the transverse part that tilts the normal.
        """
        x = np.asarray(x, dtype=float)
        radius = np.linalg.norm(x, axis=-1)
        safe = np.where(radius == 0.0, 1.0, radius)
        direction = x / safe[..., None]
        tide_axis, spin_axis = self.axes()

        total = direction / safe[..., None] ** 2
        if self.tide:
            cosines = direction @ tide_axis
            total = total - self.tide * (
                2.0 * safe[..., None] * _p2(cosines)[..., None] * direction
                + safe[..., None] * _p2_derivative(cosines)[..., None]
                * (tide_axis - cosines[..., None] * direction))
        if self.centrifugal:
            cosines = direction @ spin_axis
            total = total + self.centrifugal * (
                2.0 * safe[..., None] * (_p2(cosines) - 1.0)[..., None] * direction
                + safe[..., None] * _p2_derivative(cosines)[..., None]
                * (spin_axis - cosines[..., None] * direction))
        return total

    def level_radius(self, direction):
        """The radius of the level surface along a unit direction.

        The potential rises to a maximum and then falls, because the monopole's
        gradient is outward while the quadrupoles' is not, so a ray can cross the
        level set twice.  The star is the *inner* crossing, and the outer one belongs
        to a structure the quadrupole expansion has no business describing.  So the
        bracket runs from a small radius up to the turning point, where the ray is
        monotone and there is exactly one crossing -- or none, if the surface has
        opened.
        """
        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)

        level = self.level_value()

        def residual(radius):
            return float(self.potential(radius * direction)) - level

        def slope(radius):
            return float(np.dot(self.gradient(radius * direction), direction))

        inner, outer = 0.2, 3.0
        if slope(inner) * slope(outer) < 0.0:
            turning = brentq(slope, inner, outer, xtol=1e-13, rtol=1e-14)
        else:
            turning = outer
        if residual(inner) * residual(turning) > 0.0:
            raise ValueError("no level surface along this direction")
        return brentq(residual, inner, turning, xtol=1e-14, rtol=1e-15)

    def surface_gravity(self, direction):
        """|grad Phi| on the level surface along a unit direction."""
        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        radius = self.level_radius(direction)
        return float(np.linalg.norm(self.gradient(radius * direction)))

    def level_radii(self, directions):
        """The level radii along many directions at once, by vectorised bisection.

        Phi is monotone up to its turning point, and the deepest well has the smallest
        one, so every surface radius lies inside ``(1/(2E))^(1/3)`` with ``E`` that
        depth.  Raises if the surface has opened along any direction.
        """
        directions = np.asarray(directions, dtype=float)
        directions = directions / np.linalg.norm(directions, axis=-1, keepdims=True)
        count = directions.shape[0]

        level = self.level_value()
        # The bracket comes from the *depth* of the perturbation, not from the level: the
        # turning point of -1/r - E r^2 is at r^3 = 1/(2E) with E measured from the
        # unperturbed potential.  Using the gauged level here would give E minus its l = 0
        # part, widening the bracket past the outer root and reporting an open surface.
        depths = -1.0 - self.potential(directions)
        deepest = float(np.max(depths))
        upper = 3.0 if deepest <= 0.0 else (1.0 / (2.0 * deepest)) ** (1.0 / 3.0)
        if np.any(self.potential(upper * directions) < level):
            raise ValueError("the surface is open along some direction")

        lower, high = np.full(count, 0.2), np.full(count, upper)
        for _ in range(60):
            middle = 0.5 * (lower + high)
            rising = self.potential(middle[:, None] * directions) < level
            lower = np.where(rising, middle, lower)
            high = np.where(rising, high, middle)
        return 0.5 * (lower + high)

    def axis_lengths(self):
        """The semi-axes along the coordinate directions, in the order x, y, z."""
        return tuple(self.level_radius(axis) for axis in np.eye(3))

    def volume(self, n_mu=48, n_phi=96):
        """The volume enclosed by the level surface, by quadrature over the sphere."""
        nodes, weights = np.polynomial.legendre.leggauss(n_mu)
        phi = 2.0 * np.pi * np.arange(n_phi) / n_phi
        mu = nodes[:, None]
        sine = np.sqrt(1.0 - mu ** 2)
        directions = np.stack([sine * np.cos(phi)[None, :],
                               sine * np.sin(phi)[None, :],
                               np.broadcast_to(mu, (n_mu, n_phi))], axis=-1)
        flat = directions.reshape(-1, 3)
        radii = np.array([self.level_radius(direction) for direction in flat])
        radii = radii.reshape(n_mu, n_phi)
        integral = np.sum(weights[:, None] * radii ** 3) * (2.0 * np.pi / n_phi)
        return integral / 3.0

    def equivalent_radius(self, **kwargs):
        """The radius of the sphere of the same volume."""
        return (3.0 * self.volume(**kwargs) / (4.0 * np.pi)) ** (1.0 / 3.0)

    def first_order_radii(self, direction):
        """The first-order predictions (delta r, delta g), for comparison.

        In these units g_0 = 1 and r_0 = 1.  The radius moves by delta r = -phi_1.
        The gravity changes for two reasons: the perturbation's own radial derivative,
        which for a quadrupole term going as r^2 is d phi_1 / d r = 2 phi_1 at the
        reference radius, and the displacement of the surface, which moves the
        monopole's 1/r^2 and contributes -2 delta r.  Both are O(eps) and both count.
        """
        direction = np.asarray(direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        tide_axis, spin_axis = self.axes()
        depth = 0.0
        if self.tide:
            depth += self.tide * _p2(float(direction @ tide_axis))
        if self.centrifugal:
            depth -= self.centrifugal * (_p2(float(direction @ spin_axis)) - 1.0)
        # the gauge: the mean of the depth is the spin's l = 0 part, so the
        # displacement is a pure quadrupole
        quadrupole = depth - self.centrifugal
        delta_r = quadrupole
        delta_g = -2.0 * self.centrifugal - 4.0 * quadrupole
        return delta_r, delta_g
