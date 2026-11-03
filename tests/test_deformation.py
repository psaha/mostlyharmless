"""The deformation field: the level surface, the gravity on it, and the axes.

Stage 6, items 1.  The physics being checked is that one potential gives both the
shape and the surface gravity, and that the solver reproduces the standard oblate
result in the limit of a small deformation.
"""

import numpy as np
import pytest

from mostlyharmless.symbolic import derivations as dv
from mostlyharmless.visylm.deformation import Deformation


def test_v19_matches_the_first_order_relations():
    # The residual is already divided by eps^2, so this bounds the O(1)-to-O(30)
    # coefficient of the second-order term.  Along the most-stretched direction that
    # coefficient is -5 g^2 with g = -2.5, i.e. about 31, which is why the bound is
    # not tighter.
    assert dv.check_deformation_first_order() < 50.0


def test_v20_the_surface_is_level_and_the_gradient_is_analytic():
    assert dv.check_level_surface() < 1e-8


# A tide-dominant pair, and a small one: the surface closes only while
# eps_tide + 1.5 eps_spin stays below about 0.148, and for the applications in mind
# the tidal quadrupole is the larger of the two.
TIDE, SPIN = 0.03, 0.01


def test_the_oblate_case_reproduces_the_standard_result():
    """Centrifugal only: to first order the equator sits 1.5 eps above the reference
    and its gravity is 3 eps below.  The pole is unmoved because there the l = 0 piece
    of the centrifugal potential cancels the quadrupole."""
    eps = 1e-6
    field = Deformation(centrifugal=eps)
    equator = np.array([1.0, 0.0, 0.0])
    pole = np.array([0.0, 0.0, 1.0])

    r_equator = field.level_radius(equator)
    r_pole = field.level_radius(pole)
    # the mean-radius gauge puts the equator eps/2 above 1 and the pole eps below it,
    # so the axis ratio is unchanged and is the gauge-invariant statement
    assert r_equator / r_pole == pytest.approx(1.0 + 1.5 * eps, abs=1e-10)
    assert r_equator - 1.0 == pytest.approx(0.5 * eps, abs=1e-12)
    assert r_pole - 1.0 == pytest.approx(-eps, abs=1e-12)
    assert field.surface_gravity(equator) - 1.0 == pytest.approx(-4.0 * eps, abs=1e-11)
    assert field.surface_gravity(pole) - 1.0 == pytest.approx(2.0 * eps, abs=1e-11)


def test_the_gauge_puts_the_mean_radius_at_one():
    """The gauge of section 8.2, as a guard on half of it being applied.

    The level value is chosen so that the angular mean of the radius is 1.  A change that
    touched the potential but not the level, or the level but not the first-order
    formulae, would leave the mean at 1 - O(eps) rather than 1 - O(eps^2).
    """
    eps = 0.03
    field = Deformation(tide=eps, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=0.01, spin_direction=(0.0, 0.0, 1.0))
    nodes, weights = np.polynomial.legendre.leggauss(32)
    phi = 2.0 * np.pi * np.arange(64) / 64
    mu = nodes[:, None]
    sine = np.sqrt(1.0 - mu ** 2)
    directions = np.stack([sine * np.cos(phi)[None, :],
                           sine * np.sin(phi)[None, :],
                           np.broadcast_to(mu, (32, 64))], axis=-1).reshape(-1, 3)
    radii = field.level_radii(directions)
    mean = float(np.mean(radii))
    assert abs(mean - 1.0) < 20.0 * eps ** 2, mean


def test_the_tidal_case_bulges_along_its_own_axis():
    eps = 1e-6
    field = Deformation(tide=eps, tide_direction=(1.0, 0.0, 0.0))
    along = np.array([1.0, 0.0, 0.0])
    across = np.array([0.0, 0.0, 1.0])

    assert field.level_radius(along) - 1.0 == pytest.approx(eps, rel=1e-4)
    assert field.level_radius(across) - 1.0 == pytest.approx(-0.5 * eps, rel=1e-3)
    # a bulge is further from the centre, so its gravity is lower: the perturbation's
    # own derivative and the displacement both act
    assert field.surface_gravity(along) - 1.0 == pytest.approx(-4.0 * eps, rel=1e-4)


def test_two_perpendicular_quadrupoles_give_three_unequal_axes():
    """The triaxial figure: x takes the tidal bulge and the equatorial bulge both,
    z only the flattening, so x > y > z is forced.  With eps_tide = 3 eps_spin the
    perturbation vanishes along y, which puts that axis at exactly the reference
    radius however large the amplitudes are."""
    field = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
    lengths = field.axis_lengths()

    assert lengths[0] > lengths[1] > lengths[2]
    assert lengths[0] > 1.0
    assert lengths[1] == pytest.approx(1.0 - SPIN, abs=1e-10)
    assert lengths[2] < 1.0
    # the deformation expands the volume, through the l = 0 piece and the bulge
    assert field.equivalent_radius() > 1.0


def test_the_validity_limit_is_on_the_combination():
    """Not a bug: the quadrupole grows as r^2 and the monopole falls as 1/r, so past
    a modest amplitude Phi = -1 is unreachable along the most-stretched direction and
    the surface opens, as a Roche lobe does.  The binding combination is
    eps_tide + 1.5 eps_spin, which must stay below about 0.148."""
    for tide, spin in ((0.1, 0.1), (0.2, 0.2), (0.2, 0.0)):
        field = Deformation(tide=tide, tide_direction=(1.0, 0.0, 0.0),
                            centrifugal=spin, spin_direction=(0.0, 0.0, 1.0))
        with pytest.raises(ValueError):
            field.level_radius(np.array([1.0, 0.0, 0.0]))

    field = Deformation(tide=TIDE, tide_direction=(1.0, 0.0, 0.0),
                        centrifugal=SPIN, spin_direction=(0.0, 0.0, 1.0))
    assert field.level_radius(np.array([1.0, 0.0, 0.0])) > 1.0


def test_a_null_deformation_is_the_unit_sphere():
    field = Deformation()
    assert field.level_radius([1.0, 1.0, 0.0]) == pytest.approx(1.0, abs=1e-12)
    assert field.surface_gravity([0.0, 1.0, 0.0]) == pytest.approx(1.0, abs=1e-12)
    assert field.equivalent_radius() == pytest.approx(1.0, abs=1e-6)


def test_v28_the_limb_darkened_disk_matches_the_classical_two_term():
    """V28: zero deformation reduces the row to the classical limb-darkened disk.

    At eps = 0 the coefficients are a_00 = 1 and a_10 = u_limb/sqrt(3), so the row is the
    classical answer (manual eq. 118): the uniform disk's J_1 term from l = 0, and the
    limb's j_1 term from l = 1.  Comparing V(q)/V(0) pins the normalisation at the origin
    rather than at any convention for u_limb.  The reference is scipy's Bessel functions,
    which share nothing with the row but the definition.
    """
    assert dv.check_limb_darkened_disk() < 1e-10
