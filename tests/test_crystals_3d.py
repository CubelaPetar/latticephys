"""Two-atom 3D crystals (rocksalt, diamond) against closed-form dynamical matrices."""

import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.reciprocal import q_grid, special_points
from latticephys.lattice.structures import diamond, hcp, rocksalt
from latticephys.phonons.dispersion import dispersion


def two_by_two(A, B, m1, m2):
    """Frequencies of D = [[A/m1, B/sqrt(m1 m2)], [B/sqrt(m1 m2), A/m2]], ascending."""
    D = np.array([[A / m1, B / np.sqrt(m1 * m2)], [B / np.sqrt(m1 * m2), A / m2]])
    return np.sqrt(np.clip(np.linalg.eigvalsh(D), 0, None))


def test_rocksalt_along_gamma_x():
    """Rocksalt, central K1 to the 6 nearest neighbours (a/2 along the axes), K2 to the 12 second ones
    (same species, a/2 (+-1, +-1, 0)). For q = (q, 0, 0), c = cos(q a / 2), each polarization is a 2 x 2 problem
    [[A/m1, B/sqrt(m1 m2)], [B/sqrt(m1 m2), A/m2]]:
    longitudinal A = 2 K1 + 4 K2 (1 - c), B = -2 K1 c; transverse (twice) A = 2 K1 + 2 K2 (1 - c), B = -2 K1."""
    a, K1, K2, m1, m2 = 5.64, 1.5, 0.3, 22.99, 35.45
    crystal = rocksalt(a, [m1, m2], ["Na", "Cl"])
    bonds = neighbour_bonds(crystal, [K1, K2])
    assert len(bonds) == 6 + 6 + 6  # per cell: 6 Na-Cl, and 12 neighbours each for Na-Na and Cl-Cl listed once

    for t in np.linspace(0, 1, 7):
        q = 2 * np.pi / a * t
        c = np.cos(q * a / 2)
        longitudinal = two_by_two(2 * K1 + 4 * K2 * (1 - c), -2 * K1 * c, m1, m2)
        transverse = two_by_two(2 * K1 + 2 * K2 * (1 - c), -2 * K1, m1, m2)
        expected = np.sort(np.concatenate([longitudinal, transverse, transverse]))
        np.testing.assert_allclose(dispersion(crystal, bonds, [[q, 0.0, 0.0]])[0], expected, atol=1e-9)


def test_diamond_optical_mode_at_gamma_and_stability():
    """Diamond, central K1 to the 4 nearest neighbours along <111>, K2 to the 12 second ones (same sublattice).
    At Γ the sublattices move rigidly against each other: sum e e^T = 4/3 over the 4 bonds, reduced mass M/2,
    so omega^2 = 8 K1 / (3 M), triply degenerate; K2 does not enter. With K2 > 0 no mode is unstable."""
    a, K1, K2, M = 5.431, 7.0, 1.0, 28.086
    si = diamond(a, [M, M], ["Si", "Si"])
    bonds = neighbour_bonds(si, [K1, K2])

    omega = dispersion(si, bonds, [np.zeros(3)])[0]
    np.testing.assert_allclose(omega, [0, 0, 0] + [np.sqrt(8 * K1 / (3 * M))] * 3, atol=1e-12)
    assert np.all(dispersion(si, bonds, q_grid(si.vectors, 8)) >= 0)


def test_diamond_ta_at_x_depends_only_on_second_neighbours():
    """At X nearest-neighbour central springs exert no transverse force, so TA is the fcc sublattice
    value omega^2 = 4 K2 / M (as for fcc with nearest neighbours), for any K1."""
    a, K2, M = 5.431, 0.6, 28.086
    si = diamond(a, [M, M], ["Si", "Si"])
    x = special_points("fcc", si.vectors)["X"]
    for K1 in (1.0, 10.0, 30.0):
        omega = dispersion(si, neighbour_bonds(si, [K1, K2]), [x])[0]
        np.testing.assert_allclose(omega[:2], np.sqrt(4 * K2 / M), rtol=1e-12)


def test_hcp_optical_modes_at_gamma():
    """hcp, central springs K to the 6 neighbours in the other sublattice (distance d^2 = a^2/3 + c^2/4) and K2 in
    the plane. At Γ the sublattices move rigidly against each other, only the 6 interlayer springs act (reduced
    mass M/2): in plane omega^2 = 2 K a^2 / (M d^2) (twice), out of plane omega^2 = 3 K c^2 / (M d^2); K2 drops out."""
    a, c, K, K2, M = 3.209, 5.211, 1.2, 0.7, 24.305
    mg = hcp(a, c, [M, M], ["Mg", "Mg"])
    bonds = neighbour_bonds(mg, [K, K2])  # out-of-plane neighbours (3.198 Å) come before in-plane (3.209 Å)
    d2 = a**2 / 3 + c**2 / 4

    omega = dispersion(mg, bonds, [np.zeros(3)])[0]
    expected = np.sqrt([0, 0, 0, 2 * K * a**2 / (M * d2), 2 * K * a**2 / (M * d2), 3 * K * c**2 / (M * d2)])
    np.testing.assert_allclose(omega, np.sort(expected), atol=1e-12)
