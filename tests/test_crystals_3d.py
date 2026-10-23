"""Two-atom 3D crystals (rocksalt, diamond) against closed-form dynamical matrices."""

import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.structures import rocksalt
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
