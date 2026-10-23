import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.structures import diamond, hcp, honeycomb, kagome, lieb, rocksalt


@pytest.mark.parametrize(
    "structure, atoms, nearest, bonds_per_cell",
    [(honeycomb, 2, 1 / np.sqrt(3), 3), (kagome, 3, 1 / 2, 6), (lieb, 3, 1 / 2, 4),
     (rocksalt, 2, 1 / 2, 6), (diamond, 2, np.sqrt(3) / 4, 4)],
)
def test_nearest_neighbours(structure, atoms, nearest, bonds_per_cell):
    """Nearest-neighbour distance (in units of a) and number of nearest-neighbour bonds per cell."""
    a = 2.0
    crystal = structure(a, [1.0] * atoms, ["A"] * atoms)
    bonds = neighbour_bonds(crystal, [1.0])
    assert len(bonds) == bonds_per_cell
    lengths = [np.linalg.norm(crystal.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds]
    np.testing.assert_allclose(lengths, nearest * a)


def test_hcp_neighbours():
    """hcp with ideal c/a = sqrt(8/3): all 12 neighbours of each atom at distance a (12 bonds per cell).
    With a larger c/a the 6 out-of-plane neighbours at sqrt(a^2/3 + c^2/4) move out: 6 in-plane bonds per cell."""
    a = 3.0
    ideal = hcp(a, a * np.sqrt(8 / 3), [1.0] * 2, ["A"] * 2)
    bonds = neighbour_bonds(ideal, [1.0])
    assert len(bonds) == 12
    np.testing.assert_allclose([np.linalg.norm(ideal.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds], a)

    tall = hcp(a, 1.9 * a, [1.0] * 2, ["A"] * 2)
    assert len(neighbour_bonds(tall, [1.0])) == 6
