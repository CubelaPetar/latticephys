import numpy as np
import pytest

from libsolidstate.lattice.bonds import neighbour_bonds
from libsolidstate.lattice.structures import honeycomb, kagome, lieb


@pytest.mark.parametrize(
    "structure, atoms, nearest, bonds_per_cell",
    [(honeycomb, 2, 1 / np.sqrt(3), 3), (kagome, 3, 1 / 2, 6), (lieb, 3, 1 / 2, 4)],
)
def test_nearest_neighbours(structure, atoms, nearest, bonds_per_cell):
    """Nearest-neighbour distance (in units of a) and number of nearest-neighbour bonds per cell."""
    a = 2.0
    crystal = structure(a, [1.0] * atoms, ["A"] * atoms)
    bonds = neighbour_bonds(crystal, [1.0])
    assert len(bonds) == bonds_per_cell
    lengths = [np.linalg.norm(crystal.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds]
    np.testing.assert_allclose(lengths, nearest * a)
