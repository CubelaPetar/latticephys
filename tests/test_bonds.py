import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.structures import honeycomb

A = 2.46


def graphene():
    return honeycomb(A, [12.0, 12.0], ["C", "C"])


def triangular(m=2.0):
    return Crystal(primitive_vectors("hexagonal", A), [[0, 0]], [m], ["A"])


@pytest.mark.parametrize(
    "crystal, spring_constants, count, length",
    [
        (triangular(), [1.0], 3, A),
        (graphene(), [1.0], 3, A / np.sqrt(3)),
        (graphene(), [1.0, 1.0], 3 + 6, None),
    ],
)
def test_neighbour_bonds(crystal, spring_constants, count, length):
    bonds = neighbour_bonds(crystal, spring_constants)
    assert len(bonds) == count
    if length is not None:
        lengths = [np.linalg.norm(crystal.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds]
        np.testing.assert_allclose(lengths, length)


def test_square_shells():
    square = Crystal(primitive_vectors("square", 1.0), [[0, 0]], [1.0], ["A"])
    bonds = neighbour_bonds(square, [1.0, 0.5])
    lengths = sorted(np.linalg.norm(square.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds)
    np.testing.assert_allclose(lengths, [1, 1, np.sqrt(2), np.sqrt(2)])
    assert sorted(K for *_, K in bonds) == [0.5, 0.5, 1.0, 1.0]


def test_neighbour_bonds_raise_when_search_too_small():
    """In a 1 x 10 cell the 4th shell (distance 4) lies beyond search = 3; the 10 Å bond must not take its place."""
    cell = Crystal([[1.0, 0.0], [0.0, 10.0]], [[0.0, 0.0]], [1.0], ["A"])
    with pytest.raises(ValueError):
        neighbour_bonds(cell, [1.0] * 4)

    bonds = neighbour_bonds(cell, [1.0] * 4, search=4)
    lengths = {round(float(np.linalg.norm(cell.bond_vector(mu, nu, n))), 6) for mu, nu, n, _ in bonds}
    assert lengths == {1.0, 2.0, 3.0, 4.0}


def test_neighbour_bonds_raise_for_missing_shells():
    chain = Crystal([[1.0]], [0.0], [1.0], ["A"])
    with pytest.raises(ValueError):
        neighbour_bonds(chain, [1.0] * 5)


def test_neighbour_bonds_do_not_depend_on_search():
    crystal = Crystal(primitive_vectors("oblique", 2.5, 3.1, 73.0), [[0.0, 0.0], [0.9, 1.3]], [1.0] * 2, ["A"] * 2)
    spring_constants = [3.0, 2.0, 1.0]
    wider = neighbour_bonds(crystal, spring_constants, search=6)
    assert sorted(neighbour_bonds(crystal, spring_constants)) == sorted(wider)


def test_supercell_layout():
    crystal = graphene()
    supercell, _ = periodic_supercell(crystal, [], 3)
    assert supercell.positions.shape == (3 * 3 * 2, 2)
    np.testing.assert_allclose(supercell.vectors, 3 * crystal.vectors)
    expected = crystal.vectors[0] + 2 * crystal.vectors[1] + crystal.positions[1]  # atom 1 of cell (1, 2)
    np.testing.assert_allclose(supercell.positions[(1 * 3 + 2) * 2 + 1], expected)
    np.testing.assert_allclose(supercell.masses, 12.0)


def test_supercell_keeps_bond_vectors():
    """Wrapping into the supercell must not change any bond vector."""
    crystal = graphene()
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, 3)

    original = [crystal.bond_vector(mu, nu, n) for mu, nu, n, _ in bonds] * 9
    wrapped = [supercell.bond_vector(i, j, n) for i, j, n, _ in supercell_bonds]
    np.testing.assert_allclose(wrapped, original, atol=1e-12)
