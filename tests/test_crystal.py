import numpy as np
import pytest

from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal, atom_positions, embed


def test_graphene_nearest_neighbours():
    """Each carbon atom has 3 neighbours at a / sqrt(3)."""
    a = 2.46
    honeycomb = [[0.0, 0.0], [a / 2, a / (2 * np.sqrt(3))]]
    graphene = Crystal(primitive_vectors("hexagonal", a), honeycomb, [12.0] * 2, ["C"] * 2)
    positions, _ = atom_positions(graphene, cells=2)

    distances = np.sort(np.linalg.norm(positions - graphene.positions[0], axis=1))
    np.testing.assert_allclose(distances[1:4], a / np.sqrt(3))
    assert distances[4] > a / np.sqrt(3) + 0.1


def test_atom_positions_count():
    crystal = Crystal(primitive_vectors("square", 1.0), [[0, 0], [0.5, 0.5]], [1, 1], ["A", "B"])
    positions, basis_index = atom_positions(crystal, cells=2)

    assert positions.shape == (5 * 5 * 2, 2)
    assert np.bincount(basis_index).tolist() == [25, 25]


def test_layer_has_dimension_3_and_periodicity_2():
    """Two lattice vectors in 3D space: d = 3 from the positions, p = 2 from the vectors."""
    layer = Crystal([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 0.0]], [1.0], ["A"])
    assert (layer.dim, layer.periodicity) == (3, 2)


def test_1d_crystal_accepts_flat_positions():
    chain = Crystal([[2.0]], [0.0, 1.0], [1.0, 3.0], ["A", "B"])
    assert chain.dim == 1
    assert chain.positions.shape == (2, 1)
    np.testing.assert_allclose(chain.bond_vector(1, 0, 1), [1.0])


@pytest.mark.parametrize(
    "vectors, positions, masses, names",
    [
        (np.eye(3)[:, :2], [[0.0, 0.0]], [1.0], ["A"]),  # 3 lattice vectors in 2D
        ([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], [1.0], ["A"]),  # dependent vectors
        (np.eye(2), [[0.0, 0.0, 0.0]], [1.0], ["A"]),  # positions in 3D, lattice in 2D
        (np.eye(2), [[0.0, 0.0]], [1.0], ["A", "B"]),  # names do not match atoms
        (np.eye(2), [[0.0, 0.0], [0.5, 0.5]], [1.0], ["A", "B"]),  # masses do not match atoms
        (np.eye(2), np.zeros((2, 2, 2)), [1.0, 1.0], ["A", "B"]),  # positions not (B, d)
        (np.eye(2), [0.0, 0.0], [1.0], ["A"]),  # flat positions outside 1D
    ],
)
def test_invalid_crystal(vectors, positions, masses, names):
    with pytest.raises(ValueError):
        Crystal(vectors, positions, masses, names)


def test_atom_count_comes_from_positions():
    """Two atoms in 2D with one mass: the positions say B = 2, so the mass count is the error."""
    with pytest.raises(ValueError, match="mass"):
        Crystal(np.eye(2), [[0.0, 0.0], [0.5, 0.5]], [1.0], ["A", "B"])


@pytest.mark.parametrize("n", [(1, 0, 0), 1, (1,)])
def test_cell_index_of_wrong_length_raises(n):
    """A layer (p = 2) needs cell indices with 2 entries; a third (or missing) entry is a typo."""
    layer = Crystal([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], [[0.0, 0.0, 0.0]], [1.0], ["A"])
    with pytest.raises(ValueError, match="cell index"):
        layer.bond_vector(0, 0, n)


def test_embed_lifts_a_chain_into_the_plane():
    """embed appends zero components: the chain keeps its geometry with y = 0, p = 1, d = 2."""
    chain = Crystal([[2.0]], [0.0, 1.0], [1.0, 3.0], ["A", "B"])
    planar = embed(chain, 2)

    np.testing.assert_array_equal(planar.vectors, [[2.0, 0.0]])
    np.testing.assert_array_equal(planar.positions, [[0.0, 0.0], [1.0, 0.0]])
    assert (planar.dim, planar.periodicity) == (2, 1)
    with pytest.raises(ValueError):
        embed(planar, 1)
