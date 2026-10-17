import numpy as np
import pytest

from libsolidstate.lattice.bravais2d import Crystal, atom_positions, primitive_vectors


def lengths_and_angle(vectors):
    a1, a2 = vectors
    cos = a1 @ a2 / (np.linalg.norm(a1) * np.linalg.norm(a2))
    return np.linalg.norm(a1), np.linalg.norm(a2), np.degrees(np.arccos(cos))


@pytest.mark.parametrize(
    "lattice, args, expected",
    [
        ("oblique", (2.0, 3.0, 75.0), (2.0, 3.0, 75.0)),
        ("rectangular", (2.0, 3.0), (2.0, 3.0, 90.0)),
        ("square", (2.0,), (2.0, 2.0, 90.0)),
        ("hexagonal", (2.0,), (2.0, 2.0, 60.0)),
    ],
)
def test_primitive_vectors(lattice, args, expected):
    np.testing.assert_allclose(lengths_and_angle(primitive_vectors(lattice, *args)), expected)


def test_centered_rectangular_has_half_the_conventional_area():
    vectors = primitive_vectors("centered rectangular", 2.0, 3.0)
    assert abs(np.linalg.det(vectors)) == pytest.approx(2.0 * 3.0 / 2)


@pytest.mark.parametrize("lattice, args", [("unknown", (1.0,)), ("rectangular", (1.0,)), ("oblique", (1.0, 2.0))])
def test_invalid_input(lattice, args):
    with pytest.raises(ValueError):
        primitive_vectors(lattice, *args)


def test_graphene_nearest_neighbours():
    """Each carbon atom has 3 neighbours at a / sqrt(3)."""
    a = 2.46
    graphene = Crystal(
        primitive_vectors("hexagonal", a),
        positions=[[0.0, 0.0], [a / 2, a / (2 * np.sqrt(3))]],
        masses=[12.011, 12.011],
        names=["C", "C"],
    )
    positions, _ = atom_positions(graphene, cells=2)

    distances = np.sort(np.linalg.norm(positions - graphene.positions[0], axis=1))
    np.testing.assert_allclose(distances[1:4], a / np.sqrt(3))
    assert distances[4] > a / np.sqrt(3) + 0.1


def test_atom_positions_count():
    crystal = Crystal(primitive_vectors("square", 1.0), [[0, 0], [0.5, 0.5]], [1, 1], ["A", "B"])
    positions, basis_index = atom_positions(crystal, cells=2)

    assert positions.shape == (5 * 5 * 2, 2)
    assert np.bincount(basis_index).tolist() == [25, 25]
