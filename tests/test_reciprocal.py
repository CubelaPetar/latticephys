import numpy as np
import pytest

from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.reciprocal import (
    brillouin_zone,
    k_path,
    q_grid,
    reciprocal_vectors,
    special_points,
    wigner_seitz_cell,
)

LATTICES = [
    ("oblique", (2.0, 3.0, 75.0)),
    ("rectangular", (2.0, 3.0)),
    ("centered rectangular", (2.0, 3.0)),
    ("square", (2.0,)),
    ("hexagonal", (2.0,)),
]


def polygon_area(vertices):
    x, y = vertices.T
    return 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1))


@pytest.mark.parametrize("lattice, args", LATTICES)
def test_reciprocal_vectors(lattice, args):
    a = primitive_vectors(lattice, *args)
    np.testing.assert_allclose(a @ reciprocal_vectors(a).T, 2 * np.pi * np.eye(2), atol=1e-12)


@pytest.mark.parametrize("lattice, args", LATTICES)
def test_cell_areas(lattice, args):
    """Wigner-Seitz cell has the primitive cell area, BZ has (2 pi)^2 / area."""
    a = primitive_vectors(lattice, *args)
    area = abs(np.linalg.det(a))
    assert polygon_area(wigner_seitz_cell(a)) == pytest.approx(area)
    assert polygon_area(brillouin_zone(a)) == pytest.approx((2 * np.pi) ** 2 / area)


def test_square_brillouin_zone():
    a = 2.0
    corners = brillouin_zone(primitive_vectors("square", a))
    assert len(corners) == 4
    np.testing.assert_allclose(np.abs(corners), np.pi / a)


def test_hexagonal_brillouin_zone():
    """Hexagon with the K points at distance 4 pi / (3 a)."""
    a = 2.46
    corners = brillouin_zone(primitive_vectors("hexagonal", a))
    assert len(corners) == 6
    np.testing.assert_allclose(np.linalg.norm(corners, axis=1), 4 * np.pi / (3 * a))


def test_hexagonal_special_points():
    a = 2.46
    points = special_points("hexagonal", primitive_vectors("hexagonal", a))
    assert np.linalg.norm(points["K"]) == pytest.approx(4 * np.pi / (3 * a))
    assert np.linalg.norm(points["M"]) == pytest.approx(2 * np.pi / (np.sqrt(3) * a))
    corners = brillouin_zone(primitive_vectors("hexagonal", a))
    assert np.min(np.linalg.norm(corners - points["K"], axis=1)) == pytest.approx(0, abs=1e-12)


def test_k_path():
    points = [[0, 0], [1, 0], [1, 1]]
    qs, distance, corners = k_path(points, n=20)
    np.testing.assert_allclose(qs[0], points[0])
    np.testing.assert_allclose(qs[-1], points[-1])
    np.testing.assert_allclose(corners, [0, 1, 2])
    assert np.all(np.diff(distance) > 0)


def test_reciprocal_vectors_of_a_tilted_layer():
    """Layer (p = 2, d = 3): a_i . b_j = 2 pi delta_ij and every b_j lies in the plane, b_j . (a_1 x a_2) = 0."""
    a = np.array([[1.0, 0.0, 1.0], [0.3, 2.0, 0.0]])
    b = reciprocal_vectors(a)

    np.testing.assert_allclose(a @ b.T, 2 * np.pi * np.eye(2), atol=1e-12)
    np.testing.assert_allclose(b @ np.cross(*a), 0, atol=1e-12)


def test_q_grid_of_a_layer():
    """Layer (p = 2, d = 3): n^2 wave vectors, all in the plane (z = 0), distinct modulo the reciprocal lattice."""
    a = np.array([[2.0, 0.0, 0.0], [1.0, np.sqrt(3), 0.0]])
    qs = q_grid(a, 5)

    assert qs.shape == (25, 3)
    np.testing.assert_allclose(qs[:, 2], 0, atol=1e-12)
    fractional = np.round(qs @ a.T / (2 * np.pi) * 5).astype(int) % 5
    assert len({tuple(f) for f in fractional}) == 25


@pytest.mark.parametrize("lattice, args", [("square", (2.0,)), ("rectangular", (2.0, 3.0)), ("hexagonal", (2.0,))])
def test_special_points_of_a_layer_are_the_2d_ones_in_its_plane(lattice, args):
    """A 2D lattice lifted into 3D keeps its labels and points, with q_z = 0."""
    flat = Crystal(primitive_vectors(lattice, *args), [[0.0, 0.0]], [1.0], ["A"])
    points_2d = special_points(lattice, flat.vectors)
    points_3d = special_points(lattice, embed(flat, 3).vectors)

    assert points_3d.keys() == points_2d.keys()
    for label, q in points_2d.items():
        np.testing.assert_allclose(points_3d[label], [*q, 0.0], atol=1e-12)


@pytest.mark.parametrize("vectors", [np.eye(3), [[2.0, 0.0, 0.0]], [2.0, 1.0]])
def test_special_points_need_two_lattice_vectors(vectors):
    """The 2D labels need exactly 2 lattice vectors (a bulk 3D lattice, a chain or a flat list raise)."""
    with pytest.raises(ValueError, match="2 lattice vectors"):
        special_points("square", vectors)
