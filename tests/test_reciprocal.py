import numpy as np
import pytest

from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.reciprocal import brillouin_zone, reciprocal_vectors, wigner_seitz_cell

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
