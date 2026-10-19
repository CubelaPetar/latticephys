import numpy as np
import pytest

from latticephys.lattice.bravais2d import primitive_vectors


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

