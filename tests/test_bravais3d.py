import numpy as np
import pytest

from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import cell_offsets
from latticephys.lattice.reciprocal import reciprocal_vectors, special_points

A, B, C = 2.0, 2.6, 3.3  # a < b < c
ALPHA, BETA, GAMMA = 75.0, 80.0, 85.0


def conventional(a, b, c, alpha=90.0):
    """Conventional cell vectors: a along x, b along y, c in the yz plane at alpha to b."""
    t = np.radians(alpha)
    return np.array([[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, c * np.cos(t), c * np.sin(t)]])


# lattice, parameters, conventional cell, centring vectors as fractions of the conventional cell
CENTRED = [
    ("cubic", dict(a=A), conventional(A, A, A), []),
    ("fcc", dict(a=A), conventional(A, A, A), [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
    ("bcc", dict(a=A), conventional(A, A, A), [[0.5, 0.5, 0.5]]),
    ("tetragonal", dict(a=A, c=C), conventional(A, A, C), []),
    ("body-centered tetragonal", dict(a=A, c=C), conventional(A, A, C), [[0.5, 0.5, 0.5]]),
    ("orthorhombic", dict(a=A, b=B, c=C), conventional(A, B, C), []),
    ("face-centered orthorhombic", dict(a=A, b=B, c=C), conventional(A, B, C),
     [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
    ("body-centered orthorhombic", dict(a=A, b=B, c=C), conventional(A, B, C), [[0.5, 0.5, 0.5]]),
    ("base-centered orthorhombic", dict(a=A, b=B, c=C), conventional(A, B, C), [[0.5, 0.5, 0]]),
    ("monoclinic", dict(a=A, b=B, c=C, alpha=ALPHA), conventional(A, B, C, ALPHA), []),
    ("base-centered monoclinic", dict(a=A, b=B, c=C, alpha=ALPHA), conventional(A, B, C, ALPHA), [[0.5, 0.5, 0]]),
]


@pytest.mark.parametrize("lattice, parameters, cell, centrings", CENTRED, ids=[c[0] for c in CENTRED])
def test_primitive_cell_generates_the_conventional_lattice(lattice, parameters, cell, centrings):
    """Volume of the primitive cell = conventional volume / lattice points per cell (1 + centrings), and the
    conventional vectors and centring vectors are integer combinations of the primitive vectors: together
    the primitive vectors generate exactly the conventional lattice."""
    vectors = primitive_vectors(lattice, **parameters)
    assert abs(np.linalg.det(vectors)) == pytest.approx(abs(np.linalg.det(cell)) / (1 + len(centrings)), rel=1e-12)

    for v in np.vstack([cell, np.reshape(centrings, (-1, 3)) @ cell]):
        n = np.linalg.solve(vectors.T, v)
        np.testing.assert_allclose(n, np.round(n), atol=1e-12)


def lengths_and_angles(vectors):
    a, b, c = np.linalg.norm(vectors, axis=1)
    angle = lambda u, v: np.degrees(np.arccos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))  # noqa: E731
    return (a, b, c), (angle(vectors[1], vectors[2]), angle(vectors[0], vectors[2]), angle(vectors[0], vectors[1]))


@pytest.mark.parametrize(
    "lattice, parameters, expected_lengths, expected_angles",
    [
        ("hexagonal", dict(a=A, c=C), (A, A, C), (90, 90, 60)),
        ("rhombohedral", dict(a=A, alpha=ALPHA), (A, A, A), (ALPHA, ALPHA, ALPHA)),
        ("triclinic", dict(a=A, b=B, c=C, alpha=ALPHA, beta=BETA, gamma=GAMMA), (A, B, C), (ALPHA, BETA, GAMMA)),
    ],
)
def test_primitive_lattices_have_their_lengths_and_angles(lattice, parameters, expected_lengths, expected_angles):
    """Primitive (1 point per cell) lattices: |a_i| and the angles (alpha between a2, a3; beta a1, a3; gamma a1, a2)."""
    lengths, angles = lengths_and_angles(primitive_vectors(lattice, **parameters))
    np.testing.assert_allclose(lengths, expected_lengths, rtol=1e-12)
    np.testing.assert_allclose(angles, expected_angles, rtol=1e-12)


def volume(a, b, c, alpha, beta, gamma):
    """Cell volume from lengths and angles: a b c sqrt(1 + 2 cos a cos b cos g - cos^2 a - cos^2 b - cos^2 g)."""
    ca, cb, cg = np.cos(np.radians([alpha, beta, gamma]))
    return a * b * c * np.sqrt(1 + 2 * ca * cb * cg - ca**2 - cb**2 - cg**2)


def test_volumes_of_the_primitive_lattices():
    """Hexagonal sqrt(3)/2 a^2 c; triclinic from its lengths and angles; rhombohedral one third of its
    hexagonal setting, a_H = 2 a sin(alpha/2), c_H = a sqrt(3 (1 + 2 cos alpha)), 3 points per hexagonal cell."""
    det = lambda lattice, **p: abs(np.linalg.det(primitive_vectors(lattice, **p)))  # noqa: E731
    assert det("hexagonal", a=A, c=C) == pytest.approx(np.sqrt(3) / 2 * A**2 * C, rel=1e-12)
    assert det("triclinic", a=A, b=B, c=C, alpha=ALPHA, beta=BETA, gamma=GAMMA) == pytest.approx(
        volume(A, B, C, ALPHA, BETA, GAMMA), rel=1e-12)
    a_h, c_h = 2 * A * np.sin(np.radians(ALPHA) / 2), A * np.sqrt(3 * (1 + 2 * np.cos(np.radians(ALPHA))))
    assert det("rhombohedral", a=A, alpha=ALPHA) == pytest.approx(np.sqrt(3) / 2 * a_h**2 * c_h / 3, rel=1e-12)


@pytest.mark.parametrize(
    "lattice, parameters",
    [
        ("tetragonal", dict(a=A)),  # missing c
        ("hexagonal", dict(a=A)),
        ("orthorhombic", dict(a=A, c=C)),
        ("rhombohedral", dict(a=A)),
        ("monoclinic", dict(a=A, b=B, c=C)),
        ("base-centered monoclinic", dict(a=A, c=C, alpha=ALPHA)),
        ("triclinic", dict(a=A, b=B, c=C, alpha=ALPHA, beta=BETA)),
        ("diamond", dict(a=A)),  # unknown lattice
        ("cubic", dict(a=A, c=C)),  # c is not a parameter of cubic
        ("rhombohedral", dict(a=A, alpha=130.0)),  # alpha >= 120° has no cell
        ("triclinic", dict(a=A, b=B, c=C, alpha=10.0, beta=10.0, gamma=100.0)),  # angles cannot close a cell
    ],
)
def test_missing_unused_or_impossible_parameters_raise(lattice, parameters):
    """Every input that would give a wrong or NaN cell raises instead."""
    with pytest.raises(ValueError):
        primitive_vectors(lattice, **parameters)


LABELS = {
    "cubic": "Γ X M R",
    "fcc": "Γ X W K L U",
    "bcc": "Γ H N P",
    "tetragonal": "Γ X M Z R A",
    "hexagonal": "Γ M K A L H",
}
NEEDS_C = {"tetragonal", "hexagonal"}


@pytest.mark.parametrize("c", [C, 0.6 * A], ids=["c > a", "c < a"])
@pytest.mark.parametrize("lattice", LABELS)
def test_special_points_lie_on_the_zone_boundary(lattice, c):
    """a_i . b_j = 2 pi delta_ij; every labelled point q except Γ is on the Brillouin-zone boundary:
    no reciprocal lattice point G is closer than Γ, |q - G| >= |q|, with equality for some G != 0."""
    vectors = primitive_vectors(lattice, A, c=c if lattice in NEEDS_C else None)
    b = reciprocal_vectors(vectors)
    np.testing.assert_allclose(vectors @ b.T, 2 * np.pi * np.eye(3), atol=1e-12)

    points = special_points(lattice, vectors)
    assert set(points) == set(LABELS[lattice].split())
    G = (cell_offsets(5, 3) - 2) @ b
    G = G[np.any(G != 0, axis=1)]
    for label, q in points.items():
        if label == "Γ":
            continue
        distance_to_G = np.linalg.norm(q - G, axis=1)
        assert distance_to_G.min() == pytest.approx(np.linalg.norm(q), rel=1e-12), label


def test_fcc_and_hexagonal_points_in_cartesian_form():
    """fcc: X = 2 pi/a (1, 0, 0) up to symmetry, |L| = sqrt(3) pi/a; hexagonal: |K| = 4 pi / 3a, A = (0, 0, pi/c)."""
    fcc = special_points("fcc", primitive_vectors("fcc", A))
    assert np.linalg.norm(fcc["X"]) == pytest.approx(2 * np.pi / A)
    assert np.linalg.norm(fcc["L"]) == pytest.approx(np.sqrt(3) * np.pi / A)
    hexagonal = special_points("hexagonal", primitive_vectors("hexagonal", A, c=C))
    assert np.linalg.norm(hexagonal["K"]) == pytest.approx(4 * np.pi / (3 * A))
    np.testing.assert_allclose(hexagonal["A"], [0, 0, np.pi / C], atol=1e-12)


@pytest.mark.parametrize(
    "lattice, parameters",
    [("orthorhombic", dict(b=B, c=C)), ("rhombohedral", dict(alpha=ALPHA)),
     ("triclinic", dict(b=B, c=C, alpha=ALPHA, beta=BETA, gamma=GAMMA))],
)
def test_special_points_raise_for_lattices_without_fixed_labels(lattice, parameters):
    """Their labels depend on the parameters (Setyawan-Curtarolo variants); special_points raises."""
    vectors = primitive_vectors(lattice, A, **parameters)
    with pytest.raises(ValueError, match="no standard labels"):
        special_points(lattice, vectors)
