import numpy as np
import pytest

from libsolidstate.lattice.bravais2d import Crystal, primitive_vectors
from libsolidstate.lattice.reciprocal import brillouin_zone, k_path, reciprocal_vectors, special_points
from libsolidstate.lattice.springs import force_constants, neighbour_bonds
from libsolidstate.spectra.dispersion import crystal_dispersion
from libsolidstate.spectra.modes import normal_modes

A = 2.46


def graphene():
    return Crystal(primitive_vectors("hexagonal", A), [[0, 0], [A / 2, A / (2 * np.sqrt(3))]], [12.0, 12.0], ["C", "C"])


def triangular(m=2.0):
    return Crystal(primitive_vectors("hexagonal", A), [[0, 0]], [m], ["A"])


@pytest.mark.parametrize(
    "crystal, ks, count, length",
    [
        (triangular(), [1.0], 3, A),
        (graphene(), [1.0], 3, A / np.sqrt(3)),
        (graphene(), [1.0, 1.0], 3 + 6, None),
    ],
)
def test_neighbour_bonds(crystal, ks, count, length):
    bonds = neighbour_bonds(crystal, ks)
    assert len(bonds) == count
    if length is not None:
        np.testing.assert_allclose([np.linalg.norm(r) for *_, r in bonds], length)


def test_square_shells():
    square = Crystal(primitive_vectors("square", 1.0), [[0, 0]], [1.0], ["A"])
    bonds = neighbour_bonds(square, [1.0, 0.5])
    lengths = sorted(np.linalg.norm(r) for *_, r in bonds)
    np.testing.assert_allclose(lengths, [1, 1, np.sqrt(2), np.sqrt(2)])
    assert sorted(k for _, _, _, k, _ in bonds) == [0.5, 0.5, 1.0, 1.0]


def test_triangular_analytic():
    """Bravais lattice: D(q) = sum over bonds of k/m e e^T 2 (1 - cos q.R)."""
    k, m = 3.0, 2.0
    crystal = triangular(m)
    bonds = neighbour_bonds(crystal, [k])
    qs = np.random.default_rng(0).uniform(-2, 2, size=(20, 2))

    omega = crystal_dispersion(crystal, bonds, qs)
    for q, w in zip(qs, omega):
        D = sum(k / m * np.outer(r, r) / (r @ r) * 2 * (1 - np.cos(q @ r)) for *_, r in bonds)
        np.testing.assert_allclose(w**2, np.linalg.eigvalsh(D), atol=1e-12)


def test_acoustic_modes_vanish_at_gamma():
    omega = crystal_dispersion(graphene(), neighbour_bonds(graphene(), [10.0, 2.0]), [[0.0, 0.0]])
    np.testing.assert_allclose(omega[0, :2], 0, atol=1e-12)
    assert np.all(omega[0, 2:] > 0)


def test_unstable_lattice_gives_negative_frequencies():
    omega = crystal_dispersion(triangular(), neighbour_bonds(triangular(), [-1.0]), [[0.5, 0.0]])
    assert np.all(omega < 0)


def test_matches_finite_periodic_supercell():
    """An N x N periodic supercell has exactly the frequencies at q = (i b1 + j b2) / N."""
    crystal, cells = graphene(), 4
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    n_basis = len(crystal.masses)

    def index(c1, c2, mu):
        return ((c1 % cells) * cells + c2 % cells) * n_basis + mu

    supercell_bonds = [
        (index(c1, c2, mu), index(c1 + n[0], c2 + n[1], nu), k, r)
        for c1 in range(cells)
        for c2 in range(cells)
        for mu, nu, n, k, r in bonds
    ]
    n_atoms = cells * cells * n_basis
    K = force_constants(n_atoms, supercell_bonds, dim=2)
    masses = np.repeat(np.tile(crystal.masses, cells * cells), 2)
    omega_squared, _ = normal_modes(K, masses)

    b1, b2 = reciprocal_vectors(crystal.vectors)
    qs = [(i * b1 + j * b2) / cells for i in range(cells) for j in range(cells)]
    omega = crystal_dispersion(crystal, bonds, qs)
    np.testing.assert_allclose(np.sort(omega.ravel() ** 2), omega_squared, atol=1e-10)


def test_hexagonal_special_points():
    points = special_points("hexagonal", primitive_vectors("hexagonal", A))
    assert np.linalg.norm(points["K"]) == pytest.approx(4 * np.pi / (3 * A))
    assert np.linalg.norm(points["M"]) == pytest.approx(2 * np.pi / (np.sqrt(3) * A))
    corners = brillouin_zone(primitive_vectors("hexagonal", A))
    assert np.min(np.linalg.norm(corners - points["K"], axis=1)) == pytest.approx(0, abs=1e-12)


def test_k_path():
    points = [[0, 0], [1, 0], [1, 1]]
    qs, distance, corners = k_path(points, n=20)
    np.testing.assert_allclose(qs[0], points[0])
    np.testing.assert_allclose(qs[-1], points[-1])
    np.testing.assert_allclose(corners, [0, 1, 2])
    assert np.all(np.diff(distance) > 0)
