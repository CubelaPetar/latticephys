import numpy as np
import pytest

from libsolidstate.dynamics.integrate import simulate
from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.reciprocal import brillouin_zone, k_path, reciprocal_vectors, special_points
from libsolidstate.lattice.springs import force_constants, neighbour_bonds, periodic_supercell
from libsolidstate.spectra.dispersion import bloch_mode, dispersion
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
        lengths = [np.linalg.norm(crystal.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds]
        np.testing.assert_allclose(lengths, length)


def test_square_shells():
    square = Crystal(primitive_vectors("square", 1.0), [[0, 0]], [1.0], ["A"])
    bonds = neighbour_bonds(square, [1.0, 0.5])
    lengths = sorted(np.linalg.norm(square.bond_vector(mu, nu, n)) for mu, nu, n, _ in bonds)
    np.testing.assert_allclose(lengths, [1, 1, np.sqrt(2), np.sqrt(2)])
    assert sorted(k for *_, k in bonds) == [0.5, 0.5, 1.0, 1.0]


def test_triangular_analytic():
    """Bravais lattice: D(q) = sum over bonds of k/m e e^T 2 (1 - cos q.R)."""
    k, m = 3.0, 2.0
    crystal = triangular(m)
    bonds = neighbour_bonds(crystal, [k])
    qs = np.random.default_rng(0).uniform(-2, 2, size=(20, 2))

    omega = dispersion(crystal, bonds, qs)
    vectors = [crystal.bond_vector(mu, nu, n) for mu, nu, n, _ in bonds]
    for q, w in zip(qs, omega):
        D = sum(k / m * np.outer(r, r) / (r @ r) * 2 * (1 - np.cos(q @ r)) for r in vectors)
        np.testing.assert_allclose(w**2, np.linalg.eigvalsh(D), atol=1e-12)


def test_acoustic_modes_vanish_at_gamma():
    omega = dispersion(graphene(), neighbour_bonds(graphene(), [10.0, 2.0]), [[0.0, 0.0]])
    np.testing.assert_allclose(omega[0, :2], 0, atol=1e-12)
    assert np.all(omega[0, 2:] > 0)


def test_unstable_lattice_gives_negative_frequencies():
    omega = dispersion(triangular(), neighbour_bonds(triangular(), [-1.0]), [[0.5, 0.0]])
    assert np.all(omega < 0)


def test_matches_finite_periodic_supercell():
    """An N x N periodic supercell has exactly the frequencies at q = (i b1 + j b2) / N."""
    crystal, cells = graphene(), 4
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, cells)
    omega_squared, _ = normal_modes(force_constants(supercell, supercell_bonds), np.repeat(supercell.masses, 2))

    b1, b2 = reciprocal_vectors(crystal.vectors)
    qs = [(i * b1 + j * b2) / cells for i in range(cells) for j in range(cells)]
    omega = dispersion(crystal, bonds, qs)
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


def test_supercell_layout():
    crystal = graphene()
    supercell, _ = periodic_supercell(crystal, [], 3)
    assert supercell.positions.shape == (3 * 3 * 2, 2)
    np.testing.assert_allclose(supercell.vectors, 3 * crystal.vectors)
    expected = crystal.vectors[0] + 2 * crystal.vectors[1] + crystal.positions[1]  # atom 1 of cell (1, 2)
    np.testing.assert_allclose(supercell.positions[(1 * 3 + 2) * 2 + 1], expected)
    np.testing.assert_allclose(supercell.masses, 12.0)


def test_bloch_mode_oscillates_in_supercell():
    """Optical Bloch mode at M, set up in a supercell, evolves as u0 cos(omega t)."""
    crystal, cells = graphene(), 4
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, cells)
    K = force_constants(supercell, supercell_bonds)

    q = reciprocal_vectors(crystal.vectors)[0] / 2
    omega, u0 = bloch_mode(crystal, bonds, q, branch=3, cells=cells)

    dt, steps = 0.001, 3000
    u, _ = simulate(K, np.repeat(supercell.masses, 2), u0, np.zeros_like(u0), dt, steps)
    t = dt * np.arange(steps + 1)
    np.testing.assert_allclose(u, np.outer(np.cos(omega * t), u0), atol=1e-6)


def test_bloch_mode_rejects_q_not_fitting_supercell():
    crystal = graphene()
    q = reciprocal_vectors(crystal.vectors)[0] / 3
    with pytest.raises(ValueError):
        bloch_mode(crystal, neighbour_bonds(crystal, [1.0]), q, branch=0, cells=4)


def test_supercell_keeps_bond_vectors():
    """Wrapping into the supercell must not change any bond vector."""
    crystal = graphene()
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, 3)

    original = [crystal.bond_vector(mu, nu, n) for mu, nu, n, _ in bonds] * 9
    wrapped = [supercell.bond_vector(i, j, n) for i, j, n, _ in supercell_bonds]
    np.testing.assert_allclose(wrapped, original, atol=1e-12)


@pytest.mark.parametrize("point", ["Γ", "M", "K"])
@pytest.mark.parametrize("branch", range(4))
def test_bloch_mode_start_is_never_zero(point, branch):
    """The largest polarization component is made real, so max |u0| >= 1 / sqrt(B d m_max)."""
    crystal = graphene()
    q = special_points("hexagonal", crystal.vectors)[point]
    _, u0 = bloch_mode(crystal, neighbour_bonds(crystal, [10.0, 2.0]), q, branch=branch, cells=6)
    assert np.abs(u0).max() >= 1 / np.sqrt(4 * crystal.masses.max()) - 1e-12


@pytest.mark.parametrize("branch", [0, 1])
def test_bloch_mode_acoustic_at_gamma_is_exactly_zero(branch):
    omega, _ = bloch_mode(graphene(), neighbour_bonds(graphene(), [10.0, 2.0]), [0.0, 0.0], branch=branch, cells=6)
    assert omega == 0.0


@pytest.mark.parametrize(
    "vectors, ks",
    [
        (primitive_vectors("square", 2.5), [10.0, 3.0]),
        (primitive_vectors("oblique", 2.5, 3.1, 73.0), [10.0, 3.0]),
        (primitive_vectors("hexagonal", 2.46), [10.0, 4.0, 1.0]),
    ],
)
def test_acoustic_modes_at_gamma_alone_are_exactly_zero(vectors, ks):
    """Round-off at Gamma is judged against the spring scale, also when Gamma is the only q-point."""
    crystal = Crystal(vectors, [[0.0, 0.0]], [12.011], ["C"])
    bonds = neighbour_bonds(crystal, ks)
    np.testing.assert_array_equal(dispersion(crystal, bonds, [[0.0, 0.0]]), 0.0)
    assert bloch_mode(crystal, bonds, [0.0, 0.0], branch=0, cells=2)[0] == 0.0


def kagome():
    vectors = primitive_vectors("hexagonal", 4.0)
    return Crystal(vectors, [[0.0, 0.0], vectors[0] / 2, vectors[1] / 2], [12.0] * 3, ["A"] * 3)


def lieb():
    a = 4.0
    positions = [[0, 0], [a / 2, 0], [0, a / 2]]
    return Crystal(primitive_vectors("square", a), positions, [63.5, 16.0, 16.0], ["Cu", "O", "O"])


def test_kagome_floppy_along_gamma_m():
    """Nearest neighbours only: the lowest branch vanishes on all three Gamma-M lines, not at K."""
    crystal = kagome()
    bonds = neighbour_bonds(crystal, [1.0])
    m = special_points("hexagonal", crystal.vectors)["M"]
    rotation = np.array([[0.5, -np.sqrt(3) / 2], [np.sqrt(3) / 2, 0.5]])
    t = np.linspace(0, 1, 11)[:, None]
    for turn in range(3):
        line = t * (np.linalg.matrix_power(rotation, turn) @ m)
        np.testing.assert_allclose(dispersion(crystal, bonds, line)[:, 0], 0, atol=1e-6)
    assert dispersion(crystal, bonds, [special_points("hexagonal", crystal.vectors)["K"]])[0, 0] > 0.1


def test_lieb_zero_modes():
    """1st shell: two branches zero everywhere; 2nd shell: only the rigid square rotation at M; 3rd: none."""
    crystal = lieb()
    points = special_points("square", crystal.vectors)
    qs, _, _ = k_path([points["Γ"], points["X"], points["M"], points["Γ"]], 60)

    nearest = dispersion(crystal, neighbour_bonds(crystal, [1.0]), qs)
    np.testing.assert_allclose(nearest[:, :2], 0, atol=1e-6)

    second = dispersion(crystal, neighbour_bonds(crystal, [1.0, 0.3]), [points["M"], points["X"]])
    assert second[0, 0] == pytest.approx(0, abs=1e-6) and second[1, 0] > 0.01

    third = dispersion(crystal, neighbour_bonds(crystal, [1.0, 0.3, 0.1]), [points["M"]])
    assert third[0, 0] > 0.01


def test_free_strip_bends_quadratically():
    """A strip periodic in x with free edges (vacuum in y): two zero modes at q = 0, lowest branch ~ q^2."""
    a, rows = 2.46, 6
    positions = [[(j * a / 2) % a, j * a * np.sqrt(3) / 2] for j in range(rows)]
    strip = Crystal([[a, 0.0], [0.0, rows * a + 10.0]], positions, [12.0] * rows, ["C"] * rows)
    bonds = neighbour_bonds(strip, [10.0])
    assert all(n[1] == 0 for _, _, n, _ in bonds)  # no bonds across the vacuum

    omega = dispersion(strip, bonds, [[0.0, 0.0], [0.01, 0.0], [0.02, 0.0]])
    np.testing.assert_allclose(omega[0, :2], 0, atol=1e-12)
    assert omega[0, 2] > 0.1
    assert omega[2, 0] / omega[1, 0] == pytest.approx(4.0, rel=0.01)
