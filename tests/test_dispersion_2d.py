import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.reciprocal import k_path, q_grid, reciprocal_vectors, special_points
from latticephys.lattice.structures import honeycomb, kagome, lieb
from latticephys.phonons.dispersion import bloch_mode, dispersion, group_velocities
from latticephys.phonons.integrate import simulate
from latticephys.phonons.modes import normal_modes
from latticephys.phonons.springs import force_constants

A = 2.46


def graphene():
    return honeycomb(A, [12.0, 12.0], ["C", "C"])


def triangular(m=2.0):
    return Crystal(primitive_vectors("hexagonal", A), [[0, 0]], [m], ["A"])


def test_triangular_analytic():
    """Bravais lattice: D(q) = sum over bonds of K/m e e^T 2 (1 - cos q.R)."""
    K, m = 3.0, 2.0
    crystal = triangular(m)
    bonds = neighbour_bonds(crystal, [K])
    qs = np.random.default_rng(0).uniform(-2, 2, size=(20, 2))

    omega = dispersion(crystal, bonds, qs)
    vectors = [crystal.bond_vector(mu, nu, n) for mu, nu, n, _ in bonds]
    for q, w in zip(qs, omega):
        D = sum(K / m * np.outer(r, r) / (r @ r) * 2 * (1 - np.cos(q @ r)) for r in vectors)
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
    omega_finite, _ = normal_modes(force_constants(supercell, supercell_bonds), np.repeat(supercell.masses, 2))

    b1, b2 = reciprocal_vectors(crystal.vectors)
    qs = [(i * b1 + j * b2) / cells for i in range(cells) for j in range(cells)]
    omega = dispersion(crystal, bonds, qs)
    np.testing.assert_allclose(np.sort(omega.ravel()), omega_finite, atol=1e-10)


def test_bloch_mode_oscillates_in_supercell():
    """Optical Bloch mode at M, set up in a supercell, evolves as u0 cos(omega t)."""
    crystal, cells = graphene(), 4
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, cells)
    Phi = force_constants(supercell, supercell_bonds)

    q = reciprocal_vectors(crystal.vectors)[0] / 2
    omega, u0 = bloch_mode(crystal, bonds, q, branch=3, cells=cells)

    dt, steps = 0.001, 3000
    u, _ = simulate(Phi, np.repeat(supercell.masses, 2), u0, np.zeros_like(u0), dt, steps)
    t = dt * np.arange(steps + 1)
    np.testing.assert_allclose(u, np.outer(np.cos(omega * t), u0), atol=1e-6)


def test_bloch_mode_rejects_q_not_fitting_supercell():
    crystal = graphene()
    q = reciprocal_vectors(crystal.vectors)[0] / 3
    with pytest.raises(ValueError):
        bloch_mode(crystal, neighbour_bonds(crystal, [1.0]), q, branch=0, cells=4)


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
    "vectors, spring_constants",
    [
        (primitive_vectors("square", 2.5), [10.0, 3.0]),
        (primitive_vectors("oblique", 2.5, 3.1, 73.0), [10.0, 3.0]),
        (primitive_vectors("hexagonal", 2.46), [10.0, 4.0, 1.0]),
    ],
)
def test_acoustic_modes_at_gamma_alone_are_exactly_zero(vectors, spring_constants):
    """Round-off at Gamma is judged against the spring scale, also when Gamma is the only q-point."""
    crystal = Crystal(vectors, [[0.0, 0.0]], [12.011], ["C"])
    bonds = neighbour_bonds(crystal, spring_constants)
    np.testing.assert_array_equal(dispersion(crystal, bonds, [[0.0, 0.0]]), 0.0)
    assert bloch_mode(crystal, bonds, [0.0, 0.0], branch=0, cells=2)[0] == 0.0


def test_kagome_floppy_along_gamma_m():
    """Nearest neighbours only: the lowest branch vanishes on all three Gamma-M lines, not at K."""
    crystal = kagome(4.0, [12.0] * 3, ["A"] * 3)
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
    crystal = lieb(4.0, [63.5, 16.0, 16.0], ["Cu", "O", "O"])
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


def test_square_with_transverse_springs_analytic():
    """Born-von Karman square lattice, nearest neighbours (K_L, K_T): D is diagonal with
    D_xx = 2/m [K_L (1 - cos q_x a) + K_T (1 - cos q_y a)] and x <-> y for D_yy."""
    a, m, K_L, K_T = 2.0, 3.0, 4.0, 1.5
    crystal = Crystal(primitive_vectors("square", a), [[0, 0]], [m], ["A"])
    qs = np.random.default_rng(1).uniform(-2, 2, size=(20, 2))

    omega = dispersion(crystal, neighbour_bonds(crystal, [(K_L, K_T)]), qs)
    cx, cy = 1 - np.cos(qs[:, 0] * a), 1 - np.cos(qs[:, 1] * a)
    expected = np.sort([2 / m * (K_L * cx + K_T * cy), 2 / m * (K_T * cx + K_L * cy)], axis=0).T
    np.testing.assert_allclose(omega**2, expected, atol=1e-12)


def test_group_velocities_match_finite_differences():
    crystal = graphene()
    bonds = neighbour_bonds(crystal, [(10.0, 1.0), 2.0])
    qs = np.random.default_rng(2).uniform(-1, 1, size=(5, 2))
    h = 1e-6

    velocities = group_velocities(crystal, bonds, qs)
    for q, v in zip(qs, velocities):
        slopes = [(dispersion(crystal, bonds, [q + h * e])[0] - dispersion(crystal, bonds, [q - h * e])[0]) / (2 * h)
                  for e in np.eye(2)]
        np.testing.assert_allclose(v, np.array(slopes).T, atol=1e-7)


def test_group_velocities_undefined_at_degeneracies():
    """At Gamma the acoustic branches are zero and the optical ones degenerate; at K some branches are degenerate."""
    crystal = graphene()
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    points = special_points("hexagonal", crystal.vectors)
    at_gamma, at_K = group_velocities(crystal, bonds, [points["Γ"], points["K"]])
    assert np.all(np.isnan(at_gamma))
    omega = dispersion(crystal, bonds, [points["K"]])[0]
    degenerate = np.isclose(omega[:, None], omega[None, :]).sum(axis=1) > 1
    assert degenerate.any()
    assert np.all(np.isnan(at_K[degenerate])) and not np.any(np.isnan(at_K[~degenerate]))


def test_graphene_layer_with_central_springs():
    """Graphene lifted into 3D (p = 2, d = 3) with in-plane central springs: the bonds are those
    of 2D graphene, nothing resists out-of-plane motion, so at every q != 0 exactly B = 2 branches
    are 0 and the other 4 are the in-plane 2D branches."""
    flat = graphene()
    layer = embed(flat, 3)
    bonds = neighbour_bonds(layer, [10.0, 2.0])
    assert bonds == neighbour_bonds(flat, [10.0, 2.0])

    qs = q_grid(layer.vectors, 6)[1:]
    omega = dispersion(layer, bonds, qs)
    omega_2d = dispersion(flat, bonds, qs[:, :2])

    assert np.all((omega == 0).sum(axis=1) == 2)
    np.testing.assert_allclose(omega, np.sort(np.column_stack([omega_2d, np.zeros((len(qs), 2))]), axis=1), atol=1e-12)


def test_graphene_layer_with_transverse_springs():
    """Nearest-neighbour (K_L, K_T) springs in the plane of a graphene layer: z decouples, and
    every bond resists z motion with K_T, so ZA, ZO: omega^2 = K_T (3 -/+ |1 + e^(i q.a1) + e^(i q.a2)|) / m
    (ZA linear near Gamma, omega = a/2 sqrt(K_T/m) q); the other 4 branches are the 2D in-plane ones."""
    K_L, K_T, m = 10.0, 2.0, 12.0
    flat = graphene()
    layer = embed(flat, 3)
    bonds = neighbour_bonds(layer, [(K_L, K_T)])
    qs = q_grid(layer.vectors, 6)

    a1, a2 = layer.vectors
    f = np.abs(1 + np.exp(1j * qs @ a1) + np.exp(1j * qs @ a2))
    out_of_plane = np.sqrt(np.clip(K_T / m * np.column_stack([3 - f, 3 + f]), 0, None))
    in_plane = dispersion(flat, bonds, qs[:, :2])

    expected = np.sort(np.column_stack([in_plane, out_of_plane]), axis=1)
    np.testing.assert_allclose(dispersion(layer, bonds, qs), expected, atol=1e-7)
