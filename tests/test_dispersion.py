import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import Crystal, embed
from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.dispersion import atom_weights, bloch_mode, dispersion, dynamical_matrix, group_velocities
from latticephys.phonons.integrate import simulate
from latticephys.phonons.springs import force_constants

A = 1.0
QS = np.linspace(-np.pi / A, np.pi / A, 101)
DIATOMIC_BONDS = [(0, 1, 0, 1.0), (1, 0, 1, 1.0)]


def monatomic(m):
    return Crystal([[A]], [0.0], [m], ["A"])


def diatomic(m1, m2):
    return Crystal([[A]], [0.0, A / 2], [m1, m2], ["A", "B"])


def test_monatomic_chain():
    """omega = 2 sqrt(K/m) |sin(q a / 2)|."""
    K, m = 1.5, 2.0
    omega = dispersion(monatomic(m), [(0, 0, 1, K)], QS)

    np.testing.assert_allclose(omega[:, 0], 2 * np.sqrt(K / m) * np.abs(np.sin(QS * A / 2)), atol=1e-7)


def test_diatomic_chain():
    """omega^2 = K S -/+ K sqrt(S^2 - 4 sin^2(q a / 2) / (m1 m2)), S = 1/m1 + 1/m2."""
    K, m1, m2 = 1.0, 1.0, 3.0
    omega = dispersion(diatomic(m1, m2), DIATOMIC_BONDS, QS)

    S = 1 / m1 + 1 / m2
    root = np.sqrt(S**2 - 4 * np.sin(QS * A / 2) ** 2 / (m1 * m2))
    np.testing.assert_allclose(omega[:, 0], np.sqrt(np.clip(K * (S - root), 0, None)), atol=1e-7)
    np.testing.assert_allclose(omega[:, 1], np.sqrt(K * (S + root)), atol=1e-12)


def test_neighbour_bonds_find_diatomic_springs():
    crystal = diatomic(1.0, 3.0)
    bonds = neighbour_bonds(crystal, [1.0])
    assert sorted((mu, nu, n, K) for mu, nu, (n,), K in bonds) == [(0, 1, -1, 1.0), (0, 1, 0, 1.0)]
    np.testing.assert_allclose(dispersion(crystal, bonds, QS), dispersion(crystal, DIATOMIC_BONDS, QS))


def test_matches_finite_periodic_chain():
    """A periodic chain of N cells has exactly the frequencies at q = 2 pi j / (N a)."""
    K, m1, m2, cells = 1.0, 1.0, 3.0, 6
    qs = 2 * np.pi * np.arange(cells) / (cells * A)
    omega = dispersion(diatomic(m1, m2), DIATOMIC_BONDS, qs)

    Phi = chain_force_constants(2 * cells, K, boundary="periodic")
    masses = np.tile([m1, m2], cells)
    omega_squared = np.linalg.eigvalsh(Phi / np.sqrt(np.outer(masses, masses)))

    np.testing.assert_allclose(np.sort(omega.ravel() ** 2), omega_squared, atol=1e-12)


def test_fixed_chain_lies_on_dispersion():
    """Fixed chain of N atoms: omega_j = omega(q_j) with q_j = j pi / ((N + 1) a)."""
    K, m, n = 1.0, 2.0, 10
    omega_squared = np.linalg.eigvalsh(chain_force_constants(n, K, boundary="fixed") / m)

    qs = np.arange(1, n + 1) * np.pi / ((n + 1) * A)
    omega = dispersion(monatomic(m), [(0, 0, 1, K)], qs)[:, 0]
    np.testing.assert_allclose(np.sqrt(omega_squared), omega, atol=1e-12)


def test_simple_cubic():
    """Nearest-neighbour central springs: omega_alpha^2 = 4 K/m sin^2(q_alpha a / 2) for alpha = x, y, z."""
    K, m = 2.0, 3.0
    cubic = Crystal(A * np.eye(3), [[0.0, 0.0, 0.0]], [m], ["A"])
    qs = np.random.default_rng(0).uniform(-np.pi, np.pi, size=(20, 3)) / A
    omega = dispersion(cubic, neighbour_bonds(cubic, [K]), qs)

    expected = np.sort(4 * K / m * np.sin(qs * A / 2) ** 2, axis=1)
    np.testing.assert_allclose(omega**2, expected, atol=1e-12)


def test_wrong_q_dimension():
    with pytest.raises(ValueError):
        dynamical_matrix(monatomic(1.0), [(0, 0, 1, 1.0)], [0.1, 0.2])


def test_bloch_mode_oscillates_in_periodic_chain():
    """Optical Bloch mode at q a = pi / 2, set up in a periodic chain, evolves as u0 cos(omega t)."""
    crystal, cells = diatomic(1.0, 3.0), 8
    supercell, supercell_bonds = periodic_supercell(crystal, DIATOMIC_BONDS, cells)
    omega, u0 = bloch_mode(crystal, DIATOMIC_BONDS, np.pi / (2 * A), branch=1, cells=cells)

    dt, steps = 0.001, 5000
    u, _ = simulate(force_constants(supercell, supercell_bonds), supercell.masses, u0, np.zeros_like(u0), dt, steps)
    t = dt * np.arange(steps + 1)
    np.testing.assert_allclose(u, np.outer(np.cos(omega * t), u0), atol=1e-6)


@pytest.mark.parametrize("branch, resting", [(0, 0), (1, 1)])
def test_zone_edge_modes_move_one_sublattice(branch, resting):
    """At q = pi/a the acoustic mode moves only the heavy atoms (1), the optical mode only the light ones (0)."""
    _, u0 = bloch_mode(diatomic(1.0, 3.0), DIATOMIC_BONDS, np.pi / A, branch=branch, cells=2)
    u0 = u0.reshape(-1, 2)
    np.testing.assert_allclose(u0[:, resting], 0, atol=1e-12)
    assert np.all(np.abs(u0[:, 1 - resting]) > 0.1)


def test_bloch_mode_rejects_q_not_fitting_chain():
    with pytest.raises(ValueError):
        bloch_mode(diatomic(1.0, 3.0), DIATOMIC_BONDS, np.pi / (3 * A), branch=0, cells=8)


def test_bloch_mode_zero_mode_is_exactly_zero():
    """Round-off must not turn the acoustic mode at q = 0 into a tiny positive omega (endless period)."""
    omega, _ = bloch_mode(diatomic(1.0, 3.0), DIATOMIC_BONDS, 0.0, branch=0, cells=8)
    assert omega == 0.0


def test_group_velocity_monatomic_chain():
    """v = a sqrt(K/m) cos(q a / 2) for q > 0; undefined at q = 0."""
    K, m = 1.5, 2.0
    qs = np.linspace(0, np.pi / A, 11)
    v = group_velocities(monatomic(m), [(0, 0, 1, K)], qs)[:, 0, 0]
    assert np.isnan(v[0])
    np.testing.assert_allclose(v[1:], A * np.sqrt(K / m) * np.cos(qs[1:] * A / 2), atol=1e-12)


def test_atom_weights_at_zone_edge():
    """At q = pi/a only the heavy atom moves on the acoustic branch, only the light one on the optical."""
    weights = atom_weights(diatomic(1.0, 3.0), DIATOMIC_BONDS, [np.pi / A, 0.7])
    np.testing.assert_allclose(weights[0], [[0, 1], [1, 0]], atol=1e-12)
    np.testing.assert_allclose(weights.sum(axis=2), 1)


def test_wave_vector_outside_the_lattice_plane_raises():
    """A square layer in the xy plane of 3D space has no wave vector component along z."""
    layer = Crystal([[A, 0.0, 0.0], [0.0, A, 0.0]], [[0.0, 0.0, 0.0]], [1.0], ["A"])
    bonds = neighbour_bonds(layer, [1.0])
    dynamical_matrix(layer, bonds, [0.3, -1.0, 0.0])
    with pytest.raises(ValueError, match="plane|span"):
        dynamical_matrix(layer, bonds, [0.3, -1.0, 1e-3])


def test_embedded_chain_keeps_its_frequencies():
    """Diatomic chain lifted into the plane, central springs along x: for q along x the
    longitudinal branches are the 1D ones, the B transverse (y) branches are exactly 0."""
    planar = embed(diatomic(1.0, 3.0), 2)
    omega_1d = dispersion(diatomic(1.0, 3.0), DIATOMIC_BONDS, QS)
    omega_2d = dispersion(planar, DIATOMIC_BONDS, np.column_stack([QS, np.zeros_like(QS)]))

    expected = np.sort(np.column_stack([omega_1d, np.zeros_like(omega_1d)]), axis=1)
    np.testing.assert_allclose(omega_2d, expected, atol=1e-12)


def test_fcc_nearest_neighbours_match_the_analytic_dynamical_matrix():
    """fcc, central springs K to the 12 neighbours delta = a/2 (+-1, +-1, 0) and permutations:
    D(q) = K/m sum_delta delta delta^T / |delta|^2 (1 - cos q.delta); compare the eigenvalues at random q."""
    a, K, m = 3.6, 2.0, 63.5
    fcc = Crystal(primitive_vectors("fcc", a), [[0.0, 0.0, 0.0]], [m], ["Cu"])
    bonds = neighbour_bonds(fcc, [K])
    assert len(bonds) == 6  # 12 neighbours, each bond listed once

    deltas = [a / 2 * np.roll([s1, s2, 0], shift) for s1 in (1, -1) for s2 in (1, -1) for shift in range(3)]
    qs = np.random.default_rng(11).uniform(-2, 2, size=(20, 3))
    for q, omega in zip(qs, dispersion(fcc, bonds, qs)):
        D = sum(K / m * np.outer(d, d) / (d @ d) * (1 - np.cos(q @ d)) for d in deltas)
        np.testing.assert_allclose(omega, np.sqrt(np.linalg.eigvalsh(D)), atol=1e-10)


def test_bcc_first_and_second_neighbours_match_the_analytic_dynamical_matrix():
    """bcc, central K1 to the 8 neighbours a/2 (+-1, +-1, +-1) and K2 to the 6 at a (+-1, 0, 0) and permutations:
    D(q) = sum_delta K/m delta delta^T / |delta|^2 (1 - cos q.delta); compare the eigenvalues at random q."""
    a, K1, K2, m = 2.87, 2.0, 0.8, 55.8
    bcc = Crystal(primitive_vectors("bcc", a), [[0.0, 0.0, 0.0]], [m], ["Fe"])
    bonds = neighbour_bonds(bcc, [K1, K2])
    assert len(bonds) == 4 + 3  # 8 and 6 neighbours, each bond listed once

    first = [a / 2 * np.array([s1, s2, s3]) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
    second = [a * s * np.eye(3)[i] for s in (1, -1) for i in range(3)]
    qs = np.random.default_rng(13).uniform(-2, 2, size=(20, 3))
    for q, omega in zip(qs, dispersion(bcc, bonds, qs)):
        D = sum(K / m * np.outer(d, d) / (d @ d) * (1 - np.cos(q @ d)) for K, ds in ((K1, first), (K2, second))
                for d in ds)
        np.testing.assert_allclose(omega, np.sqrt(np.clip(np.linalg.eigvalsh(D), 0, None)), atol=1e-10)
