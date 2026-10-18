import numpy as np
import pytest

from libsolidstate.lattice.chain import chain_force_constants
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.springs import neighbour_bonds
from libsolidstate.spectra.dispersion import dispersion, dynamical_matrix

A = 1.0
QS = np.linspace(-np.pi / A, np.pi / A, 101)
DIATOMIC_BONDS = [(0, 1, 0, 1.0), (1, 0, 1, 1.0)]


def monatomic(m):
    return Crystal([[A]], [0.0], [m], ["A"])


def diatomic(m1, m2):
    return Crystal([[A]], [0.0, A / 2], [m1, m2], ["A", "B"])


def test_monatomic_chain():
    """omega = 2 sqrt(k/m) |sin(q a / 2)|."""
    k, m = 1.5, 2.0
    omega = dispersion(monatomic(m), [(0, 0, 1, k)], QS)

    np.testing.assert_allclose(omega[:, 0], 2 * np.sqrt(k / m) * np.abs(np.sin(QS * A / 2)), atol=1e-7)


def test_diatomic_chain():
    """omega^2 = k S -/+ k sqrt(S^2 - 4 sin^2(q a / 2) / (m1 m2)), S = 1/m1 + 1/m2."""
    k, m1, m2 = 1.0, 1.0, 3.0
    omega = dispersion(diatomic(m1, m2), DIATOMIC_BONDS, QS)

    S = 1 / m1 + 1 / m2
    root = np.sqrt(S**2 - 4 * np.sin(QS * A / 2) ** 2 / (m1 * m2))
    np.testing.assert_allclose(omega[:, 0], np.sqrt(np.clip(k * (S - root), 0, None)), atol=1e-7)
    np.testing.assert_allclose(omega[:, 1], np.sqrt(k * (S + root)), atol=1e-12)


def test_neighbour_bonds_find_diatomic_springs():
    crystal = diatomic(1.0, 3.0)
    bonds = neighbour_bonds(crystal, [1.0])
    assert sorted((mu, nu, n, k) for mu, nu, (n,), k in bonds) == [(0, 1, -1, 1.0), (0, 1, 0, 1.0)]
    np.testing.assert_allclose(dispersion(crystal, bonds, QS), dispersion(crystal, DIATOMIC_BONDS, QS))


def test_matches_finite_periodic_chain():
    """A periodic chain of N cells has exactly the frequencies at q = 2 pi j / (N a)."""
    k, m1, m2, cells = 1.0, 1.0, 3.0, 6
    qs = 2 * np.pi * np.arange(cells) / (cells * A)
    omega = dispersion(diatomic(m1, m2), DIATOMIC_BONDS, qs)

    K = chain_force_constants(2 * cells, k, boundary="periodic")
    masses = np.tile([m1, m2], cells)
    omega_squared = np.linalg.eigvalsh(K / np.sqrt(np.outer(masses, masses)))

    np.testing.assert_allclose(np.sort(omega.ravel() ** 2), omega_squared, atol=1e-12)


def test_fixed_chain_lies_on_dispersion():
    """Fixed chain of N atoms: omega_j = omega(q_j) with q_j = j pi / ((N + 1) a)."""
    k, m, n = 1.0, 2.0, 10
    omega_squared = np.linalg.eigvalsh(chain_force_constants(n, k, boundary="fixed") / m)

    qs = np.arange(1, n + 1) * np.pi / ((n + 1) * A)
    omega = dispersion(monatomic(m), [(0, 0, 1, k)], qs)[:, 0]
    np.testing.assert_allclose(np.sqrt(omega_squared), omega, atol=1e-12)


def test_simple_cubic():
    """Nearest-neighbour central springs: omega_alpha^2 = 4 k/m sin^2(q_alpha a / 2) for alpha = x, y, z."""
    k, m = 2.0, 3.0
    cubic = Crystal(A * np.eye(3), [[0.0, 0.0, 0.0]], [m], ["A"])
    qs = np.random.default_rng(0).uniform(-np.pi, np.pi, size=(20, 3)) / A
    omega = dispersion(cubic, neighbour_bonds(cubic, [k]), qs)

    expected = np.sort(4 * k / m * np.sin(qs * A / 2) ** 2, axis=1)
    np.testing.assert_allclose(omega**2, expected, atol=1e-12)


def test_wrong_q_dimension():
    with pytest.raises(ValueError):
        dynamical_matrix(monatomic(1.0), [(0, 0, 1, 1.0)], [0.1, 0.2])
