import numpy as np

from libsolidstate.lattice.chain import chain_force_constants
from libsolidstate.spectra.dispersion import dispersion

A = 1.0
QS = np.linspace(-np.pi / A, np.pi / A, 101)


def test_monatomic_chain():
    """omega = 2 sqrt(k/m) |sin(q a / 2)|."""
    k, m = 1.5, 2.0
    omega = dispersion(QS, [m], [(0, 0, 1, k)], A)

    np.testing.assert_allclose(omega[:, 0], 2 * np.sqrt(k / m) * np.abs(np.sin(QS * A / 2)), atol=1e-7)


def test_diatomic_chain():
    """omega^2 = k S -/+ k sqrt(S^2 - 4 sin^2(q a / 2) / (m1 m2)), S = 1/m1 + 1/m2."""
    k, m1, m2 = 1.0, 1.0, 3.0
    omega = dispersion(QS, [m1, m2], [(0, 1, 0, k), (1, 0, 1, k)], A)

    S = 1 / m1 + 1 / m2
    root = np.sqrt(S**2 - 4 * np.sin(QS * A / 2) ** 2 / (m1 * m2))
    np.testing.assert_allclose(omega[:, 0], np.sqrt(np.clip(k * (S - root), 0, None)), atol=1e-7)
    np.testing.assert_allclose(omega[:, 1], np.sqrt(k * (S + root)), atol=1e-12)


def test_matches_finite_periodic_chain():
    """A periodic chain of N cells has exactly the frequencies at q = 2 pi j / (N a)."""
    k, m1, m2, cells = 1.0, 1.0, 3.0, 6
    qs = 2 * np.pi * np.arange(cells) / (cells * A)
    omega = dispersion(qs, [m1, m2], [(0, 1, 0, k), (1, 0, 1, k)], A)

    K = chain_force_constants(2 * cells, k, boundary="periodic")
    masses = np.tile([m1, m2], cells)
    omega_squared = np.linalg.eigvalsh(K / np.sqrt(np.outer(masses, masses)))

    np.testing.assert_allclose(np.sort(omega.ravel() ** 2), omega_squared, atol=1e-12)
