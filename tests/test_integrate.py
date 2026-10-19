import numpy as np

from libsolidstate.phonons.chain import chain_force_constants
from libsolidstate.phonons.integrate import simulate
from libsolidstate.phonons.modes import normal_modes


def energy(K, masses, u, p):
    return 0.5 * np.sum(p**2 / masses, axis=-1) + 0.5 * np.einsum("ti,ij,tj->t", u, K, u)


def test_energy_conserved():
    K = chain_force_constants(5, 1.0, boundary="fixed")
    masses = np.array([1.0, 2.0, 1.0, 2.0, 1.0])
    u, p = simulate(K, masses, [0.3, 0, 0, 0, 0], np.zeros(5), dt=0.01, steps=100_000)

    E = energy(K, masses, u, p)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-4


def test_normal_mode_oscillates_harmonically():
    """Starting in a normal mode, u(t) = u0 cos(omega t)."""
    K = chain_force_constants(3, 1.0, boundary="free")
    masses = np.array([1.0, 2.0, 1.0])
    omegas, vectors = normal_modes(K, masses)
    omega = omegas[2]
    u0 = 0.1 * vectors[:, 2] / np.sqrt(masses)

    dt, steps = 0.001, 10_000
    u, _ = simulate(K, masses, u0, np.zeros(3), dt, steps)

    t = dt * np.arange(steps + 1)
    np.testing.assert_allclose(u, np.outer(np.cos(omega * t), u0), atol=1e-5)
