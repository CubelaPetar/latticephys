import numpy as np
import pytest

from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.integrate import simulate
from latticephys.phonons.modes import normal_modes


def energy(Phi, masses, u, p):
    return 0.5 * np.sum(p**2 / masses, axis=-1) + 0.5 * np.einsum("ti,ij,tj->t", u, Phi, u)


def test_energy_conserved():
    Phi = chain_force_constants(5, 1.0, boundary="fixed")
    masses = np.array([1.0, 2.0, 1.0, 2.0, 1.0])
    u, p = simulate(Phi, masses, [0.3, 0, 0, 0, 0], np.zeros(5), dt=0.01, steps=100_000)

    E = energy(Phi, masses, u, p)
    assert np.max(np.abs(E - E[0])) / E[0] < 1e-4


def test_normal_mode_oscillates_harmonically():
    """Starting in a normal mode, u(t) = u0 cos(omega t)."""
    Phi = chain_force_constants(3, 1.0, boundary="free")
    masses = np.array([1.0, 2.0, 1.0])
    omegas, vectors = normal_modes(Phi, masses)
    omega = omegas[2]
    u0 = 0.1 * vectors[:, 2] / np.sqrt(masses)

    dt, steps = 0.001, 10_000
    u, _ = simulate(Phi, masses, u0, np.zeros(3), dt, steps)

    t = dt * np.arange(steps + 1)
    np.testing.assert_allclose(u, np.outer(np.cos(omega * t), u0), atol=1e-5)


@pytest.mark.parametrize("q0a", [0.3 * np.pi, 0.8 * np.pi])
def test_wave_packet_moves_at_the_group_velocity_and_spreads(q0a):
    """Gaussian packet u = envelope cos(q0 x - omega0 t), envelope width w, in a periodic chain: its energy centre
    moves at v_g = a sqrt(K/m) cos(q0 a / 2) (slow near the zone edge), and dispersion widens the envelope to
    sigma(t)^2 = w^2 + (omega'' t / w)^2 with omega'' = -a^2 / 2 sqrt(K/m) sin(q0 a / 2); the energy width is
    sigma / sqrt(2). Measured within 40 sites of the centre: the start is not exactly one-way and leaves a
    small tail behind."""
    n, K, m, width, dt, steps = 300, 10.0, 12.0, 10, 0.05, 2000
    sites, start = np.arange(n), n // 4
    omega0 = 2 * np.sqrt(K / m) * np.sin(q0a / 2)
    envelope = np.exp(-((sites - start) ** 2) / (2 * width**2))
    u, p = simulate(chain_force_constants(n, K, boundary="periodic"), np.full(n, m),
                    envelope * np.cos(q0a * sites), m * omega0 * envelope * np.sin(q0a * sites), dt, steps)

    weight = u**2 + (p / (m * omega0)) ** 2  # local energy, smooth over the carrier oscillation
    centre = (weight * sites).sum(axis=1) / weight.sum(axis=1)
    t = dt * steps
    assert centre[-1] - centre[0] == pytest.approx(np.sqrt(K / m) * np.cos(q0a / 2) * t, rel=0.01)

    near = np.abs(sites - centre[-1]) < 40
    spread = np.sqrt((weight[-1, near] * (sites[near] - centre[-1]) ** 2).sum() / weight[-1, near].sum())
    curvature = -0.5 * np.sqrt(K / m) * np.sin(q0a / 2)
    assert spread == pytest.approx(np.sqrt(width**2 + (curvature * t / width) ** 2) / np.sqrt(2), rel=0.01)
