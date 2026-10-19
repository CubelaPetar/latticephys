import numpy as np
import pytest

from libsolidstate.lattice.bonds import neighbour_bonds
from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.reciprocal import q_grid
from libsolidstate.phonons.dispersion import dispersion, sound_velocities
from libsolidstate.phonons.thermal import heat_capacity, mean_square_displacement
from libsolidstate.units import K_B, OMEGA_TO_MEV

A = 2.46
K, M = 10.0, 12.011


def triangular():
    return Crystal(primitive_vectors("hexagonal", A), [[0.0, 0.0]], [M], ["C"])


def test_q_grid_matches_supercell_wave_vectors():
    vectors = primitive_vectors("hexagonal", A)
    qs = q_grid(vectors, 4)
    assert qs.shape == (16, 2)
    np.testing.assert_allclose(4 * qs @ vectors.T / (2 * np.pi), np.round(4 * qs @ vectors.T / (2 * np.pi)), atol=1e-12)


def test_heat_capacity_limits():
    """Einstein mode: C = x^2 e^x / (e^x - 1)^2 with x = hbar omega / k_B T; 1 for x -> 0, 0 for x -> inf."""
    omega = np.array([1.0])
    x = OMEGA_TO_MEV / 1000 / (K_B * 300.0)
    assert heat_capacity(omega, [300.0])[0] == pytest.approx(x**2 * np.exp(x) / np.expm1(x) ** 2)
    assert heat_capacity([0.0], [300.0])[0] == 1.0
    assert heat_capacity(omega, [1e6])[0] == pytest.approx(1.0, rel=1e-6)
    assert heat_capacity(omega, [0.1])[0] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError):
        heat_capacity([-1.0], [300.0])


def test_2d_debye_law():
    """Low T: C ~ T^2 in 2D; high T: Dulong-Petit, k_B per mode.

    T must stay above the lowest mode of the finite grid (~ 1 meV ~ 10 K for 100 x 100 cells).
    """
    crystal = triangular()
    omega = dispersion(crystal, neighbour_bonds(crystal, [K]), q_grid(crystal.vectors, 100))
    low = heat_capacity(omega, [10.0, 20.0])
    assert low[1] / low[0] == pytest.approx(4.0, rel=0.02)
    assert heat_capacity(omega, [1e5])[0] == pytest.approx(1.0, rel=1e-4)


def test_triangular_cauchy_relation():
    """Nearest-neighbour central springs on the triangular lattice: isotropic, c_L / c_T = sqrt(3)."""
    crystal = triangular()
    angles = np.linspace(0, np.pi, 7)
    v = sound_velocities(crystal, neighbour_bonds(crystal, [K]), np.column_stack([np.cos(angles), np.sin(angles)]))

    np.testing.assert_allclose(v[:, 1] / v[:, 0], np.sqrt(3), rtol=1e-6)
    np.testing.assert_allclose(v, np.tile(v[0], (len(v), 1)), rtol=1e-6)  # same in every direction
    # long-wavelength limit of D(q) for the three bond directions: c_T = a sqrt(3 k / (8 m))
    assert v[0, 0] == pytest.approx(A * np.sqrt(3 * K / (8 * M)), rel=1e-6)


def test_mean_square_displacement_grows_logarithmically():
    """2D: <u^2> = k_B T A_cell / (2 pi m) (1/c_L^2 + 1/c_T^2) ln(N) + const (Mermin-Wagner)."""
    crystal, temperature = triangular(), 300.0
    bonds = neighbour_bonds(crystal, [K])
    c_t, c_l = sound_velocities(crystal, bonds, [[1.0, 0.0]])[0]
    area = abs(np.linalg.det(crystal.vectors))
    expected = K_B * temperature * area / (2 * np.pi * M) * (1 / c_l**2 + 1 / c_t**2)

    u2 = [mean_square_displacement(crystal, bonds, n, temperature)[0] for n in (32, 64)]
    assert (u2[1] - u2[0]) / np.log(2) == pytest.approx(expected, rel=0.02)


def test_mean_square_displacement_rejects_floppy_lattice():
    square = Crystal(np.eye(2), [[0.0, 0.0]], [1.0], ["A"])
    with pytest.raises(ValueError):
        mean_square_displacement(square, neighbour_bonds(square, [1.0]), 8, 300.0)


def test_mean_square_displacement_1d_exact():
    """Chain of N atoms: <u^2> = k_B T (N^2 - 1) / (12 N k), linear in N."""
    chain, temperature, n = Crystal([[1.5]], [0.0], [M], ["C"]), 300.0, 64
    u2 = mean_square_displacement(chain, [(0, 0, 1, K)], n, temperature)[0]
    assert u2 == pytest.approx(K_B * temperature * (n**2 - 1) / (12 * n * K), rel=1e-12)


def test_mean_square_displacement_3d_converges():
    """Simple cubic with 1st + 2nd neighbours: <u^2> approaches a finite limit, differences shrink."""
    cubic = Crystal(2.5 * np.eye(3), [[0.0, 0.0, 0.0]], [M], ["C"])
    bonds = neighbour_bonds(cubic, [K, 3.0])
    u2 = [mean_square_displacement(cubic, bonds, n, 300.0)[0] for n in (8, 16, 24)]
    assert 0 < u2[2] - u2[1] < 0.5 * (u2[1] - u2[0])
