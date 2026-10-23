import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.reciprocal import reciprocal_vectors
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import atom_weights, dispersion
from latticephys.phonons.modes import normal_modes
from latticephys.phonons.spectra import density_of_states, structure_factor
from latticephys.phonons.springs import force_constants

A = 1.0
DIATOMIC_BONDS = [(0, 1, 0, 1.0), (1, 0, 1, 1.0)]


def chain_modes(crystal, bonds, n=40_000):
    """Frequencies on n midpoints of the Brillouin zone."""
    qs = ((np.arange(n) + 0.5) / n - 0.5) * 2 * np.pi / A
    return dispersion(crystal, bonds, qs)


def test_dos_monatomic_chain():
    """g(omega) = 2 / (pi sqrt(omega_max^2 - omega^2)), normalized to 1."""
    K, m = 1.5, 2.0
    omega = chain_modes(Crystal([[A]], [0.0], [m], ["A"]), [(0, 0, 1, K)])
    omega_max = 2 * np.sqrt(K / m)
    grid = np.linspace(0.1, 0.8, 8) * omega_max

    dos = density_of_states(omega, grid, sigma=1e-3 * omega_max)
    np.testing.assert_allclose(dos, 2 / (np.pi * np.sqrt(omega_max**2 - grid**2)), rtol=1e-3)

    wide = np.linspace(-0.5, 1.5, 4001) * omega_max
    assert np.trapezoid(density_of_states(omega, wide, sigma=0.01 * omega_max), wide) == pytest.approx(1, rel=1e-6)


def test_projected_dos_adds_up_and_has_the_gap():
    """Per-atom DOS sums to the total; no modes between the acoustic top sqrt(2K/m2) and optical bottom sqrt(2K/m1)."""
    crystal = Crystal([[A]], [0.0, A / 2], [1.0, 3.0], ["A", "B"])
    qs = ((np.arange(4000) + 0.5) / 4000 - 0.5) * 2 * np.pi / A
    omega = dispersion(crystal, DIATOMIC_BONDS, qs)
    weights = atom_weights(crystal, DIATOMIC_BONDS, qs)
    grid = np.linspace(0, 2.2, 200)

    total = density_of_states(omega, grid, sigma=0.01)
    projected = density_of_states(omega, grid, sigma=0.01, weights=weights)
    np.testing.assert_allclose(projected.sum(axis=1), total, atol=1e-12)

    gap = (grid > np.sqrt(2 / 3) + 0.1) & (grid < np.sqrt(2) - 0.1)
    assert np.all(total[gap] < 1e-6)


def hbn():
    return honeycomb(2.5, [10.81, 14.007], ["B", "N"])


def test_structure_factor_f_sum_rule():
    """sum_s omega_s I_s(Q) = sum_mu b_mu^2 |Q|^2 / m_mu by completeness of the eigenvectors."""
    crystal, b = hbn(), np.array([0.665, 0.936])
    Qs = np.random.default_rng(3).uniform(-4, 4, size=(10, 2))
    omega, intensity = structure_factor(crystal, neighbour_bonds(crystal, [10.0, 2.0]), Qs, b)

    expected = (b**2 / crystal.masses).sum() * (Qs**2).sum(axis=1)
    np.testing.assert_allclose((omega * intensity).sum(axis=1), expected, rtol=1e-10)


def test_structure_factor_matches_supercell():
    """Independent of the Bloch phase convention: normal modes of a periodic supercell, with the
    amplitude sum_j b_j / sqrt(m_j) (Q . e_j) exp(-i Q . x_j) over all atoms at their absolute
    positions x_j, carry N |amplitude|^2 at each Bloch frequency, for Q outside the first zone."""
    crystal, b, cells = hbn(), np.array([0.665, 0.936]), 3
    bonds = neighbour_bonds(crystal, [10.0, 2.0])
    b1, b2 = reciprocal_vectors(crystal.vectors)
    Q = b1 / cells + b1 + b2

    omega, intensity = structure_factor(crystal, bonds, [Q], b)
    bloch = omega[0] * intensity[0]  # |amplitude|^2 per branch

    supercell, supercell_bonds = periodic_supercell(crystal, bonds, cells)
    omega_sc, vectors = normal_modes(force_constants(supercell, supercell_bonds), np.repeat(supercell.masses, 2))
    prefactor = np.tile(b, cells**2) / np.sqrt(supercell.masses) * np.exp(-1j * supercell.positions @ Q)
    amplitude_sc = np.abs(np.einsum("j,i,jis->s", prefactor, Q, vectors.reshape(-1, 2, len(omega_sc)))) ** 2

    for w in omega[0]:
        same_sc = np.isclose(omega_sc, w, atol=1e-8)
        same_bloch = np.isclose(omega[0], w, atol=1e-8)
        assert amplitude_sc[same_sc].sum() == pytest.approx(cells**2 * bloch[same_bloch].sum(), rel=1e-8)
    assert amplitude_sc.sum() == pytest.approx(cells**2 * bloch.sum(), rel=1e-8)


def test_acoustic_extinction_in_odd_zones():
    """Diatomic chain with equal scattering lengths, atoms at 0 and a/2: near Q = 2 pi / a the two atoms
    scatter in antiphase and the acoustic branch is invisible; near Q = 4 pi / a it is strong."""
    crystal = Crystal([[A]], [0.0, A / 2], [1.0, 3.0], ["A", "B"])
    delta = 0.05
    _, intensity = structure_factor(crystal, DIATOMIC_BONDS, [2 * np.pi / A + delta, 4 * np.pi / A + delta])
    assert intensity[0, 0] < 1e-3 * intensity[1, 0]


def test_structure_factor_rejects_unstable_modes():
    crystal = Crystal([[A]], [0.0], [1.0], ["A"])
    with pytest.raises(ValueError):
        structure_factor(crystal, [(0, 0, 1, -1.0)], [1.0])
