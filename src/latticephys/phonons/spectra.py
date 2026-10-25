"""Phonon spectra: density of states and the one-phonon structure factor S(Q, omega) of scattering experiments."""

import numpy as np

from latticephys.lattice.reciprocal import lattice_projection
from latticephys.phonons.modes import phonon_modes


def density_of_states(omega, grid, sigma, weights=None):
    """Gaussian-broadened density of states g(grid) = 1/N sum_i w_i G_sigma(grid - omega_i).

    omega:   N mode frequencies of any shape, e.g. dispersion() on a q_grid()
    grid:    frequencies where g is evaluated, same units as omega
    sigma:   Gaussian width; should exceed the level spacing of the sampled modes
    weights: optional w_i, shape of omega, or omega.shape + (B,) for one DOS per
             atom, e.g. atom_weights().sum(axis=-1); without weights (w_i = 1) g integrates to 1
    Returns shape (len(grid),), or (len(grid), B) with per-atom weights.
    """
    omega = np.asarray(omega, dtype=float)
    flat = omega.ravel()
    w = np.ones(flat.size) if weights is None else np.asarray(weights, dtype=float)
    w = w.reshape(flat.size, *w.shape[omega.ndim :])
    grid = np.asarray(grid, dtype=float)

    result = np.empty((len(grid),) + w.shape[1:])
    for start in range(0, len(grid), 64):  # chunks keep the (grid, modes) matrix small
        offset = grid[start : start + 64, None] - flat[None, :]
        result[start : start + 64] = np.exp(-(offset**2) / (2 * sigma**2)) @ w
    return result / (flat.size * np.sqrt(2 * np.pi) * sigma)


def structure_factor(crystal, bonds, Qs, scattering_lengths=None):
    """Frequencies and one-phonon intensities at momentum transfers Q (1/Å), each shape (len(Qs), B d).

    Intensity of branch s: |sum_mu b_mu / sqrt(m_mu) (Q . e_mu,s) exp(-i Q . r_mu)|^2 / omega_s,
    with scattering lengths b_mu (default 1) and polarizations e in the cell-phase convention of
    force_constants(), hence the atom phase exp(-i Q . r_mu) here. Q is not reduced to the first zone, so the
    intensity tells in which Brillouin zone a branch is seen. Zero-temperature phonon
    creation in arbitrary units; multiply by n(omega) + 1 for temperature T. The Debye-Waller
    factor is left out. Zero modes (elastic, Bragg) get intensity 0; unstable modes raise.
    For p < d (e.g. a layer) Q may leave the lattice plane: the phonon has the in-plane part
    of Q, Q . e and the phases use the full Q.
    """
    d, n_basis = crystal.dim, len(crystal.masses)
    b = np.ones(n_basis) if scattering_lengths is None else np.asarray(scattering_lengths, dtype=float)

    Qs = [np.atleast_1d(np.asarray(Q, dtype=float)) for Q in Qs]
    frequencies, all_vectors = phonon_modes(crystal, bonds, [lattice_projection(crystal.vectors, Q) for Q in Qs])
    intensities = []
    for Q, omega, vectors in zip(Qs, frequencies, all_vectors):
        if np.any(omega < 0):
            raise ValueError(f"unstable mode at Q = {Q}, no structure factor")

        polarizations = vectors.reshape(n_basis, d, -1)
        prefactor = b / np.sqrt(crystal.masses) * np.exp(-1j * crystal.positions @ Q)
        amplitude = np.einsum("m,i,mis->s", prefactor, Q, polarizations)
        with np.errstate(divide="ignore", invalid="ignore"):
            intensity = np.where(omega > 0, np.abs(amplitude) ** 2 / omega, 0.0)
        intensities.append(intensity)
    return frequencies, np.array(intensities)


def intensity_map(frequencies, intensities, grid, sigma):
    """S(Q, omega) on a frequency grid: per Q the intensities broadened by Gaussians of width sigma.

    frequencies, intensities: the output of structure_factor(), shape (len(Qs), B d)
    Returns shape (len(grid), len(Qs)), for imshow with origin="lower"; each column is
    density_of_states() of that Q with the intensities as weights.
    """
    return np.array([density_of_states(w, grid, sigma, weights=i) for w, i in zip(frequencies, intensities)]).T
