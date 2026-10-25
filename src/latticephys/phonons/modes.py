"""Dynamical matrix of a crystal, its frequencies, polarizations and atom weights; normal modes of finite systems."""

import numpy as np

from latticephys.phonons.springs import force_constants


def dynamical_matrix(crystal, bonds, q):
    """Bloch dynamical matrix D(q) = M^-1/2 Phi(q) M^-1/2, shape (B d, B d), with Phi(q) from force_constants().

    bonds: list of (mu, nu, n, K), e.g. from neighbour_bonds()
    q:     wave vector (1/Å), length d (a number in 1D), in the span of the lattice vectors
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    return _mass_weighted(force_constants(crystal, bonds, q), np.repeat(crystal.masses, crystal.dim))


def phonon_modes(crystal, bonds, qs):
    """Frequencies omega (len(qs), B d) and polarizations (len(qs), B d, B d) of D(q) at each q.

    Polarizations are the orthonormal eigenvectors of D(q), modes as columns, in the cell-phase
    convention of force_constants(); omega as in diagonalize(). For a molecule or a finite
    supercell use qs = [np.zeros(d)].
    """
    modes = [diagonalize(dynamical_matrix(crystal, bonds, q)) for q in qs]
    return np.array([omega for omega, _ in modes]), np.array([vectors for _, vectors in modes])


def diagonalize(dynamical):
    """Frequencies omega (ascending) and eigenvectors (columns) of a Hermitian dynamical matrix.

    omega = sign(omega^2) sqrt(|omega^2|), so negative values are unstable modes
    with imaginary frequency. Round-off of zero modes (|omega^2| below 1e-10 of
    the largest possible eigenvalue, the maximum absolute row sum) is set to 0.
    """
    omega_squared, vectors = np.linalg.eigh(dynamical)
    omega_squared[np.abs(omega_squared) < 1e-10 * np.abs(dynamical).sum(axis=1).max()] = 0.0
    return np.sign(omega_squared) * np.sqrt(np.abs(omega_squared)), vectors


def atom_weights(vectors, dim):
    """Weights |e_mu,i|^2 of atom mu and direction i in each mode, shape (..., B d, B, d); sum to 1 per mode.

    vectors: polarizations (..., B d, B d), modes as columns, e.g. from phonon_modes()
    dim:     dimension d of the crystal
    Sum over the last axis for one weight per atom, e.g. for density_of_states().
    """
    weights = np.swapaxes(np.abs(np.asarray(vectors)) ** 2, -1, -2)
    return weights.reshape(*weights.shape[:-1], -1, dim)


def normal_modes(force_constants, masses):
    """Normal modes of the dynamical matrix D = M^-1/2 Phi M^-1/2.

    force_constants: symmetric (N, N) matrix Phi of the potential V = 1/2 u^T Phi u
    masses:          (N,) masses

    Returns omega (ascending, see diagonalize()) and the eigenvectors of D as
    columns. The displacement pattern of mode j is masses**-0.5 * vectors[:, j].
    """
    return diagonalize(_mass_weighted(np.asarray(force_constants, dtype=float), masses))


def _mass_weighted(Phi, masses):
    """M^-1/2 Phi M^-1/2 for masses per degree of freedom; Phi of shape (..., N, N)."""
    inv_sqrt_m = 1 / np.sqrt(np.asarray(masses, dtype=float))
    return inv_sqrt_m[:, None] * Phi * inv_sqrt_m[None, :]
