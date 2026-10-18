"""Phonon dispersion of a Crystal with spring bonds, in any dimension."""

import numpy as np

from libsolidstate.lattice.crystal import cell_offsets
from libsolidstate.lattice.springs import _fourier_force_constants


def dynamical_matrix(crystal, bonds, q):
    """Bloch dynamical matrix D(q) = M^-1/2 Phi(q) M^-1/2, shape (B d, B d).

    bonds: list of (mu, nu, n, k), e.g. from neighbour_bonds()
    q:     wave vector (1/Å), length d (a number in 1D)
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.shape != (crystal.dim,):
        raise ValueError(f"q must have {crystal.dim} components, got {q.shape}")

    inv_sqrt_m = np.repeat(1 / np.sqrt(crystal.masses), crystal.dim)
    return inv_sqrt_m[:, None] * _fourier_force_constants(crystal, bonds, q) * inv_sqrt_m[None, :]


def dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q), shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return _frequencies([dynamical_matrix(crystal, bonds, q) for q in qs])


def bloch_mode(crystal, bonds, q, branch, cells):
    """Frequency and displacements of one Bloch mode in a supercell of cells^d cells.

    q must fit the supercell, q = (i b1 + j b2 + ...) / cells. For degenerate
    branches the returned mode is one arbitrary combination. Displacements
    (B d cells^d,) start the standing wave Re[e_mu exp(i q.R)] / sqrt(m_mu) for
    simulate(), in the atom order of periodic_supercell().
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    fit = cells * crystal.vectors @ q / (2 * np.pi)
    if not np.allclose(fit, np.round(fit), atol=1e-8):
        raise ValueError("q does not fit the supercell, use q = (i b1 + j b2 + ...) / cells")

    d = crystal.dim
    omega_squared, vectors = np.linalg.eigh(dynamical_matrix(crystal, bonds, q))
    polarization = (vectors[:, branch] / np.repeat(np.sqrt(crystal.masses), d)).reshape(-1, d)

    cell_vectors = cell_offsets(cells, d) @ crystal.vectors
    u0 = np.real(np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None])
    omega = np.sign(omega_squared[branch]) * np.sqrt(np.abs(omega_squared[branch]))
    return omega, u0.ravel()


def _frequencies(matrices):
    """omega = sign(omega^2) sqrt(|omega^2|), with round-off around zero set to 0."""
    omega_squared = np.array([np.linalg.eigvalsh(D) for D in matrices])
    omega_squared[np.abs(omega_squared) < 1e-10 * np.abs(omega_squared).max()] = 0.0
    return np.sign(omega_squared) * np.sqrt(np.abs(omega_squared))
