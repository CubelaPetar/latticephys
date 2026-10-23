"""Phonon dispersion of a Crystal with spring bonds, in any dimension."""

import numpy as np

from latticephys.lattice.crystal import cell_offsets
from latticephys.lattice.reciprocal import reciprocal_vectors
from latticephys.phonons.modes import diagonalize
from latticephys.phonons.springs import _fourier_force_constants


def dynamical_matrix(crystal, bonds, q):
    """Bloch dynamical matrix D(q) = M^-1/2 Phi(q) M^-1/2, shape (B d, B d).

    bonds: list of (mu, nu, n, K), e.g. from neighbour_bonds()
    q:     wave vector (1/Å), length d (a number in 1D), in the span of the lattice vectors
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.shape != (crystal.dim,):
        raise ValueError(f"q must have {crystal.dim} components, got {q.shape}")
    in_plane = crystal.vectors @ q / (2 * np.pi) @ reciprocal_vectors(crystal.vectors)
    if np.linalg.norm(q - in_plane) > 1e-8 * np.linalg.norm(q):
        raise ValueError(f"q = {q} leaves the span of the lattice vectors")

    inv_sqrt_m = np.repeat(1 / np.sqrt(crystal.masses), crystal.dim)
    return inv_sqrt_m[:, None] * _fourier_force_constants(crystal, bonds, q) * inv_sqrt_m[None, :]


def dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q), shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return np.array([diagonalize(dynamical_matrix(crystal, bonds, q))[0] for q in qs])


def sound_velocities(crystal, bonds, directions, dq=1e-4):
    """Velocities omega / q of the d acoustic branches for small q along each direction.

    directions: (M, d) directions (normalized here)
    dq:         |q| used for the slope (1/Å); relative error ~ (dq a)^2
    Returns (M, d) velocities in Å per internal time unit, ascending per direction
    (transverse before longitudinal); see units.VELOCITY_TO_KM_S.
    """
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    return dispersion(crystal, bonds, dq * directions)[:, : crystal.dim] / dq


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
    omega, vectors = diagonalize(dynamical_matrix(crystal, bonds, q))
    e = vectors[:, branch]
    # eigh returns e with an arbitrary complex phase; with e.g. e = (0, i) the real part
    # below would vanish everywhere, so make the largest component real
    e = e * np.exp(-1j * np.angle(e[np.argmax(np.abs(e))]))
    polarization = (e / np.repeat(np.sqrt(crystal.masses), d)).reshape(-1, d)

    cell_vectors = cell_offsets(cells, d) @ crystal.vectors
    u0 = np.real(np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None])
    return omega[branch], u0.ravel()


def group_velocities(crystal, bonds, qs):
    """Group velocities d omega / dq of all branches, shape (len(qs), B d, d), in Å per internal time unit.

    From the eigenvectors (Hellmann-Feynman): d omega / dq_i = e^H (dD / dq_i) e / (2 omega).
    The velocity is not defined for zero modes and for degenerate branches (it depends on
    the chosen combination there), so these entries are nan.
    """
    d = crystal.dim
    inv_sqrt_m = np.repeat(1 / np.sqrt(crystal.masses), d)
    velocities = []
    for q in qs:
        q = np.atleast_1d(np.asarray(q, dtype=float))
        omega, vectors = diagonalize(dynamical_matrix(crystal, bonds, q))
        gradient = inv_sqrt_m[:, None] * _fourier_force_constants(crystal, bonds, q, gradient=True) * inv_sqrt_m
        slope = np.einsum("as,iab,bs->si", vectors.conj(), gradient, vectors).real

        gap = np.abs(omega[:, None] - omega[None, :])
        np.fill_diagonal(gap, np.inf)
        undefined = (omega == 0) | (gap.min(axis=1) < 1e-6 * np.abs(omega).max())
        with np.errstate(divide="ignore", invalid="ignore"):
            velocity = slope / (2 * omega[:, None])
        velocity[undefined] = np.nan
        velocities.append(velocity)
    return np.array(velocities)


def atom_weights(crystal, bonds, qs):
    """Weight |e_mu|^2 of each basis atom in each mode, shape (len(qs), B d, B); sums to 1 per mode."""
    d, n_basis = crystal.dim, len(crystal.masses)
    weights = []
    for q in qs:
        _, vectors = diagonalize(dynamical_matrix(crystal, bonds, q))
        weights.append((np.abs(vectors) ** 2).reshape(n_basis, d, -1).sum(axis=1).T)
    return np.array(weights)
