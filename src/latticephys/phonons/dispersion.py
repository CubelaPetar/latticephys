"""Phonon dispersion of a Crystal with spring bonds, in any dimension."""

import numpy as np

from latticephys.lattice.crystal import cell_offsets
from latticephys.phonons.modes import _mass_weighted, diagonalize
from latticephys.phonons.springs import _force_constants_gradient, force_constants


def dynamical_matrix(crystal, bonds, q):
    """Bloch dynamical matrix D(q) = M^-1/2 Phi(q) M^-1/2, shape (B d, B d), with Phi(q) from force_constants().

    bonds: list of (mu, nu, n, K), e.g. from neighbour_bonds()
    q:     wave vector (1/Å), length d (a number in 1D), in the span of the lattice vectors
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    return _mass_weighted(force_constants(crystal, bonds, q), np.repeat(crystal.masses, crystal.dim))


def dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q), shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return np.array([diagonalize(dynamical_matrix(crystal, bonds, q))[0] for q in qs])


def sound_velocities(crystal, bonds, directions, dq=1e-4):
    """Velocities omega / q of the d acoustic branches for small q along each direction.

    directions: (M, d) directions (normalized here) in the span of the lattice vectors
    dq:         |q| used for the slope (1/Å); relative error ~ (dq a)^2
    Returns (M, d) velocities in Å per internal time unit, ascending per direction
    (transverse before longitudinal); see units.VELOCITY_TO_KM_S.
    """
    if crystal.periodicity == 0:
        raise ValueError("a molecule (no lattice vectors) has no sound velocities")
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    return dispersion(crystal, bonds, dq * directions)[:, : crystal.dim] / dq


def bloch_mode(crystal, bonds, q, branch, cells):
    """Frequency and displacements of one Bloch mode in a supercell of cells^p cells.

    q must fit the supercell, q = (i b1 + j b2 + ...) / cells. For degenerate
    branches the returned mode is one arbitrary combination. Displacements
    (B d cells^p,) start the standing wave Re[e_mu exp(i q.R)] / sqrt(m_mu) for
    simulate(), in the atom order of periodic_supercell(); the cell phase exp(i q.R)
    is the convention of force_constants().
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

    cell_vectors = cell_offsets(cells, crystal.periodicity) @ crystal.vectors
    u0 = np.real(np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None])
    return omega[branch], u0.ravel()


def group_velocities(crystal, bonds, qs):
    """Group velocities d omega / dq of all branches, shape (len(qs), B d, d), in Å per internal time unit.

    From the eigenvectors (Hellmann-Feynman): d omega / dq_i = e^H (dD / dq_i) e / (2 omega).
    The velocity is not defined for zero modes and for degenerate branches (it depends on
    the chosen combination there), so these entries are nan.
    """
    masses = np.repeat(crystal.masses, crystal.dim)
    velocities = []
    for q in qs:
        q = np.atleast_1d(np.asarray(q, dtype=float))
        omega, vectors = diagonalize(dynamical_matrix(crystal, bonds, q))
        gradient = _mass_weighted(_force_constants_gradient(crystal, bonds, q), masses)
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
