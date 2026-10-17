"""Phonon dispersion of lattices with a basis (1D chains and 2D crystals)."""

import numpy as np


def dynamical_matrix(q, masses, springs, a):
    """Bloch dynamical matrix D(q) of a 1D lattice with lattice constant a.

    masses:  (B,) masses of the basis atoms in one unit cell
    springs: list of (mu, nu, n, k): spring k between atom mu in cell 0
             and atom nu in cell n
    """
    bonds = [(mu, nu, [n], k, [1.0]) for mu, nu, n, k in springs]
    return _dynamical_matrix(np.atleast_1d(q), masses, bonds, [[a]])


def dispersion(qs, masses, springs, a):
    """Branch frequencies omega(q) of a 1D lattice, shape (len(qs), B), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return _frequencies([dynamical_matrix(q, masses, springs, a) for q in qs])


def crystal_dynamical_matrix(q, crystal, bonds):
    """Bloch dynamical matrix D(q) of a crystal, shape (B d, B d).

    q:     wave vector (1/Å)
    bonds: list of (mu, nu, n, k, r) as returned by neighbour_bonds()
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    return _dynamical_matrix(np.asarray(q, dtype=float), crystal.masses, bonds, crystal.vectors)


def crystal_dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q) of a crystal, shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return _frequencies([crystal_dynamical_matrix(q, crystal, bonds) for q in qs])


def bloch_mode(crystal, bonds, q, branch, cells):
    """Frequency and displacements of one Bloch mode in a cells x cells supercell.

    q must fit the supercell, q = (i b1 + j b2) / cells. For degenerate branches
    the returned mode is one arbitrary combination. Displacements (2 M,) start the
    standing wave Re[e_mu exp(i q.R)] / sqrt(m_mu) for simulate(), in the atom
    order of periodic_supercell().
    """
    fit = cells * crystal.vectors @ q / (2 * np.pi)
    if not np.allclose(fit, np.round(fit), atol=1e-8):
        raise ValueError("q does not fit the supercell, use q = (i b1 + j b2) / cells")

    omega_squared, vectors = np.linalg.eigh(crystal_dynamical_matrix(q, crystal, bonds))
    polarization = (vectors[:, branch] / np.repeat(np.sqrt(crystal.masses), 2)).reshape(-1, 2)

    cell_vectors = np.array([np.array([c1, c2]) @ crystal.vectors for c1 in range(cells) for c2 in range(cells)])
    u0 = np.real(np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None])
    omega = np.sign(omega_squared[branch]) * np.sqrt(np.abs(omega_squared[branch]))
    return omega, u0.ravel()


def _dynamical_matrix(q, masses, bonds, vectors):
    """D(q) in d dimensions; each bond adds the block k e e^T of a central spring."""
    m = np.asarray(masses, dtype=float)
    vectors = np.asarray(vectors, dtype=float)
    dim = len(vectors)

    D = np.zeros((len(m) * dim, len(m) * dim), dtype=complex)
    for mu, nu, n, k, r in bonds:
        e = np.atleast_1d(np.asarray(r, dtype=float))
        e = e / np.linalg.norm(e)
        block = k * np.outer(e, e)
        phase = np.exp(1j * q @ (np.asarray(n) @ vectors))

        smu = slice(mu * dim, (mu + 1) * dim)
        snu = slice(nu * dim, (nu + 1) * dim)
        D[smu, smu] += block / m[mu]
        D[snu, snu] += block / m[nu]
        D[smu, snu] -= block * phase / np.sqrt(m[mu] * m[nu])
        D[snu, smu] -= block * np.conj(phase) / np.sqrt(m[mu] * m[nu])
    return D


def _frequencies(matrices):
    """omega = sign(omega^2) sqrt(|omega^2|), with round-off around zero set to 0."""
    omega_squared = np.array([np.linalg.eigvalsh(D) for D in matrices])
    omega_squared[np.abs(omega_squared) < 1e-10 * np.abs(omega_squared).max()] = 0.0
    return np.sign(omega_squared) * np.sqrt(np.abs(omega_squared))
