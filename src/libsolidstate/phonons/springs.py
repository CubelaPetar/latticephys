"""Force constants of central springs between atoms of a Crystal.

A bond (mu, nu, n, k) is a spring k between atom mu in cell 0 and atom nu in
cell n (see lattice.bonds). Its direction comes from the crystal.
"""

import numpy as np


def force_constants(crystal, bonds):
    """Force-constant matrix K of the atoms of one cell, V = 1/2 u^T K u.

    Bonds to other cells (n != 0) couple to the periodic images, so for a
    periodic_supercell() this is the finite periodic system; for a molecule
    all bonds have n = 0. Shape (B d, B d), degrees of freedom ordered atom
    by atom: [u_0x, u_0y, u_1x, u_1y, ...]. Use np.repeat(masses, d) as masses.
    """
    return _fourier_force_constants(crystal, bonds, np.zeros(crystal.dim)).real


def _fourier_force_constants(crystal, bonds, q):
    """Phi(q) = sum over bonds of the central-spring block k e e^T with Bloch phases."""
    d = crystal.dim
    phi = np.zeros((len(crystal.masses) * d, len(crystal.masses) * d), dtype=complex)
    for mu, nu, n, k in bonds:
        r = crystal.bond_vector(mu, nu, n)
        block = k * np.outer(r, r) / (r @ r)
        phase = np.exp(1j * q @ (np.atleast_1d(n) @ crystal.vectors))

        smu = slice(mu * d, (mu + 1) * d)
        snu = slice(nu * d, (nu + 1) * d)
        phi[smu, smu] += block
        phi[snu, snu] += block
        phi[smu, snu] -= block * phase
        phi[snu, smu] -= block * np.conj(phase)
    return phi
