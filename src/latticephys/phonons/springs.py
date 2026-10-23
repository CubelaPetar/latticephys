"""Force constants of springs between atoms of a Crystal.

A bond (mu, nu, n, K) is a spring between atom mu in cell 0 and atom nu in
cell n (see lattice.bonds). Its direction e comes from the crystal. K is either
a number (central spring, force only along e) or a pair (K_L, K_T) with an
additional transverse stiffness K_T against displacements perpendicular to e.

Transverse springs (K_T != 0) are the Born-von Karman model of textbooks. They are
not rotationally invariant: a rigid rotation of a molecule or a ribbon costs
energy, so free systems get no zero rotation modes.
"""

import numpy as np


def force_constants(crystal, bonds):
    """Force-constant matrix Phi of the atoms of one cell, V = 1/2 u^T Phi u.

    Bonds to other cells (n != 0) couple to the periodic images, so for a
    periodic_supercell() this is the finite periodic system; for a molecule
    all bonds have n = 0. Shape (B d, B d), degrees of freedom ordered atom
    by atom: [u_0x, u_0y, u_1x, u_1y, ...]. Use np.repeat(masses, d) as masses.
    """
    return _fourier_force_constants(crystal, bonds, np.zeros(crystal.dim)).real


def _spring_block(r, K):
    """d x d force constants K_L e e^T + K_T (1 - e e^T) of a spring along r, K = K_L or (K_L, K_T)."""
    K_L, K_T = (K, 0.0) if np.ndim(K) == 0 else K
    e = r / np.linalg.norm(r)
    longitudinal = np.outer(e, e)
    return K_L * longitudinal + K_T * (np.eye(len(r)) - longitudinal)


def _fourier_force_constants(crystal, bonds, q, gradient=False):
    """Phi(q) = sum over bonds of the spring blocks with Bloch phases, shape (B d, B d).

    With gradient=True the derivatives dPhi/dq_i instead, shape (d, B d, B d).
    """
    d = crystal.dim
    size = len(crystal.masses) * d
    Phi = np.zeros((d, size, size) if gradient else (size, size), dtype=complex)
    for mu, nu, n, K in bonds:
        block = _spring_block(crystal.bond_vector(mu, nu, n), K)
        cell = crystal.cell_index(n) @ crystal.vectors
        phase = np.exp(1j * q @ cell)

        smu = slice(mu * d, (mu + 1) * d)
        snu = slice(nu * d, (nu + 1) * d)
        if gradient:
            # only the off-diagonal blocks depend on q
            Phi[:, smu, snu] -= 1j * cell[:, None, None] * block * phase
            Phi[:, snu, smu] += 1j * cell[:, None, None] * block * np.conj(phase)
        else:
            Phi[smu, smu] += block
            Phi[snu, snu] += block
            Phi[smu, snu] -= block * phase
            Phi[snu, smu] -= block * np.conj(phase)
    return Phi
