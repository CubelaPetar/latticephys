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

from latticephys.lattice.reciprocal import lattice_projection


def force_constants(crystal, bonds, q=None):
    """Force-constant matrix Phi, V = 1/2 u^T Phi u, or its lattice Fourier transform Phi(q).

    Without q: the real Phi of the atoms of one cell. Bonds to other cells (n != 0)
    couple to the periodic images, so for a periodic_supercell() this is the finite
    periodic system; for a molecule all bonds have n = ().

    With q (1/Å, length d, a number in 1D, in the span of the lattice vectors): the
    Hermitian Phi(q)_{mu nu} = sum_n Phi_{mu 0, nu n} exp(i q.R_n) with cell phases,
    R_n the lattice vector of cell n (not the atom positions). Then Phi(q + G) = Phi(q)
    for every reciprocal lattice vector G, and a Bloch mode moves atom nu in cell n
    as e_nu exp(i q.R_n).

    Shape (B d, B d), degrees of freedom ordered atom by atom: [u_0x, u_0y, u_1x, u_1y, ...].
    Use np.repeat(masses, d) as masses.
    """
    blocks = _spring_blocks(crystal, bonds)
    if q is None:
        return (_lattice_fourier(crystal, _couplings(blocks), np.zeros(crystal.dim)) + _on_site(crystal, blocks)).real
    return _lattice_fourier(crystal, _couplings(blocks), _check_q(crystal, q)) + _on_site(crystal, blocks)


def _force_constants_gradient(crystal, bonds, q):
    """dPhi(q) / dq_i, shape (d, B d, B d): the transform of the couplings times i (R_n)_i.

    The on-site term does not depend on q.
    """
    couplings = _couplings(_spring_blocks(crystal, bonds))
    gradient = []
    for i in range(crystal.dim):
        scaled = [(mu, nu, n, 1j * (crystal.cell_index(n) @ crystal.vectors)[i] * c) for mu, nu, n, c in couplings]
        gradient.append(_lattice_fourier(crystal, scaled, q))
    return np.array(gradient)


def _check_q(crystal, q):
    """q as a float array of d components; raises if it has another length or leaves the span of the lattice vectors."""
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.shape != (crystal.dim,):
        raise ValueError(f"q must have {crystal.dim} components, got {q.shape}")
    if np.linalg.norm(q - lattice_projection(crystal.vectors, q)) > 1e-8 * np.linalg.norm(q):
        raise ValueError(f"q = {q} leaves the span of the lattice vectors")
    return q


def _spring_block(r, K):
    """d x d force constants K_L e e^T + K_T (1 - e e^T) of a spring along r, K = K_L or (K_L, K_T)."""
    K_L, K_T = (K, 0.0) if np.ndim(K) == 0 else K
    e = r / np.linalg.norm(r)
    longitudinal = np.outer(e, e)
    return K_L * longitudinal + K_T * (np.eye(len(r)) - longitudinal)


def _spring_blocks(crystal, bonds):
    """Bonds (mu, nu, n, block) with the spring block of each bond."""
    return [(mu, nu, n, _spring_block(crystal.bond_vector(mu, nu, n), K)) for mu, nu, n, K in bonds]


def _couplings(blocks):
    """Bonds (mu, nu, n, coupling): the coupling Phi_{mu 0, nu n} of a spring is minus its block."""
    return [(mu, nu, n, -block) for mu, nu, n, block in blocks]


def _on_site(crystal, blocks):
    """Sum-rule term: each spring block on the diagonal of both of its atoms, so a rigid shift costs no energy."""
    d = crystal.dim
    size = len(crystal.masses) * d
    Phi = np.zeros((size, size))
    for mu, nu, _, block in blocks:
        Phi[mu * d : (mu + 1) * d, mu * d : (mu + 1) * d] += block
        Phi[nu * d : (nu + 1) * d, nu * d : (nu + 1) * d] += block
    return Phi


def _lattice_fourier(crystal, couplings, q):
    """Sum over bonds (mu, nu, n, c) of c exp(i q.R_n) at (mu, nu) and its Hermitian conjugate at (nu, mu).

    Only the terms between bonded atoms, with cell phases; shape (B d, B d).
    """
    d = crystal.dim
    size = len(crystal.masses) * d
    Phi = np.zeros((size, size), dtype=complex)
    for mu, nu, n, c in couplings:
        phase = np.exp(1j * q @ (crystal.cell_index(n) @ crystal.vectors))
        Phi[mu * d : (mu + 1) * d, nu * d : (nu + 1) * d] += c * phase
        Phi[nu * d : (nu + 1) * d, mu * d : (mu + 1) * d] += c.conj().T * np.conj(phase)
    return Phi
