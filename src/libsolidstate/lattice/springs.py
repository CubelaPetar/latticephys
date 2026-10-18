"""Central springs between atoms of a Crystal: force constants, neighbour bonds and periodic supercells.

A bond (mu, nu, n, k) is a spring k between atom mu in cell 0 and atom nu in
cell n (integer tuple, or int in 1D). Its direction comes from the crystal.
"""

import numpy as np

from libsolidstate.lattice.crystal import Crystal, cell_offsets


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


def neighbour_bonds(crystal, spring_constants, search=3):
    """Bonds of a crystal to its neighbour shells.

    spring_constants: spring constant per shell, [k1] for nearest neighbours,
                      [k1, k2] to add second neighbours, ...
    search:           cells searched in each direction

    Shells are the distinct atom distances in ascending order. Returns a list
    of bonds (mu, nu, n, k), each bond once.
    """
    origin = (0,) * crystal.dim
    candidates = []
    for mu in range(len(crystal.masses)):
        for nu in range(len(crystal.masses)):
            for n in map(tuple, cell_offsets(2 * search + 1, crystal.dim) - search):
                # count each bond only from one end: (mu, 0) -- (nu, n) equals (nu, 0) -- (mu, -n)
                if nu < mu or (nu == mu and n <= origin):
                    continue
                candidates.append((mu, nu, n, np.linalg.norm(crystal.bond_vector(mu, nu, n))))

    shells = np.unique(np.round([c[3] for c in candidates], 6))[: len(spring_constants)]

    bonds = []
    for mu, nu, n, distance in candidates:
        shell = np.flatnonzero(np.isclose(shells, distance, atol=1e-6))
        if len(shell):
            bonds.append((mu, nu, n, spring_constants[shell[0]]))
    return bonds


def periodic_supercell(crystal, bonds, cells):
    """Supercell of cells^d unit cells as a Crystal, with its bonds.

    force_constants() of the result is the finite system with periodic
    boundaries. Atom mu of cell c (integer tuple, 0 <= c_i < cells) gets
    index ravel(c) * B + mu, with the last entry of c running fastest.
    """
    d, n_basis = crystal.dim, len(crystal.masses)
    shape = (cells,) * d
    offsets = cell_offsets(cells, d)

    supercell = Crystal(
        cells * crystal.vectors,
        (offsets @ crystal.vectors)[:, None, :] + crystal.positions[None, :, :],
        np.tile(crystal.masses, len(offsets)),
        list(crystal.names) * len(offsets),
    )

    supercell_bonds = []
    for c in offsets:
        for mu, nu, n, k in bonds:
            target = c + np.atleast_1d(n)
            i = np.ravel_multi_index(c, shape) * n_basis + mu
            j = np.ravel_multi_index(target % cells, shape) * n_basis + nu
            supercell_bonds.append((i, j, tuple(target // cells), k))
    return supercell, supercell_bonds
