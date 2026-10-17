"""Force-constant matrix of atoms coupled by central springs."""

import numpy as np


def force_constants(n_atoms, bonds, dim=1):
    """Force-constant matrix K of a finite system, V = 1/2 u^T K u.

    bonds: list of (i, j, k, r): spring k between atoms i and j along the
           bond vector r (scalar in 1D, length-dim array otherwise; only
           its direction is used). Each spring only resists stretching.
    dim:   number of displacement components per atom

    Returns K of shape (n_atoms * dim, n_atoms * dim), degrees of freedom
    ordered atom by atom: [u_0x, u_0y, u_1x, u_1y, ...].
    Use np.repeat(masses, dim) as the matching masses.
    """
    K = np.zeros((n_atoms * dim, n_atoms * dim))
    for i, j, k, r in bonds:
        e = np.atleast_1d(np.asarray(r, dtype=float))
        e = e / np.linalg.norm(e)
        block = k * np.outer(e, e)

        si = slice(i * dim, (i + 1) * dim)
        sj = slice(j * dim, (j + 1) * dim)
        K[si, si] += block
        K[sj, sj] += block
        K[si, sj] -= block
        K[sj, si] -= block
    return K


def neighbour_bonds(crystal, spring_constants, search=3):
    """Periodic bonds of a crystal to its neighbour shells.

    spring_constants: spring constant per shell, [k1] for nearest neighbours,
                      [k1, k2] to add second neighbours, ...
    search:           cells searched in each direction

    Shells are the distinct atom distances in ascending order. Returns a list
    of (mu, nu, n, k, r): spring k between atom mu in cell 0 and atom nu in
    cell n (integer array), r the bond vector. Each bond appears once.
    """
    steps = range(-search, search + 1)
    candidates = []
    for mu, r_mu in enumerate(crystal.positions):
        for nu, r_nu in enumerate(crystal.positions):
            for n in ((i, j) for i in steps for j in steps):
                # count each bond only from one end: (mu, 0) -- (nu, n) equals (nu, 0) -- (mu, -n)
                if nu < mu or (nu == mu and n <= (0, 0)):
                    continue
                r = r_nu + np.array(n) @ crystal.vectors - r_mu
                candidates.append((mu, nu, np.array(n), np.linalg.norm(r), r))

    distances = np.array([c[3] for c in candidates])
    shells = np.unique(np.round(distances, 6))[: len(spring_constants)]

    bonds = []
    for mu, nu, n, distance, r in candidates:
        shell = np.flatnonzero(np.isclose(shells, distance, atol=1e-6))
        if len(shell):
            bonds.append((mu, nu, n, spring_constants[shell[0]], r))
    return bonds
