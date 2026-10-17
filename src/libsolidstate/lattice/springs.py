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
