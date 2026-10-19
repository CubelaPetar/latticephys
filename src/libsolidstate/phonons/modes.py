"""Frequencies and eigenvectors of dynamical matrices; normal modes of finite systems."""

import numpy as np


def diagonalize(dynamical):
    """Frequencies omega (ascending) and eigenvectors (columns) of a Hermitian dynamical matrix.

    omega = sign(omega^2) sqrt(|omega^2|), so negative values are unstable modes
    with imaginary frequency. Round-off of zero modes (|omega^2| below 1e-10 of
    the largest possible eigenvalue, the maximum absolute row sum) is set to 0.
    """
    omega_squared, vectors = np.linalg.eigh(dynamical)
    omega_squared[np.abs(omega_squared) < 1e-10 * np.abs(dynamical).sum(axis=1).max()] = 0.0
    return np.sign(omega_squared) * np.sqrt(np.abs(omega_squared)), vectors


def normal_modes(force_constants, masses):
    """Normal modes of the dynamical matrix D = M^-1/2 K M^-1/2.

    force_constants: symmetric (N, N) matrix K of the potential V = 1/2 u^T K u
    masses:          (N,) masses

    Returns omega (ascending, see diagonalize()) and the eigenvectors of D as
    columns. The displacement pattern of mode j is masses**-0.5 * vectors[:, j].
    """
    inv_sqrt_m = 1 / np.sqrt(np.asarray(masses, dtype=float))
    return diagonalize(inv_sqrt_m[:, None] * np.asarray(force_constants, dtype=float) * inv_sqrt_m[None, :])
