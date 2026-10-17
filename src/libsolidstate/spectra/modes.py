"""Normal modes of a finite system of masses coupled by harmonic springs."""

import numpy as np


def normal_modes(force_constants, masses):
    """Diagonalize the dynamical matrix D = M^-1/2 K M^-1/2.

    force_constants: symmetric (N, N) matrix K of the potential V = 1/2 u^T K u
    masses:          (N,) masses

    Returns omega^2 (ascending) and the eigenvectors of D as columns.
    The displacement pattern of mode j is masses**-0.5 * vectors[:, j].
    """
    inv_sqrt_m = 1 / np.sqrt(np.asarray(masses, dtype=float))
    dynamical = inv_sqrt_m[:, None] * np.asarray(force_constants) * inv_sqrt_m[None, :]
    return np.linalg.eigh(dynamical)
