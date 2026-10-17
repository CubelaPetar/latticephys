"""Phonon dispersion of a 1D lattice with a basis."""

import numpy as np


def dynamical_matrix(q, masses, springs, a):
    """Bloch dynamical matrix D(q) of a 1D lattice with lattice constant a.

    masses:  (B,) masses of the basis atoms in one unit cell
    springs: list of (mu, nu, n, k): spring k between atom mu in cell 0
             and atom nu in cell n
    """
    m = np.asarray(masses, dtype=float)
    D = np.zeros((len(m), len(m)), dtype=complex)
    for mu, nu, n, k in springs:
        phase = np.exp(1j * q * n * a)
        coupling = k / np.sqrt(m[mu] * m[nu])
        D[mu, mu] += k / m[mu]
        D[nu, nu] += k / m[nu]
        D[mu, nu] -= coupling * phase
        D[nu, mu] -= coupling * np.conj(phase)
    return D


def dispersion(qs, masses, springs, a):
    """Branch frequencies omega(q), shape (len(qs), B), ascending per q."""
    omega_squared = np.array([np.linalg.eigvalsh(dynamical_matrix(q, masses, springs, a)) for q in qs])
    # remove round-off negatives of the acoustic branch at q = 0
    return np.sqrt(np.clip(omega_squared, 0, None))
