"""Classical and quantum thermal properties from phonon frequencies."""

import numpy as np

from libsolidstate.lattice.reciprocal import q_grid
from libsolidstate.spectra.dispersion import _zero_tolerance, dynamical_matrix
from libsolidstate.units import K_B, OMEGA_TO_MEV


def heat_capacity(omega, temperatures):
    """Heat capacity per mode in units of k_B (Bose statistics), averaged over the frequencies omega.

    omega:        angular frequencies (internal units), e.g. dispersion() on a q_grid(); omega = 0
                  counts with its classical limit 1
    temperatures: K
    Multiply by d (degrees of freedom per atom) for the heat capacity per atom.
    """
    omega = np.ravel(omega)
    if np.any(omega < 0):
        raise ValueError("unstable modes (omega < 0) have no heat capacity")

    x = (omega * OMEGA_TO_MEV / 1000)[None, :] / (K_B * np.atleast_1d(temperatures)[:, None])
    safe = np.where(x > 0, x, 1.0)
    per_mode = np.where(x > 0, safe**2 * np.exp(-safe) / np.expm1(-safe) ** 2, 1.0)
    return per_mode.mean(axis=1)


def mean_square_displacement(crystal, bonds, cells, temperature):
    """Classical thermal <|u_mu|^2> (Å²) of each basis atom in a periodic crystal of cells^d cells.

    <|u_mu|^2> = k_B T / N sum over q != rigid translation, branches s of |e_mu,s(q)|^2 / (m_mu omega_s^2),
    from equipartition. Raises for floppy or unstable modes, where it diverges.
    """
    d, n_basis = crystal.dim, len(crystal.masses)
    qs = q_grid(crystal.vectors, cells)
    tolerance = _zero_tolerance(crystal, bonds)
    total = np.zeros(n_basis)
    for q in qs:
        omega_squared, vectors = np.linalg.eigh(dynamical_matrix(crystal, bonds, q))
        zero = np.abs(omega_squared) < tolerance
        if np.any(omega_squared[~zero] < 0) or (np.any(zero) and np.any(q != 0)):
            raise ValueError(f"floppy or unstable mode at q = {q}, <u^2> diverges")

        weight = np.abs(vectors[:, ~zero]) ** 2 / omega_squared[~zero]
        total += weight.reshape(n_basis, d, -1).sum(axis=(1, 2))
    return K_B * temperature / len(qs) * total / crystal.masses
