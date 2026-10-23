"""Classical and quantum thermal properties from phonon frequencies."""

import numpy as np

from latticephys.lattice.reciprocal import q_grid
from latticephys.phonons.dispersion import dynamical_matrix
from latticephys.phonons.modes import diagonalize
from latticephys.units import K_B, OMEGA_TO_MEV


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
    """Classical thermal <|u_mu|^2> (Å²) of each basis atom in a periodic crystal of cells^p cells.

    <|u_mu|^2> = k_B T / N sum over q != rigid translation, branches s of |e_mu,s(q)|^2 / (m_mu omega_s^2),
    from equipartition. The d translations at q = 0 are left out; any other zero mode (floppy, or the
    free rotation of a molecule) or unstable mode raises, since <u^2> diverges there.
    """
    d, n_basis = crystal.dim, len(crystal.masses)
    qs = q_grid(crystal.vectors, cells)
    total = np.zeros(n_basis)
    for q in qs:
        omega, vectors = diagonalize(dynamical_matrix(crystal, bonds, q))
        zero = omega == 0
        translations = 0 if np.any(q != 0) else d
        if np.any(omega < 0) or zero.sum() > translations:
            raise ValueError(f"floppy or unstable mode at q = {q}, <u^2> diverges")

        weight = np.abs(vectors[:, ~zero]) ** 2 / omega[~zero] ** 2
        total += weight.reshape(n_basis, d, -1).sum(axis=(1, 2))
    return K_B * temperature / len(qs) * total / crystal.masses
