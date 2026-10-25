"""Phonon dispersion omega(q) of a Crystal and its slopes: group and sound velocities, in any dimension."""

import numpy as np

from latticephys.phonons.modes import _mass_weighted, phonon_modes
from latticephys.phonons.springs import _force_constants_gradient


def dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q), shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    return phonon_modes(crystal, bonds, qs)[0]


def sound_velocities(crystal, bonds, directions, dq=1e-4):
    """Velocities omega / q of the d acoustic branches for small q along each direction.

    directions: (M, d) directions (normalized here) in the span of the lattice vectors
    dq:         |q| used for the slope (1/Å); relative error ~ (dq a)^2
    Returns (M, d) velocities in Å per internal time unit, ascending per direction
    (transverse before longitudinal); see units.VELOCITY_TO_KM_S.
    """
    if crystal.periodicity == 0:
        raise ValueError("a molecule (no lattice vectors) has no sound velocities")
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    return dispersion(crystal, bonds, dq * directions)[:, : crystal.dim] / dq


def group_velocities(crystal, bonds, qs):
    """Group velocities d omega / dq of all branches, shape (len(qs), B d, d), in Å per internal time unit.

    From the eigenvectors (Hellmann-Feynman): d omega / dq_i = e^H (dD / dq_i) e / (2 omega).
    The velocity is not defined for zero modes and for degenerate branches (it depends on
    the chosen combination there), so these entries are nan.
    """
    masses = np.repeat(crystal.masses, crystal.dim)
    velocities = []
    for q, omega, vectors in zip(qs, *phonon_modes(crystal, bonds, qs)):
        q = np.atleast_1d(np.asarray(q, dtype=float))
        gradient = _mass_weighted(_force_constants_gradient(crystal, bonds, q), masses)
        slope = np.einsum("as,iab,bs->si", vectors.conj(), gradient, vectors).real

        gap = np.abs(omega[:, None] - omega[None, :])
        np.fill_diagonal(gap, np.inf)
        undefined = (omega == 0) | (gap.min(axis=1) < 1e-6 * np.abs(omega).max())
        with np.errstate(divide="ignore", invalid="ignore"):
            velocity = slope / (2 * omega[:, None])
        velocity[undefined] = np.nan
        velocities.append(velocity)
    return np.array(velocities)
