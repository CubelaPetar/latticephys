"""Phonon dispersion of a Crystal with spring bonds, in any dimension."""

import numpy as np

from latticephys.lattice.crystal import cell_offsets
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


def bloch_mode(crystal, bonds, q, branch, cells):
    """Frequency omega and complex pattern U of one Bloch mode in a supercell of cells^p cells.

    U = e_mu exp(i q.R) / sqrt(m_mu), shape (B d cells^p,), in the atom order of
    periodic_supercell(); the cell phase exp(i q.R) is the convention of force_constants().
    The largest component of e is made real, so Re U is never zero. Re[U exp(-i omega t)]
    is the travelling wave (see mode_motion()), Re U cos(omega t) the standing wave.
    q must fit the supercell, q = (i b1 + j b2 + ...) / cells. For degenerate branches
    the mode is one arbitrary combination; zero modes have omega exactly 0.
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    fit = cells * crystal.vectors @ q / (2 * np.pi)
    if not np.allclose(fit, np.round(fit), atol=1e-8):
        raise ValueError("q does not fit the supercell, use q = (i b1 + j b2 + ...) / cells")

    d = crystal.dim
    omega, vectors = phonon_modes(crystal, bonds, [q])
    omega, e = omega[0, branch], vectors[0, :, branch]
    # eigh returns e with an arbitrary complex phase; with e.g. e = (0, i) the real part
    # of U would vanish everywhere, so make the largest component real
    e = e * np.exp(-1j * np.angle(e[np.argmax(np.abs(e))]))
    polarization = (e / np.repeat(np.sqrt(crystal.masses), d)).reshape(-1, d)

    cell_vectors = cell_offsets(cells, crystal.periodicity) @ crystal.vectors
    U = np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None]
    return omega, U.ravel()


def mode_motion(omega, U, times):
    """Exact motion Re[U exp(-i omega t)] of a mode with pattern U, shape (len(times), N).

    U from bloch_mode() gives the travelling wave; U.real gives the standing wave
    Re U cos(omega t). Use it like simulate() output, e.g. for animate().
    """
    return np.real(np.exp(-1j * omega * np.asarray(times, dtype=float))[:, None] * np.asarray(U)[None, :])


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
