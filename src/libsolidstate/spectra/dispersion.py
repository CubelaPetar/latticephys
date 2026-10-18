"""Phonon dispersion of a Crystal with spring bonds, in any dimension."""

import numpy as np

from libsolidstate.lattice.crystal import cell_offsets
from libsolidstate.lattice.springs import _fourier_force_constants


def dynamical_matrix(crystal, bonds, q):
    """Bloch dynamical matrix D(q) = M^-1/2 Phi(q) M^-1/2, shape (B d, B d).

    bonds: list of (mu, nu, n, k), e.g. from neighbour_bonds()
    q:     wave vector (1/Å), length d (a number in 1D)
    Degrees of freedom are ordered atom by atom: [u_0x, u_0y, u_1x, ...].
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    if q.shape != (crystal.dim,):
        raise ValueError(f"q must have {crystal.dim} components, got {q.shape}")

    inv_sqrt_m = np.repeat(1 / np.sqrt(crystal.masses), crystal.dim)
    return inv_sqrt_m[:, None] * _fourier_force_constants(crystal, bonds, q) * inv_sqrt_m[None, :]


def dispersion(crystal, bonds, qs):
    """Branch frequencies omega(q), shape (len(qs), B d), ascending per q.

    Negative values are unstable modes with imaginary frequency.
    """
    omega_squared = np.array([np.linalg.eigvalsh(dynamical_matrix(crystal, bonds, q)) for q in qs])
    return _frequencies(omega_squared, _zero_tolerance(crystal, bonds))


def sound_velocities(crystal, bonds, directions, dq=1e-4):
    """Velocities omega / q of the d acoustic branches for small q along each direction.

    directions: (M, d) directions (normalized here)
    dq:         |q| used for the slope (1/Å); relative error ~ (dq a)^2
    Returns (M, d) velocities in Å per internal time unit, ascending per direction
    (transverse before longitudinal); see units.VELOCITY_TO_KM_S.
    """
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    return dispersion(crystal, bonds, dq * directions)[:, : crystal.dim] / dq


def bloch_mode(crystal, bonds, q, branch, cells):
    """Frequency and displacements of one Bloch mode in a supercell of cells^d cells.

    q must fit the supercell, q = (i b1 + j b2 + ...) / cells. For degenerate
    branches the returned mode is one arbitrary combination. Displacements
    (B d cells^d,) start the standing wave Re[e_mu exp(i q.R)] / sqrt(m_mu) for
    simulate(), in the atom order of periodic_supercell().
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    fit = cells * crystal.vectors @ q / (2 * np.pi)
    if not np.allclose(fit, np.round(fit), atol=1e-8):
        raise ValueError("q does not fit the supercell, use q = (i b1 + j b2 + ...) / cells")

    d = crystal.dim
    omega_squared, vectors = np.linalg.eigh(dynamical_matrix(crystal, bonds, q))
    e = vectors[:, branch]
    # eigh returns e with an arbitrary complex phase; with e.g. e = (0, i) the real part
    # below would vanish everywhere, so make the largest component real
    e = e * np.exp(-1j * np.angle(e[np.argmax(np.abs(e))]))
    polarization = (e / np.repeat(np.sqrt(crystal.masses), d)).reshape(-1, d)

    cell_vectors = cell_offsets(cells, d) @ crystal.vectors
    u0 = np.real(np.exp(1j * cell_vectors @ q)[:, None, None] * polarization[None])
    return _frequencies(omega_squared, _zero_tolerance(crystal, bonds))[branch], u0.ravel()


def _zero_tolerance(crystal, bonds):
    """|omega^2| below this is round-off of a zero mode; relative to the scale sum |k| / m_min of D."""
    return 1e-10 * sum(abs(k) for *_, k in bonds) / crystal.masses.min()


def _frequencies(omega_squared, tolerance):
    """omega = sign(omega^2) sqrt(|omega^2|), with round-off below tolerance set to 0."""
    omega_squared = np.array(omega_squared, dtype=float)
    omega_squared[np.abs(omega_squared) < tolerance] = 0.0
    return np.sign(omega_squared) * np.sqrt(np.abs(omega_squared))
