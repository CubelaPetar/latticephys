"""Time integration of the harmonic equations of motion m u'' = -K u."""

import numpy as np


def simulate(force_constants, masses, u0, p0, dt, steps):
    """Integrate with velocity Verlet (symplectic, all atoms updated together).

    force_constants: (N, N) matrix K
    masses:          (N,) masses, one per degree of freedom (np.repeat(masses, 2) in 2D)
    u0, p0:          initial displacements and momenta
    Returns displacements u and momenta p, both of shape (steps + 1, N).
    """
    m = np.asarray(masses, dtype=float)
    K = np.asarray(force_constants, dtype=float)

    u = np.empty((steps + 1, len(m)))
    p = np.empty((steps + 1, len(m)))
    u[0] = u0
    p[0] = p0

    force = -K @ u[0]
    for t in range(steps):
        p_half = p[t] + 0.5 * dt * force
        u[t + 1] = u[t] + dt * p_half / m
        force = -K @ u[t + 1]
        p[t + 1] = p_half + 0.5 * dt * force

    return u, p
