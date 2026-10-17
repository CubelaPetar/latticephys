"""Animate a monatomic chain between two fixed walls; the first atom starts displaced."""

import numpy as np

from libsolidstate.dynamics.display import animate
from libsolidstate.dynamics.integrate import simulate
from libsolidstate.lattice.chain import chain_force_constants

N_ATOMS = 11
SPRING_CONSTANT = 1.0
MASS = 1.0
LATTICE_CONSTANT = 1.0
TIMESTEP = 0.01
STEPS = 5000

RED = (188, 39, 50)


def main():
    K = chain_force_constants(N_ATOMS, SPRING_CONSTANT, boundary="fixed")
    masses = np.full(N_ATOMS, MASS)

    u0 = np.zeros(N_ATOMS)
    u0[0] = 0.3
    p0 = np.zeros(N_ATOMS)
    u, _ = simulate(K, masses, u0, p0, TIMESTEP, STEPS)

    x0 = LATTICE_CONSTANT * (np.arange(N_ATOMS) - (N_ATOMS - 1) / 2)
    animate(x0, u, radii=[15] * N_ATOMS, colors=[RED] * N_ATOMS, scale=60)


if __name__ == "__main__":
    main()
