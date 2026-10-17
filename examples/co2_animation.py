"""Animate a linear m-M-m molecule (CO2-like) starting in the symmetric stretch."""

import numpy as np

from libsolidstate.dynamics.display import animate
from libsolidstate.dynamics.integrate import simulate
from libsolidstate.lattice.chain import chain_force_constants

SPRING_CONSTANT = 1.0
MASS_SMALL = 1.0
MASS_LARGE = 2.0
TIMESTEP = 0.01
STEPS = 6283  # about 10 periods of the symmetric stretch (omega = 1)

BLUE = (100, 149, 237)
RED = (188, 39, 50)


def main():
    K = chain_force_constants(3, SPRING_CONSTANT, boundary="free")
    masses = np.array([MASS_SMALL, MASS_LARGE, MASS_SMALL])

    u0 = np.array([-0.3, 0.0, 0.3])
    p0 = np.zeros(3)
    u, _ = simulate(K, masses, u0, p0, TIMESTEP, STEPS)

    animate([-1.0, 0.0, 1.0], u, radii=[20, 30, 20], colors=[BLUE, RED, BLUE])


if __name__ == "__main__":
    main()
