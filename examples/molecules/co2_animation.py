"""Animate CO2 (O=C=O) starting in the symmetric stretch."""

import numpy as np

from libsolidstate.phonons.chain import chain_force_constants
from libsolidstate.phonons.display import animate
from libsolidstate.phonons.integrate import simulate
from libsolidstate.units import TIME_FS

SPRING_CONSTANT = 100.0  # eV/Å²
MASS_O = 15.999  # amu
MASS_C = 12.011  # amu
BOND_LENGTH = 1.16  # Å
AMPLITUDE = 0.1  # Å, exaggerated for visibility

TIMESTEP = 0.01  # internal units (≈ 0.1 fs)
PERIODS = 10  # simulated periods of the symmetric stretch, played in a loop
SECONDS_PER_PERIOD = 2.0  # playback speed

BLUE = (100, 149, 237)
RED = (188, 39, 50)


def main():
    K = chain_force_constants(3, SPRING_CONSTANT, boundary="free")
    masses = np.array([MASS_O, MASS_C, MASS_O])

    period = 2 * np.pi / np.sqrt(SPRING_CONSTANT / MASS_O)
    print(f"symmetric stretch period: {period * TIME_FS:.1f} fs")
    steps = round(PERIODS * period / TIMESTEP)

    u0 = AMPLITUDE * np.array([-1.0, 0.0, 1.0])
    p0 = np.zeros(3)
    u, _ = simulate(K, masses, u0, p0, TIMESTEP, steps)

    x0 = BOND_LENGTH * np.array([-1.0, 0.0, 1.0])
    animate(x0, u, TIMESTEP, period / SECONDS_PER_PERIOD, radii=[20, 15, 20], colors=[RED, BLUE, RED], scale=200)


if __name__ == "__main__":
    main()
