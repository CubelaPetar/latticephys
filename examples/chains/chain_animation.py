"""Animate a periodic monatomic chain (last atom coupled to the first); the first atom starts displaced."""

import numpy as np

from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.display import animate
from latticephys.phonons.integrate import simulate
from latticephys.phonons.modes import normal_modes
from latticephys.units import TIME_FS

N_ATOMS = 11
SPRING_CONSTANT = 10.0  # eV/Å²
MASS = 12.011  # amu (carbon)
LATTICE_CONSTANT = 1.5  # Å
AMPLITUDE = 0.3  # Å, exaggerated for visibility

TIMESTEP = 0.05  # internal units (≈ 0.5 fs)
PERIODS = 5  # simulated periods of the slowest mode, played in a loop
SECONDS_PER_PERIOD = 2.0  # playback speed

RED = (188, 39, 50)


def main():
    Phi = chain_force_constants(N_ATOMS, SPRING_CONSTANT, boundary="periodic")
    masses = np.full(N_ATOMS, MASS)

    omega, _ = normal_modes(Phi, masses)
    slowest = omega[omega > 0].min()
    period = 2 * np.pi / slowest
    print(f"slowest mode period: {period * TIME_FS:.1f} fs")
    steps = round(PERIODS * period / TIMESTEP)

    u0 = np.zeros(N_ATOMS)
    u0[0] = AMPLITUDE
    p0 = np.zeros(N_ATOMS)
    u, _ = simulate(Phi, masses, u0, p0, TIMESTEP, steps)

    x0 = LATTICE_CONSTANT * np.arange(N_ATOMS)
    animate(x0, u, TIMESTEP, period / SECONDS_PER_PERIOD, radii=[15] * N_ATOMS, colors=[RED] * N_ATOMS)


if __name__ == "__main__":
    main()
