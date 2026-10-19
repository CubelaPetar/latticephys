"""Animate one Bloch mode of the diatomic chain: acoustic or optical branch at wave vector q.

In the optical mode at q = 0 the two atoms of each cell move against each other
(the acoustic mode there is a rigid translation and does not oscillate). At the
zone edge q = pi/a only the heavy atoms move on the acoustic branch and only the
light atoms on the optical branch.
"""

import numpy as np

from latticephys.lattice.bonds import periodic_supercell
from latticephys.lattice.crystal import Crystal
from latticephys.phonons.dispersion import bloch_mode
from latticephys.phonons.display import animate
from latticephys.phonons.integrate import simulate
from latticephys.phonons.springs import force_constants
from latticephys.units import OMEGA_TO_THZ

BRANCH = "optical"  # "acoustic" or "optical"
QA = np.pi / 2  # q times a, in [-pi, pi]; must be a multiple of 2 pi / CELLS

CELLS = 8  # periodic chain of CELLS unit cells
SPRING_CONSTANT = 10.0  # eV/Å²
MASS_1 = 12.011  # amu
MASS_2 = 15.999  # amu
LATTICE_CONSTANT = 3.0  # Å, two atoms per cell
AMPLITUDE = 0.2  # Å, exaggerated for visibility

TIMESTEP = 0.02  # internal units
PERIODS = 5
SECONDS_PER_PERIOD = 2.0

RED = (188, 39, 50)
BLUE = (100, 149, 237)


def main():
    a = LATTICE_CONSTANT
    crystal = Crystal([[a]], [0.0, a / 2], [MASS_1, MASS_2], ["C", "O"])
    bonds = [(0, 1, 0, SPRING_CONSTANT), (1, 0, 1, SPRING_CONSTANT)]
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, CELLS)
    K = force_constants(supercell, supercell_bonds)
    positions, masses = supercell.positions, supercell.masses

    branch = {"acoustic": 0, "optical": 1}[BRANCH]
    omega, u0 = bloch_mode(crystal, bonds, QA / a, branch, CELLS)
    if omega <= 0:
        raise ValueError("this mode does not oscillate (omega = 0), choose another branch or QA")
    print(f"{BRANCH} mode at qa = {QA:.3f}: nu = {omega * OMEGA_TO_THZ:.2f} THz")

    u0 *= AMPLITUDE / np.abs(u0).max()
    period = 2 * np.pi / omega
    steps = round(PERIODS * period / TIMESTEP)
    u, _ = simulate(K, np.repeat(masses, crystal.dim), u0, np.zeros_like(u0), TIMESTEP, steps)

    centred = positions - positions.mean(axis=0)
    animate(centred, u, TIMESTEP, period / SECONDS_PER_PERIOD,
            radii=[12, 15] * CELLS, colors=[RED, BLUE] * CELLS, scale=30)


if __name__ == "__main__":
    main()
