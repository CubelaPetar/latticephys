"""Animate one Bloch mode of the diatomic chain: acoustic or optical branch at wave vector q.

The mode is the exact travelling wave Re[U exp(-i omega t)], U = e_mu exp(i q R) / sqrt(m_mu):
every cell moves alike, shifted in phase by q R, so for 0 < |q| < pi/a the pattern runs
along the chain; at q = 0 and pi/a the phases are ±1 and the wave stands.
In the optical mode at q = 0 the two atoms of each cell move against each other
(the acoustic mode there is a rigid translation and does not oscillate). At the
zone edge q = pi/a only the heavy atoms move on the acoustic branch and only the
light atoms on the optical branch.
"""

import numpy as np

from latticephys.lattice.bonds import periodic_supercell
from latticephys.lattice.crystal import Crystal
from latticephys.phonons.dispersion import bloch_mode, mode_motion
from latticephys.phonons.display import animate
from latticephys.units import OMEGA_TO_THZ

BRANCH = "optical"  # "acoustic" or "optical"
QA = np.pi / 2  # q times a, in [-pi, pi]; must be a multiple of 2 pi / CELLS

CELLS = 8  # periodic chain of CELLS unit cells
SPRING_CONSTANT = 10.0  # eV/Å²
MASS_1 = 12.011  # amu
MASS_2 = 15.999  # amu
LATTICE_CONSTANT = 3.0  # Å, two atoms per cell
AMPLITUDE = 0.2  # Å, exaggerated for visibility

TIMESTEP = 0.02  # internal units, time between frames
PERIODS = 5
SECONDS_PER_PERIOD = 2.0

RED = (188, 39, 50)
BLUE = (100, 149, 237)


def main():
    a = LATTICE_CONSTANT
    crystal = Crystal([[a]], [0.0, a / 2], [MASS_1, MASS_2], ["C", "O"])
    bonds = [(0, 1, 0, SPRING_CONSTANT), (1, 0, 1, SPRING_CONSTANT)]
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, CELLS)
    positions = supercell.positions

    branch = {"acoustic": 0, "optical": 1}[BRANCH]
    omega, U = bloch_mode(crystal, bonds, QA / a, branch, CELLS)
    if omega <= 0:
        raise ValueError("this mode does not oscillate (omega = 0), choose another branch or QA")
    print(f"{BRANCH} mode at qa = {QA:.3f}: nu = {omega * OMEGA_TO_THZ:.2f} THz")

    U *= AMPLITUDE / np.abs(U).max()
    period = 2 * np.pi / omega
    u = mode_motion(omega, U, TIMESTEP * np.arange(round(PERIODS * period / TIMESTEP)))

    animate(positions, u, TIMESTEP, period / SECONDS_PER_PERIOD,
            radii=[12, 15] * CELLS, colors=[RED, BLUE] * CELLS)


if __name__ == "__main__":
    main()
