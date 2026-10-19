"""Animate one in-plane Bloch phonon of graphene in a periodic supercell.

Choose the k-point (Γ, M or K) and the branch (0, 1 acoustic; 2, 3 optical,
ascending in frequency). At Γ the optical modes move the two sublattices against
each other; degenerate modes (Γ optical, K) show one arbitrary combination.
"""

import numpy as np

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.reciprocal import special_points
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import bloch_mode
from latticephys.phonons.display import animate
from latticephys.phonons.integrate import simulate
from latticephys.phonons.springs import force_constants
from latticephys.units import OMEGA_TO_THZ

POINT = "M"  # "Γ", "M" or "K"
BRANCH = 3  # 0 .. 3

A = 2.46  # Å
MASS = 12.011  # amu
SPRING_CONSTANTS = [25.0, 5.0]  # eV/Å², first and second neighbours
CELLS = 6  # supercell of CELLS x CELLS cells; multiple of 6 fits Γ, M and K
AMPLITUDE = 0.25  # Å, exaggerated for visibility

TIMESTEP = 0.01  # internal units
PERIODS = 5
SECONDS_PER_PERIOD = 2.0

RED = (188, 39, 50)
BLUE = (100, 149, 237)


def main():
    crystal = honeycomb(A, [MASS, MASS], ["C", "C"])
    bonds = neighbour_bonds(crystal, SPRING_CONSTANTS)
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, CELLS)
    K = force_constants(supercell, supercell_bonds)
    positions, masses = supercell.positions, supercell.masses

    q = special_points("hexagonal", crystal.vectors)[POINT]
    omega, u0 = bloch_mode(crystal, bonds, q, BRANCH, CELLS)
    if omega <= 0:
        raise ValueError("this mode does not oscillate (omega = 0), choose another branch or k-point")
    print(f"{POINT}, branch {BRANCH}: nu = {omega * OMEGA_TO_THZ:.2f} THz")

    u0 *= AMPLITUDE / np.abs(u0).max()
    period = 2 * np.pi / omega
    steps = round(PERIODS * period / TIMESTEP)
    u, _ = simulate(K, np.repeat(masses, 2), u0, np.zeros_like(u0), TIMESTEP, steps)

    centred = positions - positions.mean(axis=0)
    scale = 350 / np.abs(centred).max()
    animate(centred, u, TIMESTEP, period / SECONDS_PER_PERIOD,
            radii=[6] * len(positions), colors=[RED, BLUE] * (len(positions) // 2), scale=scale)


if __name__ == "__main__":
    main()
