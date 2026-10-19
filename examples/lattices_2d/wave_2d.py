"""Displace one atom of a triangular lattice along x and watch the waves spread.

Along x the disturbance travels as a longitudinal wave, along y as a slower
transverse wave, so two fronts with different speeds appear. Atoms are coloured by
their displacement, so the fronts stay visible where the motion is small. The lattice is
periodic, so the waves re-enter from the opposite side.
"""

import numpy as np

from libsolidstate.lattice.bonds import neighbour_bonds, periodic_supercell
from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.phonons.display import animate
from libsolidstate.phonons.integrate import simulate
from libsolidstate.phonons.springs import force_constants
from libsolidstate.units import TIME_FS

A = 2.46  # Å
MASS = 12.011  # amu
SPRING_CONSTANT = 10.0  # eV/Å², nearest neighbours
CELLS = 24  # supercell of CELLS x CELLS atoms
KICK = 0.1  # Å, initial x displacement of the centre atom
EXAGGERATION = 25.0  # displacements are drawn this many times larger
COLOR_SCALE = 0.3  # Å of drawn displacement shown in full colour (yellow); atoms at rest are dark

TIMESTEP = 0.02  # internal units
DURATION = 15.0  # internal units, played in a loop
SECONDS = 8.0  # playback time of DURATION

def main():
    crystal = Crystal(primitive_vectors("hexagonal", A), [[0.0, 0.0]], [MASS], ["C"])
    bonds = neighbour_bonds(crystal, [SPRING_CONSTANT])
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, CELLS)
    K = force_constants(supercell, supercell_bonds)
    positions, masses = supercell.positions, supercell.masses

    centred = positions - positions.mean(axis=0)
    centre = np.argmin(np.linalg.norm(centred, axis=1))
    u0 = np.zeros(2 * len(positions))
    u0[2 * centre] = KICK

    print(f"simulating {DURATION * TIME_FS:.0f} fs")
    steps = round(DURATION / TIMESTEP)
    u, _ = simulate(K, np.repeat(masses, 2), u0, np.zeros_like(u0), TIMESTEP, steps)

    scale = 350 / np.abs(centred).max()
    animate(centred, EXAGGERATION * u, TIMESTEP, DURATION / SECONDS,
            radii=[4] * len(positions), colors=None, scale=scale, color_scale=COLOR_SCALE)


if __name__ == "__main__":
    main()
