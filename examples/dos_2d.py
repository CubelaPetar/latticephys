"""In-plane phonon density of states and heat capacity of graphene (spring model).

In 2D the density of states has logarithmic peaks at saddle points of omega(q)
(van Hove singularities, here at M) and steps at band edges. The heat capacity
grows as T^2 at low T (2D Debye law) and approaches 2 k_B per atom (Dulong-Petit
for two in-plane degrees of freedom) at high T.
"""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.reciprocal import q_grid, special_points
from libsolidstate.lattice.springs import neighbour_bonds
from libsolidstate.spectra.dispersion import dispersion
from libsolidstate.spectra.thermal import heat_capacity
from libsolidstate.units import OMEGA_TO_THZ

A = 2.46  # Å
MASS = 12.011  # amu
SPRING_CONSTANTS = [25.0, 5.0]  # eV/Å², first and second neighbours
GRID = 150  # GRID x GRID wave vectors in the Brillouin zone
TEMPERATURES = np.geomspace(5, 3000, 200)  # K


def main():
    honeycomb = [[0.0, 0.0], [A / 2, A / (2 * np.sqrt(3))]]
    crystal = Crystal(primitive_vectors("hexagonal", A), honeycomb, [MASS, MASS], ["C", "C"])
    bonds = neighbour_bonds(crystal, SPRING_CONSTANTS)

    omega = dispersion(crystal, bonds, q_grid(crystal.vectors, GRID))
    at_m = dispersion(crystal, bonds, [special_points("hexagonal", crystal.vectors)["M"]])[0]
    per_atom = crystal.dim * heat_capacity(omega, TEMPERATURES)

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.5))
    left.hist(omega.ravel() * OMEGA_TO_THZ, bins=200, density=True)
    for nu in at_m * OMEGA_TO_THZ:
        left.axvline(nu, color="grey", ls="--", lw=0.8)
    left.set_xlabel(r"$\nu$ (THz)")
    left.set_ylabel("density of states (per THz)")
    left.set_title("DOS; dashed: frequencies at M")

    right.loglog(TEMPERATURES, per_atom, label="graphene (in-plane)")
    reference = TEMPERATURES[10]
    right.loglog(TEMPERATURES[:60], per_atom[10] * (TEMPERATURES[:60] / reference) ** 2, "k:", label=r"$\propto T^2$")
    right.axhline(crystal.dim, color="grey", ls="--", lw=0.8, label="Dulong-Petit")
    right.set_xlabel("T (K)")
    right.set_ylabel(r"$C_V$ per atom ($k_B$)")
    right.legend()

    for axis in (left, right):
        axis.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
