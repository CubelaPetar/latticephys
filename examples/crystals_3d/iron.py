"""Body-centred cubic iron (α-Fe) with central springs: why bcc needs second neighbours.

One atom per primitive cell, 8 nearest neighbours at a sqrt(3)/2 along the cube diagonals
and 6 second neighbours at a along the cube axes. With nearest-neighbour central springs
only, bcc is floppy: some shears stretch no diagonal spring to first order, and a branch
is omega = 0 along the whole line Γ-N. Springs to the second neighbours make it stable.

Path through the bcc zone: Γ-H-N-Γ-P-H (H at the cube face of the zone, N at an edge
centre, P at the corner where three faces meet). Spring constants are illustrative.
"""

import matplotlib.pyplot as plt

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_bands, plot_crystal
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A = 2.866  # Å, lattice constant of α-Fe
MASS = 55.845  # amu
SPRING_CONSTANTS = [3.0, 1.0]  # eV/Å², first (cube diagonals) and second (cube axes) neighbours
LABELS = ["Γ", "H", "N", "Γ", "P", "H"]


def main():
    iron = Crystal(primitive_vectors("bcc", A), [[0.0, 0.0, 0.0]], [MASS], ["Fe"])
    points = special_points("bcc", iron.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for axis, spring_constants in zip(axes, (SPRING_CONSTANTS[:1], SPRING_CONSTANTS)):
        frequencies = dispersion(iron, neighbour_bonds(iron, spring_constants), qs) * OMEGA_TO_THZ
        plot_bands(distance, frequencies, corners, LABELS, ax=axis)
        axis.set_title(f"bcc Fe, K = {spring_constants} eV/Å²")
    axes[0].set_ylabel(r"$\nu$ (THz)")
    fig.tight_layout()

    ax = plot_crystal(iron, radii=[0.4], colors=["tab:gray"], cells=2)
    ax.set_title("bcc Fe (atoms not to scale) and its primitive cell")
    plt.show()


if __name__ == "__main__":
    main()
