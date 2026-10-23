"""Hexagonal close-packed magnesium with central springs: a hexagonal lattice with two atoms.

hcp is the hexagonal Bravais lattice with a two-atom basis: layers A and B stacked ABAB
along c. Two atoms per primitive cell give 3 acoustic and 3 optical branches. At Γ the
optical modes move the two layers rigidly against each other, in the plane (twice) or
along c; only the 6 springs between the layers act, so with d the interlayer neighbour
distance omega^2 = 2 K a^2 / (M d^2) in the plane and 3 K c^2 / (M d^2) along c.

Mg has c/a = 1.624, just below the ideal sqrt(8/3) = 1.633 of touching spheres, so the 6
neighbours in the next layer (3.198 Å) are a little closer than the 6 in the same layer
(3.209 Å) and form their own shell; both get the same spring constant here. Along Γ-A
(q along c) the branches pair up at A, where the cell is doubled by the ABAB stacking.
Spring constants are illustrative.
"""

import matplotlib.pyplot as plt

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.plot import plot_bands, plot_crystal
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import hcp
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A = 3.209  # Å
C = 5.211  # Å
MASS = 24.305  # amu
SPRING_CONSTANT = 0.6  # eV/Å², the 12 nearest neighbours (two shells, same K)
LABELS = ["Γ", "M", "K", "Γ", "A", "L", "H", "A"]


def main():
    magnesium = hcp(A, C, [MASS, MASS], ["Mg (layer A)", "Mg (layer B)"])
    bonds = neighbour_bonds(magnesium, [SPRING_CONSTANT] * 2)
    points = special_points("hexagonal", magnesium.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])
    frequencies = dispersion(magnesium, bonds, qs) * OMEGA_TO_THZ

    ax = plot_bands(distance, frequencies, corners, LABELS)
    ax.set_ylabel(r"$\nu$ (THz)")
    ax.set_title(f"hcp Mg, nearest neighbours, K = {SPRING_CONSTANT} eV/Å²")
    plt.tight_layout()

    ax = plot_crystal(magnesium, radii=[0.4, 0.4], colors=["tab:olive", "tab:pink"], cells=2)
    ax.set_title("hcp Mg: layers A and B stacked ABAB (atoms not to scale); primitive cell")
    plt.show()


if __name__ == "__main__":
    main()
