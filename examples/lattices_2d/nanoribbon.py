"""In-plane phonons of a zigzag graphene nanoribbon, a strip periodic along x with free edges.

The ribbon is a crystal with one lattice vector in 2D (periodicity p = 1, dimension
d = 2): periodic along x, free across. Its bands are plotted along Γ-X (q along x, X = pi/a)
and coloured by the weight of each mode on the edge atoms. The lowest branch rises as q_x^2 (in-plane
bending of the strip). Some modes carry most of their weight on the edges: with
8 rows, up to about 70 % on the outer 25 % of the atoms.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_bands, plot_crystal
from latticephys.lattice.reciprocal import k_path
from latticephys.phonons.modes import atom_weights, phonon_modes
from latticephys.units import OMEGA_TO_THZ

A = 2.46  # Å, graphene lattice constant = period of the zigzag ribbon
MASS = 12.011  # amu
SPRING_CONSTANTS = [25.0, 5.0]  # eV/Å², first and second neighbours
ROWS = 8  # zigzag chains across the ribbon
EDGE_ROWS = 1  # zigzag chains counted as edge on each side


def zigzag_ribbon():
    """Two atoms per zigzag chain; chain j is shifted by j a2 of graphene, folded into the cell."""
    height = A * np.sqrt(3) / 2
    positions = []
    for j in range(ROWS):
        x = (j * A / 2) % A
        positions += [[x, j * height], [(x + A / 2) % A, j * height + A / (2 * np.sqrt(3))]]
    return Crystal([[A, 0.0]], positions, [MASS] * 2 * ROWS, ["C"] * 2 * ROWS)


def main():
    ribbon = zigzag_ribbon()
    bonds = neighbour_bonds(ribbon, SPRING_CONSTANTS)
    edge = np.zeros(2 * ROWS, dtype=bool)
    edge[: 2 * EDGE_ROWS] = edge[-2 * EDGE_ROWS :] = True

    qs, distance, corners = k_path([[0.0, 0.0], [np.pi / A, 0.0]], n=200)
    omega, vectors = phonon_modes(ribbon, bonds, qs)
    frequencies = omega * OMEGA_TO_THZ
    weights = atom_weights(vectors, ribbon.dim)[:, :, edge].sum(axis=(2, 3))

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5), width_ratios=[2, 1])
    plot_bands(distance, frequencies, corners, ["Γ", "X"], ax=left, weights=weights,
               weight_label="weight on edge atoms")
    left.set_ylabel(r"$\nu$ (THz)")
    left.set_title(f"zigzag ribbon, {ROWS} rows; edge atoms are {edge.mean():.0%} of all")

    plot_crystal(ribbon, radii=[0.35] * 2 * ROWS, colors=["tab:orange" if e else "tab:blue" for e in edge], ax=right)
    right.legend().remove()
    right.set_title("ribbon (orange: edge atoms), one cell outlined")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
