"""In-plane phonons of a zigzag graphene nanoribbon, a strip periodic along x with free edges.

The ribbon is a crystal with one lattice vector in 2D (periodicity p = 1, dimension
d = 2): periodic along x, free across. Its bands are plotted against q_x and coloured by the
weight of each mode on the edge atoms. The lowest branch rises as q_x^2 (in-plane
bending of the strip). Some modes carry most of their weight on the edges: with
8 rows, up to about 70 % on the outer 25 % of the atoms.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_crystal
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

    qs = np.linspace(0, np.pi / A, 200)
    q_vectors = np.column_stack([qs, np.zeros_like(qs)])
    omega, vectors = phonon_modes(ribbon, bonds, q_vectors)
    frequencies = omega * OMEGA_TO_THZ
    weights = atom_weights(vectors, ribbon.dim)[:, :, edge].sum(axis=(2, 3))

    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5), width_ratios=[2, 1])
    q_all = np.repeat(qs * A, frequencies.shape[1])
    points = left.scatter(q_all, frequencies.ravel(), c=weights.ravel(), s=2, cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(points, ax=left, label="weight on edge atoms")
    left.axhline(0, color="grey", lw=0.5)
    left.set_xlabel(r"$q_x a$")
    left.set_ylabel(r"$\nu$ (THz)")
    left.set_xticks([0, np.pi / 2, np.pi], ["0", r"$\pi/2$", r"$\pi$"])
    left.set_title(f"zigzag ribbon, {ROWS} rows; edge atoms are {edge.mean():.0%} of all")

    plot_crystal(ribbon, radii=[0.35] * 2 * ROWS, colors=["tab:orange" if e else "tab:blue" for e in edge], ax=right)
    right.legend().remove()
    right.set_title("ribbon (orange: edge atoms), one cell outlined")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
