"""In-plane phonon band structures of 2D crystals along high-symmetry paths.

Central springs to nearest neighbours only leave the square lattice without shear
stiffness (transverse branch at omega = 0 along Gamma-X) and graphene with a whole
zero-frequency branch. Second-neighbour springs stabilize both.
Spring constants are illustrative, not fitted to experiment; out-of-plane
modes need 3D displacements, see layers/graphene_layer.py.
"""

import matplotlib.pyplot as plt

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_bands
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A_SQUARE = 2.5  # Å
A_HEXAGONAL = 2.46  # Å
K_SIMPLE = [10.0, 3.0]  # eV/Å², first and second shell
K_GRAPHENE = [25.0, 5.0]  # eV/Å²


def main():
    square = Crystal(primitive_vectors("square", A_SQUARE), [[0.0, 0.0]], [12.011], ["C"])
    triangular = Crystal(primitive_vectors("hexagonal", A_HEXAGONAL), [[0.0, 0.0]], [12.011], ["C"])
    graphene = honeycomb(A_HEXAGONAL, [12.011, 12.011], ["C", "C"])
    hbn = honeycomb(A_HEXAGONAL, [10.81, 14.007], ["B", "N"])

    cases = [
        ("Triangular, 1st neighbours", triangular, "hexagonal", K_SIMPLE[:1]),
        ("Square, 1st neighbours", square, "square", K_SIMPLE[:1]),
        ("Square, 1st + 2nd neighbours", square, "square", K_SIMPLE),
        ("Graphene, 1st neighbours", graphene, "hexagonal", K_GRAPHENE[:1]),
        ("Graphene, 1st + 2nd neighbours", graphene, "hexagonal", K_GRAPHENE),
        ("hBN, 1st + 2nd neighbours", hbn, "hexagonal", K_GRAPHENE),
    ]
    paths = {"square": ["Γ", "X", "M", "Γ"], "hexagonal": ["Γ", "M", "K", "Γ"]}

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharey=True)
    for axis, (title, crystal, lattice, spring_constants) in zip(axes.flat, cases):
        points = special_points(lattice, crystal.vectors)
        labels = paths[lattice]
        qs, distance, corners = k_path([points[label] for label in labels])

        bonds = neighbour_bonds(crystal, spring_constants)
        frequencies = dispersion(crystal, bonds, qs) * OMEGA_TO_THZ

        plot_bands(distance, frequencies, corners, labels, ax=axis)
        axis.set_title(f"{title}\nK = {spring_constants} eV/Å²")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
