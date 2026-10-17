"""In-plane phonon band structures of 2D crystals along high-symmetry paths.

Central springs to nearest neighbours only leave the square lattice without shear
stiffness (transverse branch at omega = 0 along Gamma-X) and graphene with a whole
zero-frequency branch. Second-neighbour springs stabilize both.
Spring constants are illustrative, not fitted to experiment; out-of-plane
modes need 3D displacements and are not included.
"""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.lattice.bravais2d import Crystal, primitive_vectors
from libsolidstate.lattice.reciprocal import k_path, special_points
from libsolidstate.lattice.springs import neighbour_bonds
from libsolidstate.spectra.dispersion import crystal_dispersion
from libsolidstate.spectra.plot import plot_bands
from libsolidstate.units import OMEGA_TO_THZ

A_SQUARE = 2.5  # Å
A_HEXAGONAL = 2.46  # Å
K_SIMPLE = [10.0, 3.0]  # eV/Å², first and second shell
K_GRAPHENE = [25.0, 5.0]  # eV/Å²


def honeycomb(masses, names):
    a = A_HEXAGONAL
    return Crystal(primitive_vectors("hexagonal", a), [[0.0, 0.0], [a / 2, a / (2 * np.sqrt(3))]], masses, names)


def main():
    square = Crystal(primitive_vectors("square", A_SQUARE), [[0.0, 0.0]], [12.011], ["C"])
    triangular = Crystal(primitive_vectors("hexagonal", A_HEXAGONAL), [[0.0, 0.0]], [12.011], ["C"])
    graphene = honeycomb([12.011, 12.011], ["C", "C"])
    hbn = honeycomb([10.81, 14.007], ["B", "N"])

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
        frequencies = crystal_dispersion(crystal, bonds, qs) * OMEGA_TO_THZ

        plot_bands(distance, frequencies, corners, labels, ax=axis)
        axis.set_title(f"{title}\nk = {spring_constants} eV/Å²")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
