"""2D crystal structures drawn to scale, with the primitive cell outlined."""

import matplotlib.pyplot as plt

from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_crystal
from latticephys.lattice.structures import honeycomb, kagome, lieb

A_NACL = 5.64  # Å, conventional cell
A_GRAPHENE = 2.46  # Å, lattice constant (C-C distance a / sqrt(3) = 1.42 Å)
A_HBN = 2.50  # Å
A_TRIANGULAR = 3.0  # Å
A_KAGOME = 4.0  # Å
A_CUO2 = 3.85  # Å
A_AFM = 3.0  # Å


def main():
    structures = [
        (
            "NaCl (001) plane",
            Crystal(primitive_vectors("centered rectangular", A_NACL, A_NACL),
                    [[0.0, 0.0], [A_NACL / 2, 0.0]], [22.990, 35.45], ["Na+", "Cl-"]),
            [1.02, 1.81], ["red", "blue"],
        ),
        (
            "Graphene",
            honeycomb(A_GRAPHENE, [12.011, 12.011], ["C (A)", "C (B)"]),
            [0.4, 0.4], ["blue", "red"],
        ),
        (
            "Hexagonal boron nitride",
            honeycomb(A_HBN, [10.81, 14.007], ["B", "N"]),
            [0.4, 0.4], ["orange", "blue"],
        ),
        (
            "Triangular lattice",
            Crystal(primitive_vectors("hexagonal", A_TRIANGULAR), [[0.0, 0.0]], [1.0], ["A"]),
            [0.6], ["blue"],
        ),
        (
            "Kagome lattice",
            kagome(A_KAGOME, [1.0, 1.0, 1.0], ["A", "B", "C"]),
            [0.5, 0.5, 0.5], ["blue", "red", "green"],
        ),
        (
            "CuO2 plane (Lieb lattice)",
            lieb(A_CUO2, [63.546, 15.999, 15.999], ["Cu", "O (x)", "O (y)"]),
            [0.7, 0.5, 0.5], ["tab:brown", "red", "tab:pink"],
        ),
        (
            "Square antiferromagnet",
            Crystal(primitive_vectors("square", A_AFM),
                    [[0.0, 0.0], [A_AFM / 2, A_AFM / 2]], [54.938, 54.938], ["Mn up", "Mn down"]),
            [0.5, 0.5], ["red", "blue"],
        ),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    for axis, (title, crystal, radii, colors) in zip(axes.flat, structures):
        plot_crystal(crystal, radii, colors, cells=4, ax=axis)
        axis.set_title(title)
    axes.flat[-1].axis("off")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
