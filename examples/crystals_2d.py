"""2D crystal structures: rock-salt (001) plane, graphene and a square antiferromagnet."""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.lattice.bravais2d import Crystal, primitive_vectors
from libsolidstate.lattice.plot import plot_crystal

A_NACL = 5.64  # Å, conventional cell
A_GRAPHENE = 2.46  # Å, lattice constant (C-C distance a / sqrt(3) = 1.42 Å)
A_AFM = 3.0  # Å


def main():
    nacl = Crystal(
        primitive_vectors("centered rectangular", A_NACL, A_NACL),
        positions=[[0.0, 0.0], [A_NACL / 2, 0.0]],
        masses=[22.990, 35.45],
        names=["Na+", "Cl-"],
    )
    graphene = Crystal(
        primitive_vectors("hexagonal", A_GRAPHENE),
        positions=[[0.0, 0.0], [A_GRAPHENE / 2, A_GRAPHENE / (2 * np.sqrt(3))]],
        masses=[12.011, 12.011],
        names=["C (A)", "C (B)"],
    )
    afm = Crystal(
        primitive_vectors("square", A_AFM),
        positions=[[0.0, 0.0], [A_AFM / 2, A_AFM / 2]],
        masses=[54.938, 54.938],
        names=["Mn up", "Mn down"],
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    plot_crystal(nacl, radii=[1.02, 1.81], colors=["red", "blue"], ax=axes[0])
    plot_crystal(graphene, radii=[0.4, 0.4], colors=["blue", "red"], cells=4, ax=axes[1])
    plot_crystal(afm, radii=[0.5, 0.5], colors=["red", "blue"], ax=axes[2])
    axes[0].set_title("NaCl (001) plane")
    axes[1].set_title("Graphene")
    axes[2].set_title("Square antiferromagnet")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
