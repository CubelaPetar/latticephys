"""Plot 2D crystal structures to scale with matplotlib."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

from libsolidstate.lattice.crystal import atom_positions
from libsolidstate.lattice.reciprocal import brillouin_zone, reciprocal_vectors


def plot_crystal(crystal, radii, colors, cells=3, ax=None):
    """Draw a 2D crystal: atoms as circles with radius in Å, primitive cell outlined."""
    if crystal.dim != 2:
        raise ValueError("plot_crystal needs a 2D crystal")
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    positions, basis_index = atom_positions(crystal, cells)
    for (x, y), mu in zip(positions, basis_index):
        ax.add_patch(Circle((x, y), radii[mu], color=colors[mu], alpha=0.6))

    a1, a2 = crystal.vectors
    cell = Polygon([(0, 0), a1, a1 + a2, a2], fill=False, color="orange", ls="-.")
    ax.add_patch(cell)

    for name, color in zip(crystal.names, colors):
        ax.scatter([], [], color=color, alpha=0.6, label=name)
    ax.legend(loc="upper right")

    # show the central region, where the cells cover the plot completely
    extent = cells / 2 * np.max(np.abs(crystal.vectors).sum(axis=0))
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    return ax


def plot_brillouin_zone(vectors, ax=None):
    """Draw the first Brillouin zone, reciprocal lattice points and b1, b2."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    b1, b2 = reciprocal_vectors(vectors)
    n = np.arange(-2, 3)
    points = np.array([i * b1 + j * b2 for i in n for j in n])
    ax.scatter(points[:, 0], points[:, 1], color="black", s=10)

    ax.add_patch(Polygon(brillouin_zone(vectors), alpha=0.3, color="tab:blue"))
    for b, label in ((b1, "$b_1$"), (b2, "$b_2$")):
        ax.annotate("", xy=b, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="tab:red"))
        ax.annotate(label, xy=b, color="tab:red")
    ax.annotate(r"$\Gamma$", xy=(0, 0), xytext=(3, 3), textcoords="offset points")

    extent = 1.2 * np.max(np.linalg.norm([b1, b2, b1 + b2, b1 - b2], axis=1))
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$k_x$ (1/Å)")
    ax.set_ylabel(r"$k_y$ (1/Å)")
    return ax
