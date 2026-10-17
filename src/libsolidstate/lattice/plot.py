"""Plot 2D crystal structures to scale with matplotlib."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

from libsolidstate.lattice.bravais2d import atom_positions


def plot_crystal(crystal, radii, colors, cells=3, ax=None):
    """Draw atoms as circles with radius in Å and outline the primitive cell."""
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
