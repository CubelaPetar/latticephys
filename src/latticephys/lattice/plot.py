"""Plot 2D crystal structures, Brillouin zones and band structures with matplotlib."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

from latticephys.lattice.crystal import atom_positions
from latticephys.lattice.reciprocal import brillouin_zone, reciprocal_vectors


def plot_crystal(crystal, radii, colors, cells=3, ax=None):
    """Draw a crystal of dimension 2: atoms as circles with radius in Å, primitive cell outlined.

    For periodicity 1 (e.g. a ribbon) the outline runs along a1 from the origin to a1 and
    across the whole basis, plus one atom radius on each side.
    """
    if crystal.dim != 2 or crystal.periodicity == 0:
        raise ValueError("plot_crystal needs dimension 2 and periodicity 1 or 2")
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    positions, basis_index = atom_positions(crystal, cells)
    for (x, y), mu in zip(positions, basis_index):
        ax.add_patch(Circle((x, y), radii[mu], color=colors[mu], alpha=0.6))

    if crystal.periodicity == 2:
        a1, a2 = crystal.vectors
        outline = np.array([(0, 0), a1, a1 + a2, a2])
        # the central region, where the cells cover the plot completely
        extent = cells / 2 * np.max(np.abs(crystal.vectors).sum(axis=0))
        view = extent * np.array([[-1, -1], [1, 1]])
    else:
        a1 = crystal.vectors[0]
        normal = np.array([-a1[1], a1[0]]) / np.linalg.norm(a1)
        margin = max(radii)
        across = crystal.positions @ normal
        low, high = across.min() - margin, across.max() + margin
        outline = np.array([low * normal, low * normal + a1, high * normal + a1, high * normal])
        # the strip along a1 within cells / 2 cells of the origin, the outline just inside
        view = np.array([along * cells / 2 * a1 + offset * normal
                         for along in (-1, 1) for offset in (low - margin, high + margin)])
    ax.add_patch(Polygon(outline, fill=False, color="orange", ls="-."))

    for name, color in zip(crystal.names, colors):
        ax.scatter([], [], color=color, alpha=0.6, label=name)
    ax.legend(loc="upper right")

    ax.set_xlim(view[:, 0].min(), view[:, 0].max())
    ax.set_ylim(view[:, 1].min(), view[:, 1].max())
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


def plot_bands(distance, frequencies, corners, labels, ax=None):
    """Branches against the distance along a k-path, with the k-points as ticks."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    ax.plot(distance, frequencies, color="tab:blue")
    for corner in corners:
        ax.axvline(corner, color="grey", lw=0.5)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xticks(corners, labels)
    ax.set_xlim(corners[0], corners[-1])
    return ax
