"""Plot crystal structures (2D and 3D), 2D Brillouin zones and band structures with matplotlib."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Polygon

from latticephys.lattice.crystal import atom_positions
from latticephys.lattice.reciprocal import brillouin_zone, reciprocal_vectors


def plot_crystal(crystal, radii, colors, cells=3, ax=None):
    """Draw a crystal of dimension 2 or 3: atoms coloured per basis atom, primitive cell outlined.

    radii, colors: one per basis atom
    cells:         the view reaches cells / 2 lattice extents from the origin

    2D: atoms are circles with radius in Å. For periodicity 1 (e.g. a ribbon) the outline
    runs along a1 from the origin to a1 and across the whole basis, plus one atom radius on
    each side. A molecule (p = 0) has no outline.
    3D: ax must be a 3D axis (projection="3d") or None. Atoms in an axis-aligned box are
    scatter markers about 2 radii across (approximate, for the default 7 inch figure); the
    cell is drawn by its edges (12 for p = 3, 4 for a layer, 1 for a chain).
    """
    if len(radii) != len(crystal.positions) or len(colors) != len(crystal.positions):
        raise ValueError(f"need one radius and one colour per basis atom ({len(crystal.positions)})")
    if cells < 1:
        raise ValueError("cells must be at least 1")
    if crystal.dim == 3:
        if ax is not None and ax.name != "3d":
            raise ValueError("a 3D crystal needs a 3D axis, e.g. add_subplot(projection=\"3d\")")
        return _plot_crystal_3d(crystal, radii, colors, cells, ax)
    if crystal.dim != 2:
        raise ValueError("plot_crystal needs dimension 2 or 3")
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    positions, basis_index = atom_positions(crystal, cells)
    for (x, y), mu in zip(positions, basis_index):
        ax.add_patch(Circle((x, y), radii[mu], color=colors[mu], alpha=0.6))

    margin = max(radii)
    if crystal.periodicity == 2:
        a1, a2 = crystal.vectors
        outline = np.array([(0, 0), a1, a1 + a2, a2])
        # the central region, where the cells cover the plot completely
        extent = cells / 2 * np.max(np.abs(crystal.vectors).sum(axis=0))
        view = extent * np.array([[-1, -1], [1, 1]])
    elif crystal.periodicity == 1:
        a1 = crystal.vectors[0]
        normal = np.array([-a1[1], a1[0]]) / np.linalg.norm(a1)
        across = crystal.positions @ normal
        low, high = across.min() - margin, across.max() + margin
        outline = np.array([low * normal, low * normal + a1, high * normal + a1, high * normal])
        # the strip along a1 within cells / 2 cells of the origin, the outline just inside
        view = np.array([along * cells / 2 * a1 + offset * normal
                         for along in (-1, 1) for offset in (low - margin, high + margin)])
    else:
        outline = np.empty((0, 2))  # a molecule has no cell
        view = np.array([crystal.positions.min(axis=0) - margin, crystal.positions.max(axis=0) + margin])
    if len(outline):
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


def _plot_crystal_3d(crystal, radii, colors, cells, ax):
    """plot_crystal for d = 3: the box reaches cells / 2 lattice extents from the origin along the
    axes the lattice vectors span and covers the basis plus one atom radius along the others."""
    if ax is None:
        ax = plt.figure(figsize=(7, 7)).add_subplot(projection="3d")
    margin = max(radii)
    low, high = crystal.positions.min(axis=0) - margin, crystal.positions.max(axis=0) + margin
    if crystal.periodicity:
        spanned = np.abs(crystal.vectors).sum(axis=0) > 1e-12
        extent = cells / 2 * np.abs(crystal.vectors).sum(axis=0).max()
        low[spanned], high[spanned] = -extent, extent

    # enough cells to fill the box: the largest lattice coordinate of a box corner
    corners = np.array([np.where(bits, high, low) for bits in np.ndindex(2, 2, 2)])
    reach = np.abs(corners @ reciprocal_vectors(crystal.vectors).T / (2 * np.pi)).max(initial=0)
    positions, basis_index = atom_positions(crystal, int(np.ceil(reach)) + 1)
    inside = np.all((positions >= low - 1e-9) & (positions <= high + 1e-9), axis=1)

    points_per_angstrom = 150 / (high - low).max()  # about the box width in points on a 7 inch figure
    for mu, (name, color) in enumerate(zip(crystal.names, colors)):
        x, y, z = positions[inside & (basis_index == mu)].T
        ax.scatter(x, y, z, s=(2 * radii[mu] * points_per_angstrom) ** 2, color=color, alpha=0.8, label=name)

    vertices = {bits: np.array(bits) @ crystal.vectors for bits in np.ndindex((2,) * crystal.periodicity)}
    for bits, start in vertices.items():
        for i in np.flatnonzero(np.array(bits) == 0):
            end = vertices[bits[:i] + (1,) + bits[i + 1:]]
            ax.plot(*np.transpose([start, end]), color="orange", ls="-.")

    legend = ax.legend(loc="upper right")
    for handle in legend.legend_handles:
        handle.set_sizes([40])
    ax.set_xlim(low[0], high[0])
    ax.set_ylim(low[1], high[1])
    ax.set_zlim(low[2], high[2])
    ax.set_box_aspect(high - low)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    return ax


def plot_brillouin_zone(vectors, ax=None):
    """Draw the first Brillouin zone, reciprocal lattice points and b1, b2 (p = d = 2 only)."""
    zone = brillouin_zone(vectors)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    b1, b2 = reciprocal_vectors(vectors)
    n = np.arange(-2, 3)
    points = np.array([i * b1 + j * b2 for i in n for j in n])
    ax.scatter(points[:, 0], points[:, 1], color="black", s=10)

    ax.add_patch(Polygon(zone, alpha=0.3, color="tab:blue"))
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


def plot_bands(distance, frequencies, corners, labels, ax=None, weights=None, weight_label=None):
    """Branches against the distance along a k-path, with the k-points as ticks.

    weights: optional, shape of frequencies, values in [0, 1] (e.g. from atom_weights()); the
             branches are then dots coloured by weight, with a colour bar labelled weight_label
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    if weights is None:
        ax.plot(distance, frequencies, color="tab:blue")
    else:
        frequencies, weights = np.asarray(frequencies), np.asarray(weights)
        if weights.shape != frequencies.shape:
            raise ValueError(f"weights must have the shape of frequencies {frequencies.shape}, got {weights.shape}")
        if weights.min() < -1e-9 or weights.max() > 1 + 1e-9:
            raise ValueError("weights must lie in [0, 1], the range of the colour bar")
        dots = ax.scatter(np.repeat(distance, frequencies.shape[1]), frequencies.ravel(), c=weights.ravel(),
                          s=4, cmap="viridis", vmin=0, vmax=1)
        ax.figure.colorbar(dots, ax=ax, label=weight_label)
    for corner in corners:
        ax.axvline(corner, color="grey", lw=0.5)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xticks(corners, labels)
    ax.set_xlim(corners[0], corners[-1])
    return ax


def plot_points(ax, corners, labels, values, label="measured"):
    """Mark values (e.g. measured frequencies) at the k-path corners, as drawn by plot_bands.

    values: {point label: list of values}; drawn at every corner with that label (Γ may appear
    twice), as red dots with one legend entry.
    """
    for point, ys in values.items():
        for corner in np.asarray(corners)[[i for i, name in enumerate(labels) if name == point]]:
            ax.plot([corner] * len(ys), ys, "o", color="tab:red", clip_on=False)
    ax.plot([], [], "o", color="tab:red", label=label)
    return ax
