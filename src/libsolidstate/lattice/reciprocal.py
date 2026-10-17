"""Reciprocal lattice and Brillouin zone of 2D lattices."""

import numpy as np


def reciprocal_vectors(vectors):
    """Reciprocal vectors b1, b2 as rows, with a_i . b_j = 2 pi delta_ij (1/Å)."""
    return 2 * np.pi * np.linalg.inv(vectors).T


def wigner_seitz_cell(vectors):
    """Vertices (counter-clockwise) of the Wigner-Seitz cell of a 2D lattice.

    Clips a large square by the bisector of every lattice vector up to two
    steps along a1 and a2, which is enough for the short primitive vectors
    returned by primitive_vectors().
    """
    vectors = np.asarray(vectors, dtype=float)
    n = range(-2, 3)
    points = [i * vectors[0] + j * vectors[1] for i in n for j in n if (i, j) != (0, 0)]

    size = 2 * np.abs(vectors).sum()
    polygon = size * np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    for g in points:
        polygon = _clip(polygon, g, g @ g / 2)
    return polygon


def brillouin_zone(vectors):
    """Vertices (counter-clockwise) of the first Brillouin zone (1/Å)."""
    return wigner_seitz_cell(reciprocal_vectors(vectors))


def _clip(polygon, normal, offset):
    """Keep the part of a convex polygon with x . normal <= offset (Sutherland-Hodgman)."""
    tolerance = 1e-9 * abs(offset)
    result = []
    for p, q in zip(polygon, np.roll(polygon, -1, axis=0)):
        p_inside = p @ normal <= offset + tolerance
        q_inside = q @ normal <= offset + tolerance
        if p_inside:
            result.append(p)
        if p_inside != q_inside:
            t = (offset - p @ normal) / ((q - p) @ normal)
            result.append(p + t * (q - p))

    # drop points that coincide when the bisector passes through a vertex
    result = np.array(result)
    distinct = np.linalg.norm(result - np.roll(result, -1, axis=0), axis=1) > 1e-9 * np.abs(result).max()
    return result[distinct]
