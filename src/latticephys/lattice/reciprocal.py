"""Reciprocal lattice (any dimension); 2D Brillouin zone; special points and k-paths of 2D lattices and layers."""

import numpy as np

from latticephys.lattice.crystal import cell_offsets


def reciprocal_vectors(vectors):
    """Reciprocal vectors b_j as rows, with a_i . b_j = 2 pi delta_ij (1/Å).

    vectors: (p, d) lattice vectors as rows, p <= d. b = 2 pi (A A^T)^-1 A lies in
    the span of the a_i; for p = d this is 2 pi A^-T.
    """
    vectors = np.atleast_2d(np.asarray(vectors, dtype=float))
    return 2 * np.pi * np.linalg.solve(vectors @ vectors.T, vectors)


def lattice_projection(vectors, q):
    """Part of the wave vector q (1/Å, d components) in the span of the p lattice vectors.

    q = sum_j (a_j . q / 2 pi) b_j; equal to q for p = d, the in-plane part for a layer.
    """
    vectors = np.atleast_2d(np.asarray(vectors, dtype=float))
    return vectors @ np.asarray(q, dtype=float) / (2 * np.pi) @ reciprocal_vectors(vectors)


def q_grid(vectors, n):
    """The n^p wave vectors q = (i_1 b_1 + ... + i_p b_p) / n, 0 <= i_j < n (1/Å), for p lattice vectors.

    These are exactly the wave vectors of a periodic supercell of n^p cells,
    so sums over them are Brillouin-zone averages for a crystal of that size.
    """
    vectors = np.atleast_2d(vectors)
    return cell_offsets(n, len(vectors)) / n @ reciprocal_vectors(vectors)


def wigner_seitz_cell(vectors):
    """Vertices (counter-clockwise) of the Wigner-Seitz cell of a 2D lattice.

    Clips a large square by the bisector of every lattice vector up to two
    steps along a1 and a2, which is enough for the short primitive vectors
    returned by primitive_vectors().
    """
    vectors = np.asarray(vectors, dtype=float)
    if vectors.shape != (2, 2):
        raise ValueError("the Wigner-Seitz cell and Brillouin zone need 2 lattice vectors in 2D (p = d = 2)")
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


def special_points(lattice, vectors):
    """High-symmetry points of the Brillouin zone (1/Å) for the lattices with standard labels.

    "square":      Gamma, X (edge centre), M (corner)
    "rectangular": Gamma, X, Y (edge centres along b1, b2), S (corner)
    "hexagonal":   Gamma, M (edge centre), K (corner)
    Assumes the vectors of primitive_vectors() (for hexagonal: 60° between a1, a2).
    vectors: (2, d), e.g. a layer in 3D; the points lie in its plane.
    """
    if np.ndim(vectors) != 2 or len(vectors) != 2:
        raise ValueError("special_points needs 2 lattice vectors")
    b1, b2 = reciprocal_vectors(vectors)
    gamma = np.zeros_like(b1)
    if lattice == "square":
        return {"Γ": gamma, "X": b1 / 2, "M": (b1 + b2) / 2}
    if lattice == "rectangular":
        return {"Γ": gamma, "X": b1 / 2, "Y": b2 / 2, "S": (b1 + b2) / 2}
    if lattice == "hexagonal":
        return {"Γ": gamma, "M": b1 / 2, "K": (2 * b1 + b2) / 3}
    raise ValueError(f"no standard labels for {lattice} lattice, give the k-points explicitly")


def k_path(points, n=300):
    """About n wave vectors on straight segments through the given k-points.

    Returns qs (M, d), the distance along the path (M,) and the distance of
    each given point (len(points),), e.g. for tick positions.
    """
    points = np.asarray(points, dtype=float)
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    corners = np.concatenate([[0.0], np.cumsum(lengths)])

    qs, distance = [], []
    for start, end, d0, length in zip(points[:-1], points[1:], corners[:-1], lengths):
        t = np.linspace(0, 1, max(2, round(n * length / corners[-1])), endpoint=False)
        qs.append(start + t[:, None] * (end - start))
        distance.append(d0 + t * length)
    qs.append(points[-1:])
    distance.append(corners[-1:])
    return np.concatenate(qs), np.concatenate(distance), corners
