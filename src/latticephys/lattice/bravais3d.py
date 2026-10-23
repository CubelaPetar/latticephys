"""Primitive vectors of the 14 3D Bravais lattices.

Convention of W. Setyawan and S. Curtarolo, Comput. Mater. Sci. 49, 299 (2010), except
"hexagonal", which has 60° between a1 and a2 as in 2D. Angles in degrees: alpha between
b and c, beta between a and c, gamma between a and b.
"""

import numpy as np

# parameters each lattice takes besides a
PARAMETERS = {
    "cubic": (),
    "fcc": (),
    "bcc": (),
    "tetragonal": ("c",),
    "body-centered tetragonal": ("c",),
    "orthorhombic": ("b", "c"),
    "face-centered orthorhombic": ("b", "c"),
    "body-centered orthorhombic": ("b", "c"),
    "base-centered orthorhombic": ("b", "c"),
    "hexagonal": ("c",),
    "rhombohedral": ("alpha",),
    "monoclinic": ("b", "c", "alpha"),
    "base-centered monoclinic": ("b", "c", "alpha"),
    "triclinic": ("b", "c", "alpha", "beta", "gamma"),
}


def primitive_vectors(lattice, a, b=None, c=None, alpha=None, beta=None, gamma=None):
    """Primitive lattice vectors a1, a2, a3 as rows of a (3, 3) array in Å.

    "cubic", "fcc", "bcc":                        a
    "tetragonal", "body-centered tetragonal":     a, c
    "orthorhombic" and its "face-centered", "body-centered", "base-centered" (C) forms: a, b, c
    "hexagonal":                                  a, c
    "rhombohedral":                               a, alpha
    "monoclinic", "base-centered monoclinic" (C): a, b, c of the conventional cell, alpha
    "triclinic":                                  a, b, c, alpha, beta, gamma
    Raises for an unknown lattice, a missing or unused parameter, or angles that form no cell.
    """
    if lattice not in PARAMETERS:
        raise ValueError(f"unknown lattice: {lattice}")
    given = {"b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}
    missing = [name for name in PARAMETERS[lattice] if given[name] is None]
    unused = [name for name, value in given.items() if value is not None and name not in PARAMETERS[lattice]]
    if missing or unused:
        raise ValueError(f"{lattice} lattice takes a, {', '.join(PARAMETERS[lattice]) or 'nothing else'}; "
                         f"missing {missing}, unused {unused}")

    if lattice == "cubic":
        return a * np.eye(3)
    if lattice == "fcc":
        return a / 2 * np.array([[0.0, 1, 1], [1, 0, 1], [1, 1, 0]])
    if lattice == "bcc":
        return a / 2 * np.array([[-1.0, 1, 1], [1, -1, 1], [1, 1, -1]])
    if lattice == "tetragonal":
        return np.diag([a, a, c])
    if lattice == "body-centered tetragonal":
        return np.array([[-a, a, c], [a, -a, c], [a, a, -c]]) / 2
    if lattice == "orthorhombic":
        return np.diag([a, b, c])
    if lattice == "face-centered orthorhombic":
        return np.array([[0.0, b, c], [a, 0, c], [a, b, 0]]) / 2
    if lattice == "body-centered orthorhombic":
        return np.array([[-a, b, c], [a, -b, c], [a, b, -c]]) / 2
    if lattice == "base-centered orthorhombic":
        return np.array([[a / 2, -b / 2, 0.0], [a / 2, b / 2, 0], [0, 0, c]])
    if lattice == "hexagonal":
        return np.array([[a, 0.0, 0], [a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])

    ca, cb, cg = (np.cos(np.radians(angle)) if angle is not None else None for angle in (alpha, beta, gamma))
    if lattice == "rhombohedral":
        if not 0 < alpha < 120:
            raise ValueError("rhombohedral alpha must lie between 0 and 120 degrees")
        half = np.radians(alpha) / 2
        return np.array([
            [a * np.cos(half), -a * np.sin(half), 0.0],
            [a * np.cos(half), a * np.sin(half), 0.0],
            [a * ca / np.cos(half), 0.0, a * np.sqrt(1 - ca**2 / np.cos(half) ** 2)],
        ])
    sa = np.sin(np.radians(alpha))
    if lattice == "monoclinic":
        return np.array([[a, 0.0, 0], [0, b, 0], [0, c * ca, c * sa]])
    if lattice == "base-centered monoclinic":
        return np.array([[a / 2, b / 2, 0.0], [-a / 2, b / 2, 0], [0, c * ca, c * sa]])
    # triclinic
    height = 1 + 2 * ca * cb * cg - ca**2 - cb**2 - cg**2  # (volume / a b c)^2
    if height <= 0:
        raise ValueError(f"angles {alpha}, {beta}, {gamma} form no cell")
    sg = np.sin(np.radians(gamma))
    return np.array([
        [a, 0.0, 0.0],
        [b * cg, b * sg, 0.0],
        [c * cb, c * (ca - cb * cg) / sg, c * np.sqrt(height) / sg],
    ])
