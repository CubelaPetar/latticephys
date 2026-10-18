"""Primitive vectors of the five 2D Bravais lattices."""

import numpy as np


def primitive_vectors(lattice, a, b=None, angle=None):
    """Primitive lattice vectors a1, a2 as rows of a (2, 2) array in Å.

    "oblique":              a, b and angle (degrees) between a1 and a2
    "rectangular":          a, b
    "centered rectangular": a, b of the conventional cell
    "square":               a
    "hexagonal":            a (lattice constant, |a1| = |a2| = a, 60° between them)
    """
    if lattice in ("oblique", "rectangular", "centered rectangular") and b is None:
        raise ValueError(f"{lattice} lattice needs b")

    if lattice == "oblique":
        if angle is None:
            raise ValueError("oblique lattice needs angle")
        phi = np.radians(angle)
        return np.array([[a, 0.0], [b * np.cos(phi), b * np.sin(phi)]])
    if lattice == "rectangular":
        return np.array([[a, 0.0], [0.0, b]])
    if lattice == "centered rectangular":
        return np.array([[a / 2, b / 2], [a / 2, -b / 2]])
    if lattice == "square":
        return np.array([[a, 0.0], [0.0, a]])
    if lattice == "hexagonal":
        return np.array([[a, 0.0], [a / 2, a * np.sqrt(3) / 2]])
    raise ValueError(f"unknown lattice: {lattice}")

