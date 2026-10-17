"""The five 2D Bravais lattices and crystal structures built on them."""

from dataclasses import dataclass

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


@dataclass
class Crystal:
    """Bravais lattice plus basis.

    vectors:   (2, 2) primitive vectors as rows, Å
    positions: (B, 2) Cartesian positions of the basis atoms, Å
    masses:    (B,) amu
    names:     B labels
    """

    vectors: np.ndarray
    positions: np.ndarray
    masses: np.ndarray
    names: list

    def __post_init__(self):
        self.vectors = np.asarray(self.vectors, dtype=float)
        self.positions = np.asarray(self.positions, dtype=float)
        self.masses = np.asarray(self.masses, dtype=float)


def atom_positions(crystal, cells):
    """Positions of all atoms in the cells n1 a1 + n2 a2 with |n1|, |n2| <= cells.

    Returns positions (M, 2) and the basis index (M,) of each atom.
    """
    n = np.arange(-cells, cells + 1)
    n1, n2 = np.meshgrid(n, n)
    translations = np.column_stack([n1.ravel(), n2.ravel()]) @ crystal.vectors

    positions = (translations[:, None, :] + crystal.positions[None, :, :]).reshape(-1, 2)
    basis_index = np.tile(np.arange(len(crystal.positions)), len(translations))
    return positions, basis_index
