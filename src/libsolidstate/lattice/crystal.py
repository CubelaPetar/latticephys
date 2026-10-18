"""Crystal: Bravais lattice plus basis, in any dimension."""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass
class Crystal:
    """Bravais lattice plus basis in d dimensions.

    vectors:   (d, d) primitive vectors as rows, Å
    positions: (B, d) Cartesian positions of the basis atoms, Å
    masses:    (B,) amu
    names:     B labels

    A finite system (molecule) is a crystal whose bonds all stay in cell 0;
    its vectors are then not used.
    """

    vectors: np.ndarray
    positions: np.ndarray
    masses: np.ndarray
    names: list

    def __post_init__(self):
        self.vectors = np.atleast_2d(np.asarray(self.vectors, dtype=float))
        self.masses = np.atleast_1d(np.asarray(self.masses, dtype=float))
        self.positions = np.asarray(self.positions, dtype=float).reshape(len(self.masses), -1)

        d = len(self.vectors)
        if self.vectors.shape != (d, d):
            raise ValueError(f"vectors must be (d, d), got {self.vectors.shape}")
        if self.positions.shape[1] != d:
            raise ValueError(f"positions must be (B, {d}), got {self.positions.shape}")
        if len(self.names) != len(self.masses):
            raise ValueError("need one name per atom")

    @property
    def dim(self):
        return len(self.vectors)

    def bond_vector(self, mu, nu, n):
        """Vector from atom mu in cell 0 to atom nu in cell n (integer tuple, or int in 1D)."""
        return self.positions[nu] + np.atleast_1d(n) @ self.vectors - self.positions[mu]


def cell_offsets(cells, dim):
    """All integer cell tuples with entries 0 .. cells - 1, last index fastest."""
    return np.array(list(product(range(cells), repeat=dim)), dtype=int).reshape(-1, dim)


def atom_positions(crystal, cells):
    """Positions of all atoms in the cells n . vectors with all |n_i| <= cells.

    Returns positions (M, d) and the basis index (M,) of each atom.
    """
    offsets = cell_offsets(2 * cells + 1, crystal.dim) - cells
    translations = offsets @ crystal.vectors

    positions = (translations[:, None, :] + crystal.positions[None, :, :]).reshape(-1, crystal.dim)
    basis_index = np.tile(np.arange(len(crystal.positions)), len(translations))
    return positions, basis_index
