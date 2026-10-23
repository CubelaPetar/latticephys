"""Crystal: Bravais lattice plus basis, in any dimension."""

from dataclasses import dataclass
from itertools import product

import numpy as np


@dataclass
class Crystal:
    """Bravais lattice plus basis in d dimensions, periodic along p <= d lattice vectors.

    vectors:   (p, d) primitive vectors as rows, linearly independent, Å; p < d for
               e.g. a layer (p = 2, d = 3) or a chain moving in a plane (p = 1, d = 2)
    positions: (B, d) Cartesian positions of the basis atoms, Å
    masses:    (B,) amu
    names:     B labels

    A molecule has no lattice vectors: vectors np.zeros((0, d)), bonds with cell index ().
    """

    vectors: np.ndarray
    positions: np.ndarray
    masses: np.ndarray
    names: list

    def __post_init__(self):
        self.vectors = np.atleast_2d(np.asarray(self.vectors, dtype=float))
        self.masses = np.atleast_1d(np.asarray(self.masses, dtype=float))
        self.positions = np.asarray(self.positions, dtype=float).reshape(len(self.masses), -1)

        p, d = len(self.vectors), self.positions.shape[1]
        if self.vectors.shape != (p, d) or p > d:
            raise ValueError(f"vectors must be (p, {d}) with p <= {d} for positions in {d}D, got {self.vectors.shape}")
        if np.linalg.matrix_rank(self.vectors) < p:
            raise ValueError("lattice vectors must be linearly independent")
        if len(self.names) != len(self.masses):
            raise ValueError("need one name per atom")

    @property
    def dim(self):
        """Dimension d: components of positions, displacements and wave vectors."""
        return self.positions.shape[1]

    @property
    def periodicity(self):
        """Periodicity p: number of lattice vectors."""
        return len(self.vectors)

    def cell_index(self, n):
        """Cell index n (p integers, or an int for p = 1) as an integer array; raises for the wrong length."""
        index = np.atleast_1d(np.asarray(n, dtype=int))
        if index.shape != (self.periodicity,):
            raise ValueError(f"cell index must have {self.periodicity} entries, got {n}")
        return index

    def bond_vector(self, mu, nu, n):
        """Vector from atom mu in cell 0 to atom nu in cell n (p integers, or an int for p = 1)."""
        return self.positions[nu] + self.cell_index(n) @ self.vectors - self.positions[mu]


def embed(crystal, d):
    """The crystal in d dimensions: zero components appended to vectors and positions.

    E.g. honeycomb() lifted into 3D is a graphene layer in the xy plane. Raises if d
    is smaller than the dimension of the crystal.
    """
    extra = d - crystal.dim
    if extra < 0:
        raise ValueError(f"cannot embed a {crystal.dim}D crystal in {d}D")
    return Crystal(
        np.pad(crystal.vectors, ((0, 0), (0, extra))),
        np.pad(crystal.positions, ((0, 0), (0, extra))),
        crystal.masses,
        crystal.names,
    )


def cell_offsets(cells, dim):
    """All integer cell tuples with entries 0 .. cells - 1, last index fastest."""
    return np.array(list(product(range(cells), repeat=dim)), dtype=int).reshape(-1, dim)


def atom_positions(crystal, cells):
    """Positions of all atoms in the cells n . vectors with all |n_i| <= cells.

    Returns positions (M, d) and the basis index (M,) of each atom.
    """
    offsets = cell_offsets(2 * cells + 1, crystal.periodicity) - cells
    translations = offsets @ crystal.vectors

    positions = (translations[:, None, :] + crystal.positions[None, :, :]).reshape(-1, crystal.dim)
    basis_index = np.tile(np.arange(len(crystal.positions)), len(translations))
    return positions, basis_index
