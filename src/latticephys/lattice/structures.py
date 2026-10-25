"""Common crystal structures with several basis atoms; bilayers of any crystal."""

import numpy as np

from latticephys.lattice import bravais3d
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal


def honeycomb(a, masses, names):
    """Honeycomb lattice (graphene, hBN): hexagonal lattice constant a, two atoms a / sqrt(3) apart."""
    return Crystal(primitive_vectors("hexagonal", a), [[0.0, 0.0], [a / 2, a / (2 * np.sqrt(3))]], masses, names)


def kagome(a, masses, names):
    """Kagome lattice: hexagonal lattice constant a, atoms at 0, a1 / 2 and a2 / 2."""
    vectors = primitive_vectors("hexagonal", a)
    return Crystal(vectors, [[0.0, 0.0], vectors[0] / 2, vectors[1] / 2], masses, names)


def lieb(a, masses, names):
    """Lieb lattice (CuO2 plane): square lattice constant a, atoms at 0 and the two edge centres."""
    return Crystal(primitive_vectors("square", a), [[0.0, 0.0], [a / 2, 0.0], [0.0, a / 2]], masses, names)


def rocksalt(a, masses, names):
    """Rocksalt (NaCl): fcc lattice with cubic lattice constant a, atoms at 0 and (a/2, 0, 0)."""
    return Crystal(bravais3d.primitive_vectors("fcc", a), [[0.0, 0.0, 0.0], [a / 2, 0.0, 0.0]], masses, names)


def diamond(a, masses, names):
    """Diamond (Si, Ge, C) or zincblende with two species: fcc lattice, atoms at 0 and a/4 (1, 1, 1)."""
    return Crystal(bravais3d.primitive_vectors("fcc", a), [[0.0, 0.0, 0.0], [a / 4, a / 4, a / 4]], masses, names)


def hcp(a, c, masses, names):
    """Hexagonal close-packed (Mg, Zn, Ti): hexagonal lattice a, c, atoms at 0 and (a/2, a / (2 sqrt 3), c/2).

    For the ideal c/a = sqrt(8/3) all 12 neighbours are at distance a.
    """
    positions = [[0.0, 0.0, 0.0], [a / 2, a / (2 * np.sqrt(3)), c / 2]]
    return Crystal(bravais3d.primitive_vectors("hexagonal", a, c=c), positions, masses, names)


def bilayer(crystal, bonds, shift):
    """Two sheets of a crystal as one crystal, the second shifted by shift (Cartesian, d components), with bonds.

    The atoms of the second sheet come after the first (index mu + B), and each sheet keeps the
    given bonds; springs between the sheets are added by the caller. E.g. embed(honeycomb(...), 3)
    with shift (0, 0, c) is AA bilayer graphene, with the B - A offset added in the plane AB.
    """
    n_basis = len(crystal.positions)
    second = [(mu + n_basis, nu + n_basis, n, K) for mu, nu, n, K in bonds]
    stacked = Crystal(
        crystal.vectors,
        np.vstack([crystal.positions, crystal.positions + np.asarray(shift, dtype=float)]),
        np.tile(crystal.masses, 2),
        list(crystal.names) * 2,
    )
    return stacked, list(bonds) + second
