"""Common 2D crystal structures with several basis atoms."""

import numpy as np

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
