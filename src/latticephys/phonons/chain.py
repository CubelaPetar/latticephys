"""One-dimensional chain of masses coupled by nearest-neighbour springs."""

import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.springs import force_constants


def chain_force_constants(n, k, boundary="free"):
    """Force-constant matrix K of a chain of n atoms, V = 1/2 u^T K u.

    k:        spring constant between neighbours
    boundary: "free"     - open ends
              "fixed"    - end atoms attached by a spring k to immobile walls
              "periodic" - atom n-1 coupled to atom 0
    """
    if boundary not in ("free", "fixed", "periodic"):
        raise ValueError(f"unknown boundary condition: {boundary}")

    # the chain as one cell of length n; masses do not enter K
    chain = Crystal([[n]], np.arange(n), np.ones(n), [""] * n)
    bonds = [(i, i + 1, 0, k) for i in range(n - 1)]
    if boundary == "periodic":
        bonds.append((n - 1, 0, 1, k))

    K = force_constants(chain, bonds)

    if boundary == "fixed":
        K[0, 0] += k
        K[-1, -1] += k

    return K
