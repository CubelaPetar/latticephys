"""One-dimensional chain of masses coupled by nearest-neighbour springs."""

import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.springs import force_constants


def chain_force_constants(n, spring_constant, boundary="free"):
    """Force-constant matrix Phi of a chain of n atoms, V = 1/2 u^T Phi u.

    spring_constant: K between neighbours
    boundary:        "free"     - open ends
                     "fixed"    - end atoms attached by a spring K to immobile walls
                     "periodic" - atom n-1 coupled to atom 0
    """
    if boundary not in ("free", "fixed", "periodic"):
        raise ValueError(f"unknown boundary condition: {boundary}")

    # the chain as one cell of length n; masses do not enter Phi
    chain = Crystal([[n]], np.arange(n), np.ones(n), [""] * n)
    K = spring_constant
    bonds = [(i, i + 1, 0, K) for i in range(n - 1)]
    if boundary == "periodic":
        bonds.append((n - 1, 0, 1, K))

    Phi = force_constants(chain, bonds)

    if boundary == "fixed":
        Phi[0, 0] += K
        Phi[-1, -1] += K

    return Phi
