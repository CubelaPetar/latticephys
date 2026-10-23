"""Bonds between atoms of a Crystal: neighbour shells and periodic supercells.

A bond (mu, nu, n, K) couples atom mu in cell 0 to atom nu in cell n (p integers,
or an int for p = 1) with a parameter K, e.g. a spring constant. Its direction
comes from the crystal.
"""

import numpy as np

from latticephys.lattice.crystal import Crystal, cell_offsets
from latticephys.lattice.reciprocal import reciprocal_vectors


def neighbour_bonds(crystal, spring_constants, search=3):
    """Bonds of a crystal to its neighbour shells.

    spring_constants: bond parameter per shell, [K1] for nearest neighbours,
                      [K1, K2] to add second neighbours, ...; each entry is used
                      as given, e.g. a spring K or a pair (K_L, K_T)
    search:           cells searched along each lattice vector

    Shells are the distinct atom distances in ascending order. Returns a list
    of bonds (mu, nu, n, K), each bond once. Raises if search is too small to
    find all requested shells completely.
    """
    origin = (0,) * crystal.periodicity
    candidates = []
    for mu in range(len(crystal.masses)):
        for nu in range(len(crystal.masses)):
            for n in map(tuple, cell_offsets(2 * search + 1, crystal.periodicity) - search):
                # count each bond only from one end: (mu, 0) -- (nu, n) equals (nu, 0) -- (mu, -n)
                if nu < mu or (nu == mu and n <= origin):
                    continue
                candidates.append((mu, nu, n, np.linalg.norm(crystal.bond_vector(mu, nu, n))))

    shells = np.unique(np.round([c[3] for c in candidates], 6))[: len(spring_constants)]
    if len(shells) < len(spring_constants) or shells[-1] >= _search_radius(crystal, search) - 1e-6:
        raise ValueError(f"{len(spring_constants)} shells need more than search = {search} cells, increase search")

    bonds = []
    for mu, nu, n, distance in candidates:
        shell = np.flatnonzero(np.isclose(shells, distance, atol=1e-6))
        if len(shell):
            bonds.append((mu, nu, n, spring_constants[shell[0]]))
    return bonds


def _search_radius(crystal, search):
    """Every bond to a cell outside |n_i| <= search is at least this long.

    Along the unit normal b_i / |b_i| of the lattice planes a bond to cell n has
    the component (2 pi n_i + (r_nu - r_mu) . b_i) / |b_i|.
    """
    b = reciprocal_vectors(crystal.vectors)
    b_norm = np.linalg.norm(b, axis=1)
    spread = np.abs((crystal.positions[:, None, :] - crystal.positions[None, :, :]) @ b.T).max(axis=(0, 1))
    return np.min((2 * np.pi * (search + 1) - spread) / b_norm, initial=np.inf)  # no bound for a molecule


def periodic_supercell(crystal, bonds, cells):
    """Supercell of cells^p unit cells (cells along each lattice vector) as a Crystal, with its bonds.

    force_constants() of the result is the finite system with periodic
    boundaries. Atom mu of cell c (integer tuple, 0 <= c_i < cells) gets
    index ravel(c) * B + mu, with the last entry of c running fastest.
    """
    p, n_basis = crystal.periodicity, len(crystal.masses)
    shape = (cells,) * p
    offsets = cell_offsets(cells, p)

    supercell = Crystal(
        cells * crystal.vectors,
        (offsets @ crystal.vectors)[:, None, :] + crystal.positions[None, :, :],
        np.tile(crystal.masses, len(offsets)),
        list(crystal.names) * len(offsets),
    )

    supercell_bonds = []
    for c in offsets:
        for mu, nu, n, K in bonds:
            target = c + crystal.cell_index(n)
            i = np.ravel_multi_index(c, shape) * n_basis + mu
            j = np.ravel_multi_index(target % cells, shape) * n_basis + nu
            supercell_bonds.append((i, j, tuple(target // cells), K))
    return supercell, supercell_bonds
