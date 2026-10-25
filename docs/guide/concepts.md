# Concepts

How the main objects look in code. The terms (dimension, periodicity, bond, cell index, ...) are defined in [`GLOSSARY.md`](../../GLOSSARY.md).

## Crystal

`Crystal(vectors, positions, masses, names)`: lattice vectors (p, d) as rows, Cartesian positions (B, d) of the basis atoms in Å, masses in amu. d is the dimension, p ≤ d the periodicity.

```python
import numpy as np
from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.structures import honeycomb

molecule = Crystal(np.zeros((0, 1)), [[0.0], [1.16], [2.32]], [15.999, 12.011, 15.999], ["O", "C", "O"])  # p = 0
chain = Crystal([[2.5]], [0.0, 1.25], [12.011, 15.999], ["C", "O"])                                       # p = d = 1
graphene = honeycomb(2.46, [12.011] * 2, ["C"] * 2)                                                       # p = d = 2
layer = embed(graphene, 3)                                                       # p = 2, d = 3: atoms move out of the plane
copper = Crystal(primitive_vectors("fcc", 3.615), [[0.0, 0.0, 0.0]], [63.546], ["Cu"])                    # p = d = 3
```

Ready-made structures in `lattice.structures`: `honeycomb`, `kagome`, `lieb`, `rocksalt`, `diamond`, `hcp`, and `bilayer` to stack two sheets. Bravais lattices: `bravais2d.primitive_vectors` (5) and `bravais3d.primitive_vectors` (14).

## Bonds

A bond `(mu, nu, n, K)` is a spring from atom mu in cell 0 to atom nu in cell n, listed once. The cell index n has p entries, an int for p = 1, `()` for a molecule. K is a spring constant (eV/Å²) or a pair (K_L, K_T), see [models](models.md).

```python
from latticephys.lattice.bonds import neighbour_bonds

molecule_bonds = [(0, 1, (), 100.0), (1, 2, (), 100.0)]
chain_bonds = [(0, 1, 0, 10.0), (1, 0, 1, 10.0)]      # C-O in the cell, O to the C of the next cell
graphene_bonds = neighbour_bonds(graphene, [25.0, 5.0])  # one K per neighbour shell
```

## Wave vectors and paths

q is a Cartesian d-vector (1/Å). `special_points(lattice, vectors)` gives the labelled points of the Brillouin zone, `k_path` a path through them, `q_grid` a uniform grid for the DOS and thermal properties.

```python
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

points = special_points("hexagonal", graphene.vectors)
qs, distance, corners = k_path([points[p] for p in "ΓMKΓ"])
frequencies = dispersion(graphene, graphene_bonds, qs) * OMEGA_TO_THZ  # (len(qs), B d)
print(dispersion(chain, chain_bonds, [0.0, np.pi / 2.5]) * OMEGA_TO_THZ)  # acoustic and optical at Γ and the zone edge
```

## Finite systems

For a molecule, or a supercell of a crystal (`periodic_supercell`), the normal modes are the phonon modes at q = 0:

```python
from latticephys.phonons.modes import phonon_modes

omega, vectors = phonon_modes(molecule, molecule_bonds, [np.zeros(molecule.dim)])
omega, vectors = omega[0], vectors[0]  # the only q
```

For a force-constant matrix Φ itself, from `force_constants(crystal, bonds)` or changed by hand (walls, mass defects), use `normal_modes(Phi, masses)`.

Degrees of freedom are ordered atom by atom, `[u_0x, u_0y, u_1x, ...]`, masses per degree of freedom `np.repeat(masses, d)`. `simulate` integrates the equations of motion and `animate` shows them (pygame).
