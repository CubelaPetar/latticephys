# latticephys

Lattice physics with single-particle or harmonic models: crystals, their bonds, and the waves (phonons now) that run on them.

## Geometry

**Dimension**:
The number d of components of positions, displacements and wave vectors.
_Avoid_: space dimension, embedding dimension

**Periodicity**:
The number p ≤ d of lattice vectors of a crystal: 0 for a molecule, 1 for a chain or ribbon, 2 for a layer.
_Avoid_: lattice dimension, periodic dimension

**Molecule**:
A crystal with periodicity 0.

**Layer**:
A crystal with periodicity 2 and dimension 3, e.g. graphene, whose atoms can move out of its plane. A bilayer is a layer with two sheets in its basis.
_Avoid_: slab, 2D crystal

**Wave vector**:
k or q (1/Å); these letters never mean anything else. q for phonons.

**Cell index**:
The p integers n of a cell, located at n_1 a_1 + ... + n_p a_p; cell 0 is the home cell, and a molecule has only n = ().

**Bond**:
A link from atom mu in cell 0 to atom nu in cell n, listed once, carrying one model parameter (a spring constant for phonons).

**Mode**:
One solution at wave vector q: a frequency ω and a polarization (eigenvector of D). The modes of all q with the same index s form branch s.
_Avoid_: eigenmode; band (except in plots)

**Branch label**:
Polarization (L longitudinal, T transverse in-plane, Z out-of-plane) plus A acoustic or O optical, e.g. ZA.
_Avoid_: flexural branch (use ZA)

## Springs

**Spring constant**:
K (eV/Å²) of a central spring, or the pair K_L, K_T for longitudinal and transverse stiffness.
_Avoid_: k, C, D

**Shell**:
All neighbours at the same distance from an atom; springs are given per shell (K1 nearest, K2 second, ...).

**Force-constant matrix**:
Φ, the matrix of the harmonic potential V = ½ uᵀ Φ u; Φ(q) is its lattice Fourier transform.
_Avoid_: K (that is a spring constant)

**Lattice Fourier transform**:
The sum over bonds of one matrix block per bond times a phase: e^{iq·R_n} of the cell the bond reaches (phonons; repeats for q + G) or e^{iq·(R_n + r_ν − r_μ)} of the bond (repeats only up to a phase per atom). Gives Φ(q) for phonons, H(k) for electrons; both phases give the same frequencies or energies.
_Avoid_: Bloch sum (in electron texts a combination of atomic orbitals)

**Dynamical matrix**:
D = M^-1/2 Φ M^-1/2; its eigenvalues are ω².
