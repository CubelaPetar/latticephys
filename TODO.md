# TODO

Done: `Crystal` + bonds `(mu, nu, n, k)` in any dimension; 1D and 2D (in-plane) band structures along symmetry paths, lattice dynamics with velocity Verlet, 2D crystal and Brillouin zone plots, 1D and 2D animations; DOS, heat capacity, ⟨u²⟩ and sound velocities; JSON export.

## Next

- [ ] out-of-plane modes (3D displacements, e.g. graphene ZA branch)
- [ ] `animate()`: optionally colour atoms by displacement, so small waves are visible; draw bonds

## Later

- [ ] example: anisotropic wave fronts on the square lattice with 2nd neighbours (phonon focusing)
- [ ] example: isotope defect in 2D, localized mode above the band
- [ ] chiral phonons at K in hBN: `bloch_mode` with travelling waves (initial momenta)
- [ ] strained lattices: pre-tension gives transverse stiffness; needs an extended spring model
- [ ] ribbon / strip builder in the library if more ribbon examples follow (now in `nanoribbon.py`)
- [ ] example: spring-dimerized chain (equal masses, k1 ≠ k2), phonon SSH analogue
- [ ] example: second-neighbour springs in 1D `(0, 0, 2, k2)`, unstable for negative k2
- [ ] example: energy error vs. time step, Verlet vs. explicit Euler
- [ ] Fermi-Pasta-Ulam-Tsingou chain: `simulate()` with a force function instead of only K
- [ ] example: molecules with unequal springs (e.g. acetylene) via `force_constants`
