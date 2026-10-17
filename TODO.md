# TODO

## Library for lattice systems

- [x] band structures along symmetry points in Brillouin zones (1D, 2D in-plane)
- [x] solving coupled differential equations for lattice dynamics (1D, velocity Verlet)
- [ ] drawing a lattice / structure and simulating its oscillations (1D done)

## Next

- [ ] 2D dynamics and animation (supercell from `neighbour_bonds` + `force_constants`, pygame in 2D)
- [ ] out-of-plane modes (3D displacements, e.g. graphene ZA branch)

## Later

- [ ] example: spring-dimerized chain (equal masses, k1 ≠ k2), phonon SSH analogue
- [ ] example: second-neighbour springs in 1D `(0, 0, 2, k2)`, unstable for negative k2
- [ ] example: energy error vs. time step, Verlet vs. explicit Euler
- [ ] Fermi-Pasta-Ulam-Tsingou chain: `simulate()` with a force function instead of only K
- [ ] example: molecules with unequal springs (e.g. acetylene) via `force_constants`
