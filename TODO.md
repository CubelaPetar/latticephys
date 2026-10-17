# TODO

## Library for lattice systems

- [x] band structures along symmetry points in Brillouin zones (1D, 2D in-plane)
- [x] solving coupled differential equations for lattice dynamics (1D, velocity Verlet)
- [x] drawing a lattice / structure and simulating its oscillations (1D, 2D)

## Next

- [ ] out-of-plane modes (3D displacements, e.g. graphene ZA branch)
- [ ] `animate()`: optionally colour atoms by displacement, so small waves are visible; draw bonds

## Later

- [ ] example: spring-dimerized chain (equal masses, k1 ≠ k2), phonon SSH analogue
- [ ] example: second-neighbour springs in 1D `(0, 0, 2, k2)`, unstable for negative k2
- [ ] example: energy error vs. time step, Verlet vs. explicit Euler
- [ ] Fermi-Pasta-Ulam-Tsingou chain: `simulate()` with a force function instead of only K
- [ ] example: molecules with unequal springs (e.g. acetylene) via `force_constants`
