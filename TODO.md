# TODO

## Library for lattice systems

- [ ] drawing band structures along symmetry points in Brillouin zones (1D done, 2D zones drawn)
- [x] solving coupled differential equations for lattice dynamics (1D, velocity Verlet)
- [ ] drawing a lattice / structure and simulating its oscillations (1D done)

## Next

- [ ] 2D phonons: nearest-neighbour search on a `Crystal` → bond list for `force_constants`, 2D `dispersion()` along symmetry paths
- [ ] 2D dynamics and animation

## Later

- [ ] example: spring-dimerized chain (equal masses, k1 ≠ k2), phonon SSH analogue
- [ ] example: second-neighbour springs `(0, 0, 2, k2)`; first make `dispersion()` warn on negative ω² instead of clipping
- [ ] example: energy error vs. time step, Verlet vs. explicit Euler
- [ ] Fermi-Pasta-Ulam-Tsingou chain: `simulate()` with a force function instead of only K
- [ ] example: molecules with unequal springs (e.g. acetylene) via `force_constants`
