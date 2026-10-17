# TODO

Done: 1D and 2D (in-plane) band structures along symmetry paths, lattice dynamics with velocity Verlet, 2D crystal and Brillouin zone plots, 1D and 2D animations.

## Cleanup first (before new features)

- [ ] same argument order for `dispersion(qs, masses, springs, a)` and `crystal_dispersion(crystal, bonds, qs)`; update examples and tests
- [ ] one bond format instead of three: `(mu, nu, n, k)` 1D lattices, `(mu, nu, n, k, r)` crystals, `(i, j, k, r)` finite systems (e.g. describe 1D chains with `Crystal` too)
- [ ] `bloch_mode()` for 1D as well; `diatomic_mode_animation.py` builds its initial state by hand

## Next

- [ ] out-of-plane modes (3D displacements, e.g. graphene ZA branch)
- [ ] `animate()`: optionally colour atoms by displacement, so small waves are visible; draw bonds

## Later

- [ ] example: spring-dimerized chain (equal masses, k1 ≠ k2), phonon SSH analogue
- [ ] example: second-neighbour springs in 1D `(0, 0, 2, k2)`, unstable for negative k2
- [ ] example: energy error vs. time step, Verlet vs. explicit Euler
- [ ] Fermi-Pasta-Ulam-Tsingou chain: `simulate()` with a force function instead of only K
- [ ] example: molecules with unequal springs (e.g. acetylene) via `force_constants`
