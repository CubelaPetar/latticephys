# TODO

Done: `Crystal` + bonds `(mu, nu, n, k)` in any dimension; 1D and 2D (in-plane) band structures along symmetry paths, (k_L, k_T) springs, group velocities, DOS and S(Q, ω), lattice dynamics with velocity Verlet, 2D crystal and Brillouin zone plots, 1D and 2D animations (centred and scaled to the window, with bonds, optionally coloured by displacement); DOS, heat capacity, ⟨u²⟩ and sound velocities; JSON export.

Planned examples are in [`examples/IDEAS.md`](examples/IDEAS.md); the items below are the library features they need.

## Next

- [ ] expand library to 3D physics
- [ ] include electron physics
- [ ] include magnon physics

## Later

- [ ] `bloch_mode` with travelling waves (initial momenta), e.g. chiral phonons
- [ ] `thermal_state()`: random mode amplitudes and phases for temperature T
- [ ] `simulate()` with a force function instead of only K (anharmonic chains)
- [ ] angle-bending forces (molecule bending modes)
- [ ] extended spring model with pre-tension (strained lattices)
- [ ] sparse force-constant matrices for systems beyond ~3000 atoms
- [ ] ribbon / strip builder in the library if more ribbon examples follow (now in `examples/lattices_2d/nanoribbon.py`)
