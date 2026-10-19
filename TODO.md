# TODO

Done: `Crystal` + bonds `(mu, nu, n, k)` in any dimension; 1D and 2D (in-plane) band structures along symmetry paths, (k_L, k_T) springs, group velocities, DOS and S(Q, ω), lattice dynamics with velocity Verlet, 2D crystal and Brillouin zone plots, 1D and 2D animations (optionally coloured by displacement); DOS, heat capacity, ⟨u²⟩ and sound velocities; JSON export.

Planned examples are in [`examples/IDEAS.md`](examples/IDEAS.md); the items below are the library features they need.

## Next

- [ ] out-of-plane modes (3D displacements on a 2D lattice). With (k_L, k_T) springs the graphene ZA branch is linear; the flexural ω ∝ q² needs rotationally invariant bending terms (see angle-bending forces)
- [ ] `animate()`: draw bonds; centre and scale to the window itself (now done in each example)

## Later

- [ ] `bloch_mode` with travelling waves (initial momenta), e.g. chiral phonons
- [ ] `thermal_state()`: random mode amplitudes and phases for temperature T
- [ ] `simulate()` with a force function instead of only K (anharmonic chains)
- [ ] angle-bending forces (molecule bending modes)
- [ ] extended spring model with pre-tension (strained lattices)
- [ ] sparse force-constant matrices for systems beyond ~3000 atoms
- [ ] ribbon / strip builder in the library if more ribbon examples follow (now in `examples/lattices_2d/nanoribbon.py`)
