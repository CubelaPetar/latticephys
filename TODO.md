# TODO

Done: `Crystal` + bonds `(mu, nu, n, K)` with dimension d ≤ 3 and periodicity p ≤ d (molecules, chains in a plane, ribbons, layers, bilayers); 14 3D Bravais lattices, special points (2D, layers, cubic, fcc, bcc, tetragonal, hexagonal), structures (honeycomb, kagome, Lieb, rocksalt, diamond, hcp, bilayer); band structures along symmetry paths, (K_L, K_T) springs, group velocities, DOS and S(Q, ω), lattice dynamics with velocity Verlet, 2D and 3D crystal plots, 2D Brillouin zone plots, 1D and 2D animations (centred and scaled to the window, with bonds, optionally coloured by displacement); DOS, heat capacity, ⟨u²⟩ and sound velocities; JSON export.

Planned examples are in [`examples/IDEAS.md`](examples/IDEAS.md); the items below are the library features they need.

## Next

- [ ] verify that rule 6 in CLAUDE.md is realistic physically (some things can't be computed analytically)
- [ ] re-evaluate the README.md to properly reflect the library to make it users able to use the library efficiently. Emphasize is less on the examples than on the power of the lib
- [ ] check if the library captures all physical concepts described in chapters 1-4 of the book '2016 TheoretischeFestkörperphysik - Czycholl.pdf' in references/ (as far as possible).
- [ ] if possible, use our library to solve the exercises in chapter 2 of the book '2016 TheoretischeFestkörperphysik - Czycholl.pdf' and place them as examples (especially exe 2.4 is interesting)
- [ ] if possible, use our library to solve the exercises in chapter 4 of the book '2016 TheoretischeFestkörperphysik - Czycholl.pdf' and place them as examples (especially exe 4.8 and 4.9 are interesting)
- [ ] compare the calculated dispersion relation of this lib for NaCl with the graphs shown in 'Abb. 4.6' of '2016 TheoretischeFestkörperphysik - Czycholl.pdf'
- [ ] include electron physics
- [ ] include magnon physics

## Later

- [ ] `bloch_mode` with travelling waves (initial momenta), e.g. chiral phonons
- [ ] `thermal_state()`: random mode amplitudes and phases for temperature T
- [ ] `simulate()` with a force function instead of only Φ (anharmonic chains)
- [ ] angle-bending forces: molecule bending modes, ZA ∝ q² of graphene, flat TA and optical branches of Si
- [ ] 3D Brillouin-zone geometry and plot
- [ ] XYZ trajectory export (3D animation outside pygame)
- [ ] extended spring model with pre-tension (strained lattices)
- [ ] sparse force-constant matrices for systems beyond ~3000 atoms
- [ ] ribbon / strip builder in the library if more ribbon examples follow (now in `examples/lattices_2d/nanoribbon.py`)
