# TODO

Done: `Crystal` + bonds `(mu, nu, n, K)` with dimension d ≤ 3 and periodicity p ≤ d (molecules, chains in a plane, ribbons, layers, bilayers); 14 3D Bravais lattices, special points (2D, layers, cubic, fcc, bcc, tetragonal, hexagonal), structures (honeycomb, kagome, Lieb, rocksalt, diamond, hcp, bilayer); Φ(q) with cell phases, band structures along symmetry paths (optionally coloured by atom weights), (K_L, K_T) springs, group velocities, DOS and S(Q, ω), lattice dynamics with velocity Verlet and exact Bloch-mode motion, 2D and 3D crystal plots, 2D Brillouin zone plots, 1D and 2D animations (centred and scaled to the window, with bonds, optionally coloured by displacement); DOS, heat capacity, ⟨u²⟩ and sound velocities; JSON export.

Planned examples are in [`examples/IDEAS.md`](examples/IDEAS.md); the items below are the library features they need.

## Next

- [ ] verify that rule 6 in CLAUDE.md is realistic physically (some things can't be computed analytically)
- [ ] review the CLAUDE.md rules as a whole (after the deepening plan)
- [ ] isnt src/latticephys/lattice/structures.py already an example?
- [ ] check if the library captures all physical concepts described in chapters 1-4 of the book '2016 TheoretischeFestkörperphysik - Czycholl.pdf' in references/ (as far as possible).
- [ ] include electron physics (tight binding)
  - shared with phonons: geometry, lattice Fourier transform, Hermitian eigenproblem per k, band plot, DOS broadening
  - per package: bond block, on-site term, meaning of the eigenvalues, statistics, masses / orbitals / spins, phase convention
  - decide first: phase convention (check Vanderbilt), atom-phase conversion only where a quantity needs it, moving the transform, DOS and heat capacity to a shared module, optional masses, renaming the spring-constant argument
- [ ] include magnon physics

## Later

- [ ] `thermal_state()`: random mode amplitudes and phases for temperature T
- [ ] `simulate()` with a force function instead of only Φ (anharmonic chains)
- [ ] angle-bending forces: molecule bending modes, ZA ∝ q² of graphene, flat TA and optical branches of Si
- [ ] long-range Coulomb forces (Ewald sum): LO–TO splitting in ionic crystals such as NaCl
- [ ] 3D Brillouin-zone geometry and plot
- [ ] XYZ trajectory export (3D animation outside pygame)
- [ ] extended spring model with pre-tension (strained lattices)
- [ ] sparse force-constant matrices for systems beyond ~3000 atoms
- [ ] ribbon / strip builder in the library if more ribbon examples follow (now in `examples/lattices_2d/nanoribbon.py`)
