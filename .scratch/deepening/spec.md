# Deepening plan: Φ(q), modes, Bloch modes, band plots

Status: ready-for-agent

## Problem Statement

The phonon package works, but its core is spread thin. The Fourier-transformed force-constant matrix Φ(q) has no public name: the dispersion code imports a private function, its phase convention is undocumented although the structure factor and Bloch modes depend on it, and it mixes the spring block, the lattice Fourier transform and the sum-rule on-site term in one loop. The same "build D(q), diagonalize" step is written seven times, one example goes around the library to get eigenvectors, and another solves every q twice. The mode animations integrate numerically what is known exactly, and two examples copy the decoration of the band plot only to colour branches. Later physics packages (tight-binding electrons, magnons) need to know which parts are shared lattice physics and which are phonon-specific, and today that seam is invisible.

## Solution

Give Φ(q) a public home with a documented cell-phase convention and keep its three pieces separable. Route every eigenproblem of a crystal at many wave vectors through one function that returns frequencies and polarizations together, and derive atom weights from those polarizations. Return Bloch modes as complex patterns and give their exact motion in time. Let the band plot colour branches by a weight. Count atoms in the geometry from positions, not masses. Record the seam map for electrons and magnons in the TODO list, and release the interface changes as version 0.4.0.

## User Stories

1. As a library user, I want to compute Φ(q) of a crystal through a public function, so that I don't depend on private names.
2. As a library user, I want the same function to give the real Φ of a finite system when I give no wave vector, so that one name covers the force-constant matrix in both cases.
3. As a library user, I want the docstring of Φ(q) to state the phase convention (cell phase e^{iq·R_n}), so that I can interpret eigenvectors and compare with other codes.
4. As a library user, I want Φ(q) to raise when q has the wrong number of components or leaves the span of the lattice vectors, so that a wrong q never gives a silent result.
5. As a physicist, I want a test showing Φ(q + G) = Φ(q) for a crystal with several atoms, so that the convention cannot change silently.
6. As a maintainer, I want the spring block, the lattice Fourier transform and the sum-rule on-site term as separate pieces, so that the transform can later serve hopping and exchange blocks.
7. As a maintainer, I want dΦ/dq as its own private function instead of a flag that changes the return shape, so that each function has one job.
8. As a maintainer, I want the mass weighting M^-½ Φ M^-½ written once, so that the dynamical matrix, group velocities and normal modes cannot drift apart.
9. As a maintainer, I want a rule that private names may be imported between modules of one package but never from examples, tests or another package, so that "private" has a clear meaning.
10. As a library user, I want `phonon_modes(crystal, bonds, qs)` returning frequencies and polarizations for many wave vectors, so that I solve each q once and get everything.
11. As a library user, I want `dispersion` to agree with the frequencies of `phonon_modes`, so that both paths give the same bands.
12. As a maintainer, I want every eigenproblem of a crystal (dispersion, group velocities, atom weights, structure factor, mean square displacement, Bloch mode) to go through `phonon_modes`, so that the solve lives in one place.
13. As a library user, I want `atom_weights(vectors, dim)` giving |e_{μ,i}|² per atom and direction from polarizations, so that I can sum over edge atoms, one sublattice or the out-of-plane direction without solving again.
14. As a physicist, I want the atom weights of each mode to sum to 1, so that they read as fractions of the mode.
15. As a reader of the nanoribbon example, I want the edge weights computed from one solve, so that the example shows the efficient path.
16. As a reader of the graphene layer example, I want its out-of-plane weight from `atom_weights`, so that the example no longer goes around the library.
17. As a reader of the guide, I want the normal modes of a molecule shown as `phonon_modes` at q = 0, so that molecules use the same model-first call as crystals.
18. As a library user, I want `normal_modes(Phi, masses)` to stay for a raw Φ, so that walls, mass defects and later anharmonic models keep working.
19. As a library user, I want `bloch_mode` to return the complex pattern U on the supercell, so that I can build standing waves (Re U) and travelling waves (U e^{-iωt}).
20. As a library user, I want `mode_motion(omega, U, times)` giving the exact motion Re[U e^{-iωt}], so that mode animations have no time-step error.
21. As a reader of the mode animations (diatomic chain, graphene), I want them to use the exact motion, so that they show the mode without a numerical integrator.
22. As a physicist, I want a test showing that `simulate` on the supercell reproduces `mode_motion` of a Bloch mode, so that the integrator and the mode stay checked against each other.
23. As a future user of chiral phonons, I want the travelling wave available from the complex U, so that the hBN idea needs no new function.
24. As a maintainer of the export example, I want it to take the real part of U, so that the exported displacements stay as before.
25. As a library user, I want `plot_bands` to colour branches by a weight with a colour bar when I pass weights, so that I don't copy its decoration.
26. As a reader of the graphene layer and nanoribbon examples, I want them to use the coloured band plot, so that they are shorter.
27. As a future electron developer, I want the geometry to count atoms from positions, so that a crystal without meaningful masses still has the right number of atoms.
28. As a future electron or magnon developer, I want a TODO item with the seam map (shared: geometry, lattice Fourier transform, Hermitian eigenproblem per k, band plot, DOS broadening; per package: bond block, on-site term, eigenvalue meaning, statistics, masses or orbitals or spins, phase convention), so that I know what to move and what to write anew.
29. As a future electron developer, I want the open questions noted (phase convention after checking Vanderbilt, the atom-phase conversion only when a quantity needs it, the moves of the transform, DOS, heat capacity, optional masses, renaming the spring-constant argument), so that they are decided at the start of that work.
30. As a library user, I want the README table to list the new and changed functions, so that the front page matches the library.
31. As a library user, I want version 0.4.0 with a tag listing the changed interfaces, so that I see that upgrading needs code changes.

## Implementation Decisions

- **Φ(q):** `force_constants(crystal, bonds, q=None)`. Without q: the real Φ of the finite system (as now, the q = 0 sum). With q: the complex Φ(q) with cell phases e^{iq·R_n}; the docstring states the convention and the q + G periodicity, and the Bloch-mode and structure-factor docstrings refer to it. The q check (d components, in the span of the lattice vectors) moves from the dynamical matrix into this function.
- **Inner split of the spring module:** the spring block stays phonon-specific; a generic private lattice Fourier transform takes one block per bond (mu, nu, n, block) and adds only the terms between bonded atoms with the cell phase; the sum-rule on-site term is added separately for phonons. No phase-convention parameter now.
- **Gradient:** a private dΦ/dq function: d calls of the same transform with each block multiplied by i (R_n)_i, no on-site term. The `gradient` flag goes away.
- **Mass weighting:** one private helper in the modes module, used by the dynamical matrix, group velocities and normal modes.
- **Private-import rule:** one line under Conventions in CLAUDE.md: private names may be imported between modules of the same package, never from examples, tests or another package.
- **Modes:** `phonon_modes(crystal, bonds, qs)` in the modes module returns `omega` (len(qs), B d) and `vectors` (len(qs), B d, B d), modes as columns, frequencies as ω via `diagonalize`. It stays phonon-specific (no pluggable eigen-solver). `dispersion` returns its frequencies; group velocities, atom weights, structure factor, mean square displacement and Bloch mode use it. `dynamical_matrix` stays public.
- **Atom weights:** `atom_weights(vectors, dim)` returns |e_{μ,i}|² with shape (..., B d modes, B, d); it works on results like `density_of_states` does.
- **Finite systems:** `phonon_modes` at q = 0 is the model-first way to the modes of a molecule or supercell; `normal_modes(Phi, masses)`, `simulate(Phi, masses)` and `chain_force_constants` stay unchanged.
- **Bloch mode:** `bloch_mode(crystal, bonds, q, branch, cells)` returns `(omega, U)` with U complex, shape (B d cells^p,), U = e_μ e^{iq·R_c} / √m_μ with the largest component made real; ω = 0 stays exactly 0 for zero modes. `mode_motion(omega, U, times)` beside it returns Re[U e^{-iωt}], shape (len(times), N). The ω > 0 check and the amplitude scaling stay in the examples.
- **Band plot:** `plot_bands(distance, frequencies, corners, labels, ax=None, weights=None, weight_label=None)`: lines without weights; with weights (same shape as frequencies) coloured dots and a colour bar. The nanoribbon example switches to a k-path Γ–X with those labels.
- **Atom count:** the geometry (Crystal, bonds, supercell, crystal plot) counts atoms from positions; masses stay a required field.
- **Release:** version 0.4.0 in the project file and tag `v0.4.0` listing the changed interfaces (force_constants q, atom_weights, bloch_mode U) and the new ones (phonon_modes, mode_motion, plot_bands weights). Push only when the user says so.
- **Glossary** already updated during the grilling: Bond carries one model parameter, Lattice Fourier transform, Mode.

## Testing Decisions

- Tests check external behaviour through public functions, and each docstring states the analytic formula or independent calculation checked (repo convention). Private helpers are covered through the public functions.
- Φ(q): Φ(q + G) = Φ(q) for graphene; the q check raises (moved from the dynamical-matrix tests); existing finite-system tests unchanged. Prior art: the spring and dispersion tests.
- Gradient and mass weighting: the existing group-velocity tests against finite differences of the dispersion and the normal-mode tests.
- `phonon_modes`: frequencies equal `dispersion`; polarizations are orthonormal eigenvectors of D(q). Prior art: dispersion tests.
- `atom_weights`: sums to 1 per mode; the zone-edge diatomic case (only one sublattice moves) per atom; one direction case (out-of-plane weight of the flat layer).
- `bloch_mode` / `mode_motion`: the periodic-chain test becomes "`simulate` on the supercell equals `mode_motion`"; the zone-edge, misfit-q and zero-mode tests take the real part.
- `plot_bands` with weights: draws a colour bar and one dot per frequency. Prior art: the plot tests.
- Atom count: a Crystal whose masses are irrelevant still gives the right bonds and supercell size from positions.
- Every changed example and the guide snippet run headless through the examples and docs tests; `pytest`, `ruff check` before each commit.

## Out of Scope

- Moving any code to a shared module or to the geometry package (first ticket of the electron work).
- A phase-convention choice in the transform, the atom-phase conversion, optional masses, renaming the spring-constant argument.
- Changing `k_path`, `chain_force_constants`, `normal_modes`, `simulate`, the K fits and the model-versus-measured prints of the examples.
- Moving density of states, heat capacity, intensity map or export.
- Tests for `animate`; an ADR.

## Further Notes

Origin: the architecture review of 2026-10-24 (candidates 1–6; 7, `animate` untested, excluded). Order of the work: Φ(q) → phonon_modes → Bloch mode and band plot → release; the atom count is independent.
