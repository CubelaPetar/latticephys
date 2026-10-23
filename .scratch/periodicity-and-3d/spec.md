# Spec: periodicity and 3D crystals

Status: ready-for-agent

Vocabulary: see `GLOSSARY.md` (dimension, periodicity, layer, molecule, cell index, bond, branch label, spring constant, force-constant matrix). Decision: `docs/adr/0001-periodicity-separate-from-dimension.md`. Experimental data: `docs/research/si-phonon-frequencies.md`, `docs/research/bilayer-graphene-interlayer-modes.md`.

## Problem Statement

latticephys handles 1D and 2D bulk crystals only. A physicist cannot compute phonons of 3D crystals (Cu, NaCl, Si), nor of systems that are periodic in fewer directions than they move in: layers and bilayers with out-of-plane motion, ribbons, chains moving in a plane. Today such systems need workarounds (a long empty lattice vector with vacuum, as in the nanoribbon example; a dummy lattice vector for molecules) that introduce arbitrary lengths and meaningless wave-vector components.

Before extending, two inconsistencies need to go: every animation example centres and scales its atoms by hand, and the letter k means both spring constant and wave vector while K means both spring constant and force-constant matrix.

## Solution

Four steps, in this order:

1. **`animate()` cleanup**: `animate` centres and scales the system to its window itself and can draw bonds; examples no longer do this by hand.
2. **Symbol rename**: spring constants are K or (K_L, K_T), the force-constant matrix is Φ, the dynamical matrix stays D, and k, q only ever mean wave vectors.
3. **Periodicity**: a crystal has dimension d and periodicity p ≤ d. Molecules (p = 0), ribbons and chains in a plane (p = 1), layers and bilayers (p = 2, d = 3) become first-class crystals.
4. **3D**: named 3D Bravais lattices with standard special points, diamond and rocksalt structures, a 3D crystal plot, and examples for Cu, NaCl and Si (compared with neutron data).

## User Stories

### Animation

1. As a physicist, I want `animate` to centre my system in the window, so that I do not subtract the mean position in every example.
2. As a physicist, I want `animate` to choose a scale that fits the whole system into the window, so that I do not compute pixels per Å by hand.
3. As a physicist, I want to still pass an explicit scale if I want one, so that I can compare two animations at the same scale.
4. As a physicist, I want `animate` to draw the bonds between atoms, so that I see the structure and how springs stretch.
5. As a physicist, I want bonds to follow the moving atoms, so that compressions and shears are visible.
6. As a maintainer, I want the centring and scaling in one tested function, so that the logic is not duplicated across examples.

### Symbols

7. As a physicist reading the code, I want spring constants named K (or K_L, K_T), so that k and q are free for wave vectors as in textbooks.
8. As a physicist reading the code, I want the force-constant matrix named Φ, so that it is not confused with a spring constant K.
9. As a physicist, I want docstrings, examples, tests and `CLAUDE.md` to use the same symbols, so that I learn one convention.

### Periodicity

10. As a physicist, I want a crystal with `vectors` of shape (p, d), so that a layer has two lattice vectors in 3D space without a made-up third one.
11. As a physicist, I want a molecule to be a crystal with no lattice vectors (p = 0), so that I do not invent a dummy cell.
12. As a physicist, I want the cell index n of a bond to have p entries (an int still allowed for p = 1, `()` for a molecule), so that bonds match the lattice.
13. As a physicist, I want a bond with a cell index of the wrong length to raise an error, so that typos do not silently give wrong physics.
14. As a physicist, I want to lift a 2D structure (e.g. honeycomb) into 3D with `embed`, so that I get a graphene layer from the existing structure functions.
15. As a physicist, I want to lift a 1D chain into 2D with `embed`, so that I can study transverse motion of a chain.
16. As a physicist, I want wave vectors to stay Cartesian d-vectors, so that q of a layer is a 3-vector in its plane.
17. As a physicist, I want a q with a component outside the plane of the lattice vectors to raise an error, so that a meaningless q_z is never silently ignored.
18. As a physicist, I want `reciprocal_vectors` to work for p < d, so that a layer has in-plane reciprocal vectors with a_i · b_j = 2π δ_ij.
19. As a physicist, I want `q_grid` to give n^p wave vectors, so that Brillouin-zone averages of a layer use only its two periodic directions.
20. As a physicist, I want `special_points` and `k_path` to work for a layer, so that I plot its bands along Γ–M–K–Γ.
21. As a physicist, I want `neighbour_bonds` to search only along the p lattice vectors, so that shells of a layer or ribbon come out right.
22. As a physicist, I want `periodic_supercell`, `bloch_mode`, DOS, heat capacity, ⟨u²⟩, S(Q, ω) and group velocities to work for p < d, so that every existing tool applies to layers and ribbons.
23. As a physicist, I want `sound_velocities` to reject directions outside the lattice plane, so that I do not get a velocity along a non-periodic direction.
24. As a physicist, I want a graphene layer with central springs to have exactly B zero out-of-plane branches, so that I see why bending stiffness is needed.
25. As a physicist, I want a graphene layer example with (K_L, K_T) showing ZA, TA, LA, ZO, TO, LO, so that I see out-of-plane modes.
26. As a physicist, I want a bilayer graphene example with interlayer springs fitted to the measured shear (31 cm⁻¹) and breathing (81 cm⁻¹) modes, so that I see loosely coupled sheets.
27. As a physicist, I want the nanoribbon example as a p = 1, d = 2 crystal, so that it needs no vacuum length.
28. As a physicist, I want existing molecule results (CO₂) unchanged, so that the change breaks nothing.
29. As a physicist, I want `plot_crystal` to choose 2D or 3D drawing from d, so that one function draws any crystal.
30. As a physicist, I want `wigner_seitz_cell`, `plot_brillouin_zone` and `animate` to raise for systems they cannot show, so that I get a clear message instead of a wrong picture.

### 3D

31. As a physicist, I want primitive vectors for all 14 3D Bravais lattices, so that I build any 3D crystal in one line.
32. As a physicist, I want a lattice to require the parameters it needs (b, c, angles), so that a missing parameter raises instead of defaulting.
33. As a physicist, I want the standard special points (Setyawan–Curtarolo labels) for cubic, fcc, bcc, tetragonal and hexagonal lattices, and a clear error for the others, so that my band paths match the literature.
34. As a physicist, I want `diamond` and `rocksalt` structure functions, so that Si and NaCl are one line each.
35. As a physicist, I want a 3D crystal plot with atoms and the cell edges, so that I check my structure visually.
36. As a physicist, I want an fcc Cu example with dispersion along Γ–X–W–K–Γ–L, DOS and C_V(T), so that I see a complete 3D phonon calculation.
37. As a physicist, I want a NaCl example with nearest and second-neighbour springs, so that I see acoustic and optical branches of a two-atom 3D crystal, with a note that LO–TO splitting needs Coulomb forces.
38. As a physicist, I want a Si example with central springs to first and second neighbours, K₁ fixed by ω_Γ = 15.53 THz and K₂ fitted to TA(X), plotted against neutron data, so that I see how far a simple model gets.
39. As a physicist, I want the Si example to say that short-range springs miss the flat TA branch, so that I do not over-trust the fit.

### Process

40. As a maintainer, I want every new physical result checked against an analytic formula in a test, so that rule 6 of `CLAUDE.md` holds.
41. As a maintainer, I want every new example listed in `examples/README.md` and run headless in the tests, so that examples cannot rot.
42. As a maintainer, I want `TODO.md` *Later* to list XYZ trajectory export, bending forces and a 3D Brillouin-zone plot, so that postponed work is not lost.

## Implementation Decisions

- **Order**: animate cleanup, symbol rename, periodicity, 3D. Each step leaves `pytest`, `ruff check` and all examples green and is committed separately.
- **animate**: a new pure function `screen_coordinates(positions, size)` in the display module centres positions on the window and scales them to fit with a margin, keeping the aspect ratio. `animate` uses it unless an explicit scale is given. `animate` takes an optional list of atom-index pairs to draw as bonds (lines between current positions). Examples drop their own centring and scaling.
- **Symbols**: bond tuple `(mu, nu, n, K)`; transverse pair `(K_L, K_T)`; force-constant matrix Φ in docstrings and as variable name `Phi`; dynamical matrix D. Function names (`force_constants`, `dynamical_matrix`, `neighbour_bonds(crystal, spring_constants, ...)`) stay. `CLAUDE.md` conventions are updated to match.
- **Crystal**: `vectors` shape (p, d), `positions` (B, d); `dim` is d, a new property `periodicity` is p. Validation: p ≤ d, lattice vectors linearly independent. A molecule has `vectors` of shape (0, d).
- **Cell index**: p entries; int allowed for p = 1; `()` for p = 0. Wrong length raises.
- **embed(crystal, d)**: new function in the crystal module; appends zero components to vectors and positions; raises if d is smaller than the current dimension.
- **Reciprocal vectors** for (p, d): b = 2π (A Aᵀ)⁻¹ A, which reduces to the current formula for p = d.
- **Wave vectors**: Cartesian, d components. The dynamical matrix (and everything that takes q) raises if q has a component outside the span of the lattice vectors (relative tolerance). `q_grid` returns n^p points.
- **Out of the 2D-only tools**: `wigner_seitz_cell`, `brillouin_zone`, `plot_brillouin_zone` only for p = d = 2; `animate` only for d ≤ 2; all raise otherwise. `plot_crystal` dispatches on d (2 or 3).
- **Special points**: one `special_points(lattice, vectors)`; the label set is chosen by the lattice name and p (2D "hexagonal" Γ M K; 3D "hexagonal" Γ M K A L H). For a layer the 2D labels apply in its plane. 3D labels: cubic Γ X M R; fcc Γ X W K L U; bcc Γ H N P; tetragonal Γ X M Z R A; hexagonal Γ M K A L H. Vectors are assumed to be those of `primitive_vectors`. The other nine lattices have parameter-dependent label sets; `special_points` raises for them (give k-points explicitly).
- **3D Bravais lattices**: new module next to the 2D one, `primitive_vectors(lattice, a, b=None, c=None, alpha=None, beta=None, gamma=None)` (angles in degrees, as in 2D) for all 14 lattices, primitive vectors in the Setyawan–Curtarolo convention (Comput. Mater. Sci. 49, 299 (2010)); e.g. fcc a/2·(0,1,1), (1,0,1), (1,1,0); bcc a/2·(−1,1,1), (1,−1,1), (1,1,−1); hexagonal a(1,0,0), a(1/2, √3/2, 0), c(0,0,1). Missing parameters raise.
- **Structures**: `diamond(a, masses, names)` and `rocksalt(a, masses, names)` in the existing structures module, whose docstring drops "2D".
- **3D plot**: matplotlib 3D axis, atoms as spheres (scatter) coloured per basis atom, primitive cell edges, same arguments as the 2D `plot_crystal`.
- **Export**: JSON format unchanged; it already stores d components.
- **Si model**: central springs, first and second shells, `neighbour_bonds(si, [K1, K2])`. Data: neutron values from `docs/research/si-phonon-frequencies.md` (secondary source: Savrasov 1996 Table I quoting Nilsson & Nelin 1972), a = 5.431 Å.
- **Bilayer**: AB-stacked, interlayer distance 3.356 Å (graphite), interlayer springs fitted to shear 31 cm⁻¹ and breathing 81 cm⁻¹; the docstring notes the published α⊥ would give ≈ 89 cm⁻¹.

## Testing Decisions

- Tests check external behaviour through the public functions, never internals. Each docstring states the analytic formula or independent calculation (`CLAUDE.md` convention).
- **Seam 1, public library functions** (prior art: the 2D dispersion and bond tests):
  - `screen_coordinates`: result centred in the window, all points inside, aspect ratio kept.
  - Rename: all existing tests stay green.
  - Periodicity: reciprocal vectors of a (2, 3) lattice satisfy a_i · b_j = 2π δ_ij and lie in the plane; `q_grid` has n^p points; q outside the plane raises; wrong cell-index length raises; `embed` keeps all 2D frequencies of a crystal for in-plane q.
  - Layer with central springs: exactly B zero branches at every q, in-plane branches equal the 2D result.
  - AA bilayer with one vertical spring K⊥ per atom pair: breathing mode at Γ, ω² = 2K⊥/M.
  - Molecule with p = 0: same modes as with the current workaround.
  - 3D lattices: for all 14, the primitive cell volume equals the conventional volume divided by the lattice points per conventional cell (1, 2 or 4), and the centring vectors are integer combinations of the primitive vectors; a_i · b_j = 2π δ_ij; each special point lies on the Brillouin-zone boundary (|q| = |q − G| for some reciprocal G), Γ excepted.
  - fcc nearest-neighbour central springs: dispersion equals the eigenvalues of the analytic 3×3 dynamical matrix at random q.
  - Rocksalt with first and second shells: along Γ–X the closed-form 2×2 blocks per polarization.
  - Diamond with first and second shells: optical frequency at Γ equals √(8K₁/3M); no negative frequency on a q-grid.
- **Seam 2, `test_examples.py`**: runs every example headless; covers new examples, the migrated nanoribbon and the 3D `plot_crystal`.

## Out of Scope

- 3D Brillouin-zone geometry and plot (`TODO.md` *Later*).
- 3D animation; XYZ trajectory export (`TODO.md` *Later*).
- Bending (angle) forces, so ZA is linear not ∝ q² (`TODO.md` *Later*).
- Long-range Coulomb forces (LO–TO splitting in NaCl).
- Models beyond second neighbours for Si (flat TA branch).
- Checking that the Phonon Playground web page can show 3D data (outside this repo).
- Electrons and magnons.

## Further Notes

- Experimental Si values are secondary-source (original table paywalled); the example docstring says so.
- Bilayer interlayer distance is the graphite value; no bilayer measurement was found.
- Rule 6: each step is checked against the `CLAUDE.md` list before its commit.
