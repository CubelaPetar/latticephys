# Examples

Run from the repo root, e.g. `.venv/bin/python examples/chains/wave_packet.py`. Parameters are constants at the top of each file; the docstring explains the physics. Planned examples are in [`IDEAS.md`](IDEAS.md). Animations (pygame) loop until the window is closed; `SECONDS_PER_PERIOD` (or `SECONDS` in `wave_2d.py`) sets the playback speed, `TIMESTEP` only the accuracy.

## Molecules (`molecules/`)

| Example | Shows | Uses |
|---|---|---|
| `co2_modes.py` | CO₂ normal modes in cm⁻¹ | `normal_modes` |
| `co2_animation.py` | CO₂ symmetric stretch (pygame) | `simulate`, `animate` |

## 1D chains (`chains/`)

| Example | Shows | Uses |
|---|---|---|
| `chain_animation.py` | Periodic chain, one atom displaced (pygame) | `chain_force_constants`, `simulate`, `animate` |
| `chain_dispersion.py` | ω(q) of monatomic and diatomic chain | `dispersion` |
| `diatomic_mode_animation.py` | One acoustic or optical Bloch mode (pygame) | `periodic_supercell`, `bloch_mode`, `simulate`, `animate` |
| `finite_chain_modes.py` | Fixed-chain modes lie on ω(q) at q_j = jπ/((N+1)a) | `normal_modes`, `dispersion` |
| `mass_defect.py` | Localized mode of a light impurity above the band | `normal_modes` |
| `wave_packet.py` | Packet moving with the group velocity (space–time plot) | `simulate` |
| `heat_capacity.py` | Density of states, C_V(T) vs. Debye and Einstein | `dispersion`, `density_of_states`, `heat_capacity` |
| `structure_factor.py` | Scattering intensity S(Q, ω) over three zones; acoustic branch extinct near Q = 2π/a | `structure_factor`, `density_of_states` |

## 2D lattices (`lattices_2d/`)

| Example | Shows | Uses |
|---|---|---|
| `crystals_2d.py` | NaCl plane, graphene, hBN, triangular, kagome, CuO₂ (Lieb), AFM | `Crystal`, `honeycomb`, `kagome`, `lieb`, `plot_crystal` |
| `brillouin_zones.py` | First Brillouin zones of the five 2D Bravais lattices | `plot_brillouin_zone` |
| `phonons_2d.py` | In-plane band structures (triangular, square, graphene, hBN); floppy modes with 1st neighbours only | `neighbour_bonds`, `dispersion`, `k_path`, `plot_bands` |
| `mode_animation_2d.py` | One graphene Bloch phonon at Γ, M or K (pygame) | `periodic_supercell`, `bloch_mode`, `simulate`, `animate` |
| `wave_2d.py` | One displaced atom launches longitudinal and transverse waves, atoms coloured by displacement (pygame) | `periodic_supercell`, `force_constants`, `simulate`, `animate` |
| `dos_2d.py` | Graphene DOS with van Hove peaks, C_V ∝ T² (2D Debye) | `q_grid`, `dispersion`, `density_of_states`, `heat_capacity` |
| `transverse_springs.py` | Square lattice: central springs leave Γ–X floppy, (K_L, K_T) springs give shear stiffness | `neighbour_bonds` with (K_L, K_T), `plot_bands` |
| `sound_velocities_2d.py` | Sound velocity against direction; c_L/c_T = √3 (Cauchy relation) | `sound_velocities` |
| `nanoribbon.py` | Zigzag graphene ribbon: bands coloured by edge weight | `Crystal` with one lattice vector (p = 1), `dispersion`, `atom_weights` |
| `mermin_wagner.py` | ⟨u²⟩ against size: linear in 1D, logarithmic in 2D, finite in 3D | `mean_square_displacement` |
| `kagome_lieb.py` | Floppy modes of kagome and Lieb (CuO₂) lattices, rigid unit mode at M | `kagome`, `lieb`, `neighbour_bonds`, `dispersion` |

## Layers (`layers/`)

| Example | Shows | Uses |
|---|---|---|
| `graphene_layer.py` | Graphene in 3D: ZA, TA, LA, ZO, TO, LO along Γ–M–K–Γ, coloured by out-of-plane weight; ZA linear without bending forces | `embed`, `neighbour_bonds` with (K_L, K_T), `special_points`, `k_path`, `dynamical_matrix`, `diagonalize` |
| `bilayer_graphene.py` | AB bilayer: interlayer springs fitted to the Raman shear (31 cm⁻¹) and breathing (81 cm⁻¹) modes; published constants give ≈ 89 cm⁻¹ | `embed`, `honeycomb`, `bilayer`, `neighbour_bonds`, `special_points`, `k_path`, `dispersion` |

## 3D crystals (`crystals_3d/`)

| Example | Shows | Uses |
|---|---|---|
| `copper.py` | fcc Cu, nearest neighbours fitted to LA(X): ω(q) along Γ–X–W–K–Γ–L against neutron data at X and L, DOS, C_V ∝ T³ | `primitive_vectors` (3D), `special_points`, `neighbour_bonds`, `dispersion`, `density_of_states`, `heat_capacity` |

## Data export (`export/`)

`export_for_site.py [directory]` writes `diatomic_chain.json` and `graphene.json` (default directory `site_data/`) with `save_json`:

| Key | Content |
|---|---|
| `system`, `units` | name; Å, THz, amu |
| `crystal` | `vectors` (d×d), `positions` (B×d), `masses`, `names` |
| `dispersion` | `distance` along the path, `frequencies` (points × branches), `ticks` and `labels` of the k-points |
| `supercell` | `cells` per direction, `positions` and `names` of all atoms |
| `modes` | list of `point`, `branch`, `frequency`, `displacements` (atoms × d, max 1) |

Animate a mode as `positions + scale * displacements * cos(2π frequency t)`. Zero-frequency modes (translations) do not move.
