# Examples

Run from the repo root, e.g. `.venv/bin/python examples/wave_packet.py`. Parameters are constants at the top of each file; the docstring explains the physics. Animations (pygame) loop until the window is closed; `SECONDS_PER_PERIOD` (or `SECONDS` in `wave_2d.py`) sets the playback speed, `TIMESTEP` only the accuracy.

## Molecules

| Example | Shows | Uses |
|---|---|---|
| `co2_modes.py` | CO₂ normal modes in cm⁻¹ | `normal_modes` |
| `co2_animation.py` | CO₂ symmetric stretch (pygame) | `simulate`, `animate` |

## 1D chains

| Example | Shows | Uses |
|---|---|---|
| `chain_animation.py` | Periodic chain, one atom displaced (pygame) | `chain_force_constants`, `simulate`, `animate` |
| `chain_dispersion.py` | ω(q) of monatomic and diatomic chain | `dispersion` |
| `diatomic_mode_animation.py` | One acoustic or optical Bloch mode (pygame) | `periodic_supercell`, `bloch_mode`, `simulate`, `animate` |
| `finite_chain_modes.py` | Fixed-chain modes lie on ω(q) at q_j = jπ/((N+1)a) | `normal_modes`, `dispersion` |
| `mass_defect.py` | Localized mode of a light impurity above the band | `normal_modes` |
| `wave_packet.py` | Packet moving with the group velocity (space–time plot) | `simulate` |
| `heat_capacity.py` | Density of states, C_V(T) vs. Debye and Einstein | `dispersion`, `units` |

## 2D lattices

| Example | Shows | Uses |
|---|---|---|
| `crystals_2d.py` | NaCl plane, graphene, hBN, triangular, kagome, CuO₂ (Lieb), AFM | `Crystal`, `plot_crystal` |
| `brillouin_zones.py` | First Brillouin zones of the five 2D Bravais lattices | `plot_brillouin_zone` |
| `phonons_2d.py` | In-plane band structures (triangular, square, graphene, hBN); floppy modes with 1st neighbours only | `neighbour_bonds`, `dispersion`, `k_path`, `plot_bands` |
| `mode_animation_2d.py` | One graphene Bloch phonon at Γ, M or K (pygame) | `periodic_supercell`, `bloch_mode`, `simulate`, `animate` |
| `wave_2d.py` | One displaced atom launches longitudinal and transverse waves (pygame) | `periodic_supercell`, `force_constants`, `simulate`, `animate` |
| `dos_2d.py` | Graphene DOS with van Hove peaks, C_V ∝ T² (2D Debye) | `q_grid`, `dispersion`, `heat_capacity` |
| `sound_velocities_2d.py` | Sound velocity against direction; c_L/c_T = √3 (Cauchy relation) | `sound_velocities` |
| `nanoribbon.py` | Zigzag graphene ribbon: bands coloured by edge weight | `Crystal` with vacuum, `dynamical_matrix` |
| `mermin_wagner.py` | ⟨u²⟩ against size: linear in 1D, logarithmic in 2D, finite in 3D | `mean_square_displacement` |
| `kagome_lieb.py` | Floppy modes of kagome and Lieb (CuO₂) lattices, rigid unit mode at M | `neighbour_bonds`, `dispersion` |

## Data export

`export_for_site.py [directory]` writes `diatomic_chain.json` and `graphene.json` (default directory `site_data/`) with `save_json`:

| Key | Content |
|---|---|
| `system`, `units` | name; Å, THz, amu |
| `crystal` | `vectors` (d×d), `positions` (B×d), `masses`, `names` |
| `dispersion` | `distance` along the path, `frequencies` (points × branches), `ticks` and `labels` of the k-points |
| `supercell` | `cells` per direction, `positions` and `names` of all atoms |
| `modes` | list of `point`, `branch`, `frequency`, `displacements` (atoms × d, max 1) |

Animate a mode as `positions + scale * displacements * cos(2π frequency t)`. Zero-frequency modes (translations) do not move.
