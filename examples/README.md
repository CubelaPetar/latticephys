# Examples

Run from the repo root, e.g. `.venv/bin/python examples/wave_packet.py`. Parameters are constants at the top of each file; the docstring explains the physics. Animations (pygame) loop until the window is closed, `SECONDS_PER_PERIOD` sets the playback speed.

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
| `diatomic_mode_animation.py` | One acoustic or optical Bloch mode (pygame) | `dynamical_matrix`, `simulate`, `animate` |
| `finite_chain_modes.py` | Fixed-chain modes lie on ω(q) at q_j = jπ/((N+1)a) | `normal_modes`, `dispersion` |
| `mass_defect.py` | Localized mode of a light impurity above the band | `normal_modes` |
| `wave_packet.py` | Packet moving with the group velocity (space–time plot) | `simulate` |
| `heat_capacity.py` | Density of states, C_V(T) vs. Debye and Einstein | `dispersion`, `units` |

## 2D lattices

| Example | Shows | Uses |
|---|---|---|
| `crystals_2d.py` | NaCl plane, graphene, hBN, triangular, kagome, CuO₂ (Lieb), AFM | `Crystal`, `plot_crystal` |
| `brillouin_zones.py` | First Brillouin zones of the five 2D Bravais lattices | `plot_brillouin_zone` |
