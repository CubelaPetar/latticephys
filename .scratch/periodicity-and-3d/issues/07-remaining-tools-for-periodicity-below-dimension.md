# 07: Remaining tools for p < d

**What to build:** Every existing tool either works for p < d or raises a clear error. Working: `periodic_supercell`, `bloch_mode`, DOS, heat capacity, ⟨u²⟩, S(Q, ω), group velocities. `sound_velocities` rejects directions outside the lattice plane. `wigner_seitz_cell`, `brillouin_zone`, `plot_brillouin_zone` (only p = d = 2) and `animate` (only d ≤ 2) raise otherwise. Spec: stories 22, 23, 30.

**Blocked by:** 03

**Status:** ready-for-agent

- [ ] Tests for a layer: `periodic_supercell` of n² cells, `bloch_mode` at an allowed q, DOS normalized to 3B branches, ⟨u²⟩ and S(Q, ω) run and match the embedded 2D crystal for in-plane motion
- [ ] Test: `sound_velocities` raises for an out-of-plane direction
- [ ] Tests: 2D-only tools raise for p < d or d = 3
- [ ] `pytest`, `ruff check` and all examples green

## Comments
- From the review of 03: `dynamical_matrix` raises for q outside the lattice span, but for a layer a scattering vector Q with a z part is meaningful (it shows the out-of-plane modes). `structure_factor` should give the phonon the in-plane part of Q and keep the full Q in the Q · e factor.
- From the review of 04: `cell_offsets(n, 0)` fails (`reshape(-1, 0)`), so `neighbour_bonds`, `q_grid`, `atom_positions`, `periodic_supercell` and `bloch_mode` raise an unclear numpy error for a molecule (p = 0). It should return one empty offset, shape (1, 0); then `q_grid` gives q = 0 and `neighbour_bonds` gives bonds with `()`.
