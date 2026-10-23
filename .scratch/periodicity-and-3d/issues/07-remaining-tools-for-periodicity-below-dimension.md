# 07: Remaining tools for p < d

**What to build:** Every existing tool either works for p < d or raises a clear error. Working: `periodic_supercell`, `bloch_mode`, DOS, heat capacity, ⟨u²⟩, S(Q, ω), group velocities. `sound_velocities` rejects directions outside the lattice plane. `wigner_seitz_cell`, `brillouin_zone`, `plot_brillouin_zone` (only p = d = 2) and `animate` (only d ≤ 2) raise otherwise. Spec: stories 22, 23, 30.

**Blocked by:** 03

**Status:** ready-for-agent

- [ ] Tests for a layer: `periodic_supercell` of n² cells, `bloch_mode` at an allowed q, DOS normalized to 3B branches, ⟨u²⟩ and S(Q, ω) run and match the embedded 2D crystal for in-plane motion
- [ ] Test: `sound_velocities` raises for an out-of-plane direction
- [ ] Tests: 2D-only tools raise for p < d or d = 3
- [ ] `pytest`, `ruff check` and all examples green
