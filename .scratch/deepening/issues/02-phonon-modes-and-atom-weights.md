# 02: phonon_modes for all solves, atom weights from polarizations

**What to build:** `phonon_modes(crystal, bonds, qs)` in the modes module returns `omega` (len(qs), B d) and `vectors` (len(qs), B d, B d), modes as columns, ω via `diagonalize`; it stays phonon-specific. `dispersion`, group velocities, structure factor, mean square displacement and Bloch mode use it, so no solve of D(q) is written anywhere else. `atom_weights(vectors, dim)` returns |e_{μ,i}|² per atom and direction, shape (..., B d, B, d). The nanoribbon and graphene layer examples solve each q once and take their weights from `atom_weights`; the guide shows a molecule's normal modes as `phonon_modes` at q = 0. `normal_modes(Phi, masses)` stays for a raw Φ. Spec: `../spec.md` (stories 10–18).

**Blocked by:** 01

**Status:** done

- [x] Test: `phonon_modes` frequencies equal `dispersion`; polarizations are orthonormal eigenvectors of D(q)
- [x] Test: atom weights sum to 1 per mode; diatomic zone edge moves one sublattice; flat layer has out-of-plane weight only in the Z modes
- [x] No `diagonalize(dynamical_matrix(...))` left outside `phonon_modes` (library and examples)
- [x] Nanoribbon and graphene layer examples use `phonon_modes` + `atom_weights`
- [x] Guide molecule snippet uses `phonon_modes` at q = 0
- [x] `pytest`, `ruff check`, all examples headless
