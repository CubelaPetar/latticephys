# 06: Electron TODO, README, v0.4.0

**What to build:** one `TODO.md` item for the electron work with the seam map (shared: geometry, lattice Fourier transform, Hermitian eigenproblem per k, band plot, DOS broadening; per package: bond block, on-site term, eigenvalue meaning, statistics, masses or orbitals or spins, phase convention) and the open questions (phase convention after checking Vanderbilt; atom-phase conversion only when a quantity needs it; moving the transform, density of states, heat capacity; optional masses; renaming the spring-constant argument). README table lists `phonon_modes`, `mode_motion` and the changed functions. Version 0.4.0 and tag `v0.4.0` listing the changed interfaces (force_constants q, atom_weights, bloch_mode U, dynamical_matrix moved to `phonons.modes`) and new ones (phonon_modes, mode_motion, plot_bands weights). Spec: `../spec.md` (stories 28–31).

**Blocked by:** 01, 02, 03, 04, 05

**Status:** ready-for-agent

- [ ] TODO item written, short
- [ ] README table and guide match the code; docs test passes
- [ ] Version 0.4.0, tag `v0.4.0`; push only when the user says so
- [ ] `pytest`, `ruff check`, all examples headless
