# 03: Complex Bloch mode and exact motion

**What to build:** `bloch_mode(crystal, bonds, q, branch, cells)` returns `(omega, U)` with U the complex pattern e_μ e^{iq·R_c} / √m_μ on the supercell (largest component real; ω = 0 exactly for zero modes). `mode_motion(omega, U, times)` beside it returns Re[U e^{-iωt}], shape (len(times), N). The diatomic-chain and graphene mode animations use `mode_motion` instead of `force_constants` + `simulate`; their ω > 0 check and amplitude scaling stay in the examples. The export example takes the real part. Spec: `../spec.md` (stories 19–24).

**Blocked by:** 02

**Status:** ready-for-agent

- [ ] Test: `simulate` on the supercell equals `mode_motion` of a Bloch mode (replaces the periodic-chain test)
- [ ] Zone-edge, misfit-q and zero-mode tests pass with the real part
- [ ] Both mode animations and the export example updated; examples README lists `mode_motion`
- [ ] `pytest`, `ruff check`, all examples headless
