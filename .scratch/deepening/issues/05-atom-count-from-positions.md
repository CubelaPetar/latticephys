# 05: Atom count from positions

**What to build:** the geometry (Crystal, neighbour bonds, periodic supercell, crystal plot) counts atoms from positions instead of masses, so it no longer depends on a phonon quantity. Masses stay a required field. Spec: `../spec.md` (story 27).

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [ ] No atom count from masses left in the geometry package
- [ ] Test: a Crystal with irrelevant masses gives the right bonds and supercell size
- [ ] `pytest`, `ruff check`, all examples headless
