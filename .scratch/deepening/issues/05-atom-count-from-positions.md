# 05: Atom count from positions

**What to build:** the geometry (Crystal, neighbour bonds, periodic supercell, crystal plot) counts atoms from positions instead of masses, so it no longer depends on a phonon quantity. Masses stay a required field. Spec: `../spec.md` (story 27).

**Blocked by:** None (can start immediately)

**Status:** done

- [x] No atom count from masses left in the geometry package
- [x] Test: a Crystal with irrelevant masses gives the right bonds and supercell size
- [x] `pytest`, `ruff check`, all examples headless
