# 06: Nanoribbon as p = 1, d = 2

**What to build:** The zigzag nanoribbon example becomes a crystal with one lattice vector in 2D: no vacuum length, no second lattice vector. `plot_crystal` (2D) draws p = 1 crystals (atoms along the ribbon, cell outline adapted). Results of the example stay the same. Spec: story 27.

**Blocked by:** 03

**Status:** ready-for-agent

- [ ] Nanoribbon example has no `VACUUM` constant and uses `vectors` of shape (1, 2)
- [ ] Bands and edge weights unchanged (compare before/after numbers in the ticket's commit message)
- [ ] `plot_crystal` draws a p = 1, d = 2 crystal
- [ ] `pytest`, `ruff check` and all examples green
