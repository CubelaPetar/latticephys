# 11: NaCl

**What to build:** A `rocksalt(a, masses, names)` structure and a NaCl example with central springs to first and second neighbours: acoustic and optical branches of a two-atom 3D crystal. The docstring notes that nearest-neighbour springs alone give TA = 0 along Γ–X, and that LO–TO splitting needs long-range Coulomb forces. Spec: story 37.

**Blocked by:** 09

**Status:** ready-for-agent

- [ ] Test: along Γ–X the frequencies equal the closed-form 2×2 blocks per polarization
- [ ] Example listed in `examples/README.md`; `pytest`, `ruff check` and all examples green
