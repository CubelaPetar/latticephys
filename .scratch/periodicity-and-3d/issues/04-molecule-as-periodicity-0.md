# 04: Molecule as p = 0

**What to build:** A molecule is a crystal with no lattice vectors (`vectors` of shape (0, d)) and bonds with cell index `()`, replacing the "vectors unused" workaround in the `Crystal` docstring. `chain.py` and the CO₂ examples stay unchanged. Spec: stories 11, 28.

**Blocked by:** 03

**Status:** ready-for-agent

- [ ] `Crystal` accepts `vectors` (0, d); `force_constants` works for it
- [ ] Test: a molecule with p = 0 gives the same normal modes as the same molecule built the old way
- [ ] `Crystal` docstring and `CLAUDE.md` describe a molecule as p = 0
- [ ] `pytest`, `ruff check` and all examples green
