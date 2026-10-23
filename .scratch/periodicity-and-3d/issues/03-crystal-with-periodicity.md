# 03: Crystal with periodicity p ≤ d (tracer through a layer)

**What to build:** A crystal has dimension d and periodicity p ≤ d: `vectors` of shape (p, d), a `periodicity` property. A graphene layer (honeycomb lifted into 3D with `embed`) gives correct phonons with in-plane wave vectors. Includes: cell-index length checked against p (int allowed for p = 1); reciprocal vectors b = 2π (A Aᵀ)⁻¹ A; any function taking q raises when q leaves the span of the lattice vectors; `q_grid` returns n^p points; `neighbour_bonds` searches only along the p lattice vectors; `embed(crystal, d)`. Spec: stories 10, 12–19, 21; ADR 0001.

**Blocked by:** 02

**Status:** done

- [x] `Crystal` validates p ≤ d and independent lattice vectors; p = d crystals behave exactly as before
- [x] Test: reciprocal vectors of a (2, 3) lattice satisfy a_i · b_j = 2π δ_ij and lie in its plane
- [x] Test: q with an out-of-plane component raises; cell index of wrong length raises
- [x] Test: `q_grid` of a layer has n² points
- [x] Test: graphene layer with central springs has exactly B zero branches at every in-plane q, and its other branches equal the 2D result
- [x] Test: `embed` keeps the 2D frequencies for in-plane q; raises for a smaller d
- [x] `CLAUDE.md` Crystal convention updated to `vectors` (p, d)
- [x] `pytest`, `ruff check` and all examples green
