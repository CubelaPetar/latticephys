# 05: Graphene layer example

**What to build:** `special_points` and `k_path` work for a layer (2D labels Γ, M, K in its plane), and a new example shows the phonons of a graphene layer with (K_L, K_T) springs along Γ–M–K–Γ: ZA, TA, LA, ZO, TO, LO. The docstring explains the branch labels and that the real ZA ∝ q² needs bending forces. Spec: stories 20, 25.

**Blocked by:** 03

**Status:** ready-for-agent

- [x] Test: special points of an embedded hexagonal lattice equal the 2D ones with a zero z component
- [x] Example in its topic folder, constants at the top, physics docstring, `main()`, listed in `examples/README.md`
- [x] ZA is linear near Γ with (K_L, K_T); docstring says why and what is missing
- [x] `pytest`, `ruff check` and all examples green
