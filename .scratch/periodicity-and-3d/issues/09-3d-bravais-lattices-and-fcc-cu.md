# 09: 3D Bravais lattices and fcc Cu (tracer)

**What to build:** `primitive_vectors(lattice, a, b, c, alpha, beta, gamma)` for all 14 3D Bravais lattices (Setyawan–Curtarolo convention, missing parameters raise), 3D special points for cubic, fcc, bcc, tetragonal and hexagonal (others raise), and an fcc Cu example: dispersion along Γ–X–W–K–Γ–L, DOS and C_V(T). Spec: stories 31–33, 36.

**Blocked by:** 05 (both change `special_points`)

**Status:** done

- [x] Test, all 14 lattices: primitive cell volume = conventional volume / lattice points per conventional cell; centring vectors are integer combinations of the primitive vectors
- [x] Test: a_i · b_j = 2π δ_ij; every special point except Γ lies on the Brillouin-zone boundary (|q| = |q − G| for some reciprocal G)
- [x] Test: fcc with nearest-neighbour central springs equals the eigenvalues of the analytic 3×3 dynamical matrix at random q
- [x] Cu example listed in `examples/README.md`; `pytest`, `ruff check` and all examples green
