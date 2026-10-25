# 01: Φ(q) public with a documented convention

**What to build:** `force_constants(crystal, bonds, q=None)` is the one public name for Φ: without q the real Φ of the finite system (as now), with q the complex Φ(q) with cell phases e^{iq·R_n}, its docstring stating the convention and the q + G periodicity. The q check (d components, in the span of the lattice vectors) moves here from the dynamical matrix. Inside, the spring block, a generic private lattice Fourier transform (one block per bond, terms between bonded atoms only) and the sum-rule on-site term are separate pieces; dΦ/dq is its own private function (the same transform with blocks times i (R_n)_i, no on-site term) and the `gradient` flag goes away; M^-½ Φ M^-½ is one private helper in the modes module, used by the dynamical matrix, group velocities and normal modes. Spec: `../spec.md` (stories 1–9).

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [ ] `force_constants` with q returns Φ(q); without q unchanged; Bloch-mode and structure-factor docstrings refer to its convention
- [ ] Test: Φ(q + G) = Φ(q) for graphene
- [ ] Wrong-length or out-of-span q raises in `force_constants` (and so in `dynamical_matrix`)
- [ ] No caller outside the spring module uses the old private transform or a `gradient` flag
- [ ] Mass weighting written once; group-velocity and normal-mode tests pass unchanged
- [ ] CLAUDE.md Conventions: private names may be imported between modules of one package, never from examples, tests or another package
- [ ] `pytest`, `ruff check`, all examples headless
