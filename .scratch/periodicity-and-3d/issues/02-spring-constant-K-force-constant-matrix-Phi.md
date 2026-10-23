# 02: Spring constant K, force-constant matrix Φ

**What to build:** One symbol per concept everywhere: spring constants are K or (K_L, K_T), bonds `(mu, nu, n, K)`, the force-constant matrix is Φ (variable `Phi`), the dynamical matrix stays D, and k, q only mean wave vectors. Function names stay. Mechanical rename across library, tests, examples and the `CLAUDE.md` conventions. Spec: stories 7–9; `GLOSSARY.md`.

**Blocked by:** 01 (touches the same example files)

**Status:** ready-for-agent

- [x] No spring constant named k, k_l, k_t or C remains in docstrings, variables or examples
- [x] No force-constant matrix named K remains; it is Φ / `Phi`
- [x] `CLAUDE.md` conventions use K, (K_L, K_T) and Φ
- [x] Behaviour unchanged: existing tests pass unmodified apart from names; all examples green
