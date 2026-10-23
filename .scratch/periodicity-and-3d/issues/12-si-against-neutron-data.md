# 12: Si against neutron data

**What to build:** A `diamond(a, masses, names)` structure and a Si example with central springs to first and second neighbours: K₁ fixed by the Γ optical frequency 15.53 THz, K₂ fitted to TA(X) = 4.49 THz, dispersion along Γ–X–W–K–Γ–L with the measured points of `docs/research/si-phonon-frequencies.md` (a = 5.431 Å). Spec: stories 34, 38, 39.

**Blocked by:** 09

**Status:** ready-for-agent

- [ ] Test: optical frequency at Γ equals √(8K₁/3M); no negative frequency on a q-grid
- [ ] Example plots model against measured points at Γ, X, L
- [ ] Docstring: data are secondary-source (Savrasov 1996 quoting Nilsson & Nelin 1972); short-range springs miss the flat TA branch
- [ ] Listed in `examples/README.md`; `pytest`, `ruff check` and all examples green
