# 02: Remove experiment verdicts and example-only numbers from the prose

**What to build:** remove the claims below from docstrings; where the example prints the compared values, the docstring says what is compared, not how well it agrees. Spec: `../spec.md`.

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [ ] `crystals_3d/copper.py:6,15` "TA(X) matches, LA(L) within 1 %, TA(L) ~6 % high", "lowest mode ~0.5 THz", "log slope > 3 at 8 K"
- [ ] `crystals_3d/silicon.py:9` "LA(X) 3 % high, TO(X) 16 % high"
- [ ] `layers/bilayer_graphene.py:16` the cm⁻¹ values of the fit
- [ ] `layers/graphene_layer.py:16` "LO/TO near the Raman G band"
- [ ] `molecules/co2_modes.py:3` agreement with the measured 1340 / 2349 cm⁻¹
- [ ] `lattices_2d/nanoribbon.py:6` "up to about 70 % on the outer 25 %"
- [ ] `lattices_2d/structure_factor_2d.py:9` "b/m differ 8 %, 99.5 % extinct"
- [ ] `pytest`, `ruff check`, all examples headless
