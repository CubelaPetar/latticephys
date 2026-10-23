# 08: Bilayer graphene

**What to build:** Two loosely coupled graphene sheets as one layer crystal. A test checks the breathing mode of an AA bilayer; an example shows AB-stacked bilayer graphene with interlayer springs fitted to the measured shear (31 cm⁻¹) and breathing (81 cm⁻¹) modes at Γ, interlayer distance 3.356 Å. Data and caveats: `docs/research/bilayer-graphene-interlayer-modes.md`. Spec: story 26.

**Blocked by:** 05

**Status:** ready-for-agent

- [ ] Test: AA bilayer with one vertical spring K⊥ per atom pair has breathing mode ω² = 2K⊥/M at Γ
- [ ] Example: shear and breathing modes at Γ match 31 and 81 cm⁻¹ after the fit; plot of low-frequency branches near Γ
- [ ] Docstring: sources, that the published α⊥ would give ≈ 89 cm⁻¹, that 3.356 Å is the graphite value
- [ ] Listed in `examples/README.md`; `pytest`, `ruff check` and all examples green
