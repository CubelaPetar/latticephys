# 01: Test the analytic and symmetry claims

**What to build:** one test per claim below, against its formula, symmetry or limit (docstring states it). If a claim turns out wrong, fix the prose. Spec: `../spec.md`.

**Blocked by:** None (can start immediately)

**Status:** done

- [x] `chains/mass_defect.py:4` amplitude decays as ((1-ε)/(1+ε))^|n| with alternating sign
- [x] `chains/heat_capacity.py:4` 1D: C_V linear in T at low T, Dulong–Petit k_B per atom
- [x] `crystals_3d/copper.py:4` at X: ω_L² = 8K/m, ω_T² = 4K/m (now only indirect)
- [x] `crystals_3d/copper.py:8` T branches degenerate along Γ–X and Γ–L, split along Γ–K
- [x] `crystals_3d/copper.py:15` C ∝ T³ at low T, 3 k_B limit
- [x] `crystals_3d/iron.py:5` bcc with nearest neighbours only is floppy along Γ–N
- [x] `crystals_3d/magnesium.py:12` branches pair up at A
- [x] `crystals_3d/nacl.py:5,10` nearest neighbours only: TA = 0 along Γ–X, X–W, Γ–K; LO = TO at Γ
- [x] `lattices_2d/kagome_lieb.py:6` Lieb: third branch vanishes along Γ–X
- [x] `lattices_2d/phonons_2d.py:3` square floppy along Γ–X; graphene with nearest neighbours has a zero branch
- [x] `lattices_2d/sound_velocities_2d.py:4` square anisotropic, graphene isotropic
- [x] `lattices_2d/structure_factor_2d.py:6` transverse branches invisible along Γ–M
- [x] `lattices_2d/dos_2d.py:5` graphene: C ∝ T², 2 k_B limit (now only for the triangular lattice)
- [x] `lattices_2d/mode_animation_2d.py:6` Γ optical and K modes degenerate
- [x] `lattices_2d/nanoribbon.py:5` lowest ribbon branch ∝ q_x² (now only a generic strip)
- [x] `layers/bilayer_graphene.py:8` rigid-sheet mode ω² = K/M (now only AA)
- [x] `docs/guide/models.md:12` floppy list: square, honeycomb, rocksalt, bcc
- [x] `src/latticephys/phonons/dispersion.py:21` sound velocity relative error ~ (dq a)²
- [x] `pytest`, `ruff check`
