# 03: Qualitative claims: test or reword

**What to build:** for each claim, a test if a cheap independent check exists; otherwise the docstring says what is plotted. Spec: `../spec.md`.

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [ ] `chains/wave_packet.py:4` packet moves at v_g, slows near the zone edge, spreads (check: v_g at the packet's q)
- [ ] `lattices_2d/wave_2d.py:3` longitudinal front ahead of the transverse one (check: c_L > c_T)
- [ ] `lattices_2d/dos_2d.py:3` log van Hove peaks at M, steps at band edges
- [ ] `crystals_3d/silicon.py` "optical branches rise" away from Γ
- [ ] `pytest`, `ruff check`, all examples headless
