# 10: 3D plot_crystal

**What to build:** `plot_crystal` draws a d = 3 crystal on a matplotlib 3D axis: atoms as spheres coloured per basis atom, primitive cell edges, same arguments as in 2D. The Cu example shows it. Spec: stories 29, 35.

**Blocked by:** 09

**Status:** ready-for-agent

- [ ] `plot_crystal` dispatches on d (2 or 3) and raises otherwise
- [ ] Cu example draws the fcc crystal
- [ ] `pytest`, `ruff check` and all examples green
