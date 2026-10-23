# 01: animate() centres, scales and draws bonds

**What to build:** `animate` places any system in the middle of its window at a scale that fits, unless the caller gives an explicit scale, and can draw bonds as lines between the moving atoms. The centring and scaling live in one pure function `screen_coordinates(positions, size, scale=None)`, so the animation examples no longer compute centre and scale by hand. Spec: `../spec.md` (stories 1–6).

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent

- [x] `screen_coordinates` tested: result centred in the window, all points inside, aspect ratio kept
- [x] `animate` uses it when no scale is given; an explicit scale still works
- [x] `animate` takes optional bonds `(i, j, n, K)` and draws the n = 0 ones between current positions
- [x] All animation examples drop their own centring/scaling; 2D examples pass their bonds
- [x] `TODO.md` item for `animate()` ticked; `pytest`, `ruff check` and all examples green

