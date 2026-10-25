# 04: plot_bands coloured by weight

**What to build:** `plot_bands(distance, frequencies, corners, labels, ax=None, weights=None, weight_label=None)` draws lines as now, or with weights (same shape as frequencies) coloured dots and a colour bar labelled `weight_label`. The graphene layer example drops its copied decoration (keeps its branch annotations); the nanoribbon example uses a k-path Γ–X with those labels. Spec: `../spec.md` (stories 25–26).

**Blocked by:** 02

**Status:** ready-for-agent

- [ ] Test: with weights, one dot per frequency and a colour bar; without, unchanged
- [ ] Graphene layer and nanoribbon examples use it
- [ ] `docs/images/make_images.py` still produces the graphene layer image
- [ ] `pytest`, `ruff check`, all examples headless
