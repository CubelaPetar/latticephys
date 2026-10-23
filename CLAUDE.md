# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules Set

1. Ask, don't assume. If something is unclear, ask before writing a single line. Never make silent assumptions about intent, architecture, or requirements.
2. Simplest solution first. Always implement the simplest thing that could work. Do not add abstractions or flexibility that weren't explicitly requested.
3. Don't touch unrelated code. If a file or function is not directly part of the current task, do not modify it, even if you think it could be improved.
4. Flag uncertainty explicitly. If you are not confident about an approach or technical detail, say so before proceeding. Confidence without certainty causes more damage than admitting a gap.
5. Keep docs short. Simple, clear, concise, only the most important info. Split by topic or date rather than growing one file. Describe the current state, not the history of how we got there. When something is done, shorten its section instead of piling on notes.
6. Correct and clean before new features. The code must be correct, consistent and concise before extending it. Check every change against this list and report what fails:
   - **Physics:** each result is checked against an analytic formula, an experimental result or an independent calculation in a test.
   - **No silent errors:** inputs that would give a wrong result without notice (e.g. a wave vector that doesn't fit the supercell) raise an error.
   - **Consistent API:** new functions follow the argument order, data formats and names of existing ones. If they can't, say so before writing them.
   - **No duplication:** logic needed in two places (library, examples, tests) lives in one library function.
   - **Verify:** run `pytest`, `ruff check` and every example headless before committing.
   - Open cleanup items in `TODO.md` come before new features.

## Scope

Lattice physics with single-particle or harmonic methods: phonons now; tight-binding electrons, classical spins / linear spin waves and electron–phonon models fit later. Dimensions: 1D and 2D done, 3D next.

## Architecture

- `lattice/` is geometry only (Crystal, Bravais lattices, reciprocal lattice, k-paths, q-grids, bonds, supercells, structures, plots) and never imports a physics package.
- `phonons/` builds on `lattice/`. Later physics packages (`electrons/`, `magnons/`) sit beside it. Code moves to `lattice/` or a shared module only once a second package actually needs it; no abstract base classes in advance.
- Keep everything dimension-generic where it costs nothing; 2D-only functions check the dimension and raise.

## Conventions

- Argument order model first: `(crystal, bonds, q or qs, ...)`.
- `Crystal(vectors, positions, masses, names)`: vectors (p, d) as rows, positions (B, d) Cartesian; dimension d, periodicity p ≤ d (a layer: p = 2, d = 3; `embed` lifts a crystal into higher d). q is a Cartesian d-vector in the span of the lattice vectors. A molecule is a Crystal whose bonds all stay in cell 0.
- Bond `(mu, nu, n, K)`: atom mu in cell 0 to atom nu in cell n (p entries; an int for p = 1), each bond listed once. The last entry is the model parameter: a spring K or (K_L, K_T) now, hopping or exchange later.
- Symbols: spring constant K, force-constant matrix Φ (`Phi`), dynamical matrix D; k and q only mean wave vectors.
- Degrees of freedom atom by atom `[u_0x, u_0y, u_1x, ...]`; masses per degree of freedom `np.repeat(masses, d)`.
- Units in: Å, amu, eV (K in eV/Å²); internal time unit ≈ 10.18 fs. Convert only at output with `latticephys.units`.
- Frequencies are always ω, never ω²: negative = unstable, zero modes exactly 0. All eigenproblems of dynamical matrices go through `phonons.modes.diagonalize`.
- Tests: the docstring states the analytic formula or independent calculation being checked.
- Examples: topic folders in `examples/`, parameters as constants at the top, docstring explains the physics, `main()`. Each is listed in `examples/README.md`. Planned examples go to `examples/IDEAS.md`, the library features they need to `TODO.md`.

## Verify

`.venv/bin/pytest && .venv/bin/ruff check`. `tests/test_examples.py` runs every example headless (Agg backend, `animate()` skipped, exports into a temp dir).

Push to `origin` (Gitea); the GitHub repo `CubelaPetar/latticephys` (remote `public`) is public.

## Agent skills

### Issue tracker

Local markdown under `.scratch/<feature>/`; `TODO.md` stays the overview. See `docs/agents/issue-tracker.md`.

### Triage labels

Default five roles (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`) as a `Status:` line. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `GLOSSARY.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
