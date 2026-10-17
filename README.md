# libsolidstate

Python library for lattice systems: crystal structures, phonon spectra (band structures along symmetry paths) and time-domain simulation of lattice oscillations.

## Setup

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pytest
```

Uses `pygame-ce` (drop-in fork of pygame, `import pygame`), since pygame has no wheels for recent Python versions.

## Layout

| Path | Content |
|---|---|
| `src/libsolidstate/lattice/` | Bravais lattices, basis, crystal structures |
| `src/libsolidstate/spectra/` | Dynamical matrix, normal modes, dispersions |
| `src/libsolidstate/dynamics/` | Time integration and pygame animation |
| `examples/` | Concrete systems (CO₂, chains, 2D crystals) |
| `tests/` | Checks against analytic results |
| `notes/` | Condensed matter script (LaTeX) |
