# libsolidstate

Python library for lattice systems: crystal structures, phonon spectra (band structures along symmetry paths) and time-domain simulation of lattice oscillations.

## Setup

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pytest
```

Uses `pygame-ce` (drop-in fork of pygame, `import pygame`), since pygame has no wheels for recent Python versions.

## Units

Inputs in Å, amu and eV (spring constants in eV/Å²). Time is then in units of √(amu·Å²/eV) ≈ 10.18 fs. `libsolidstate.units` has the factors to convert ω to THz, meV and cm⁻¹.

## Examples

Run from the repo root with the venv's Python (or `source .venv/bin/activate` first, then `python ...`):

```sh
.venv/bin/python examples/co2_modes.py         # CO2 normal modes, printed
.venv/bin/python examples/co2_animation.py     # CO2 symmetric stretch, pygame window
.venv/bin/python examples/chain_animation.py   # periodic chain, pygame window
.venv/bin/python examples/chain_dispersion.py  # monatomic + diatomic omega(q), matplotlib
.venv/bin/python examples/crystals_2d.py       # NaCl plane, graphene, AFM, matplotlib
```

Parameters are constants at the top of each example. Animations loop until the window is closed; `SECONDS_PER_PERIOD` sets the playback speed, `TIMESTEP` only the accuracy.

## Layout

| Path | Content |
|---|---|
| `src/libsolidstate/lattice/` | Bravais lattices, basis, crystal structures |
| `src/libsolidstate/spectra/` | Dynamical matrix, normal modes, dispersions |
| `src/libsolidstate/dynamics/` | Time integration and pygame animation |
| `examples/` | Concrete systems (CO₂, chains, 2D crystals) |
| `tests/` | Checks against analytic results |
| `notes/` | Condensed matter script (LaTeX) |
