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

## Quickstart

```python
import numpy as np
from libsolidstate.spectra.dispersion import dispersion
from libsolidstate.units import OMEGA_TO_THZ

a = 2.5                                        # Å
masses = [12.011, 15.999]                      # amu
springs = [(0, 1, 0, 10.0), (1, 0, 1, 10.0)]   # (atom, atom, cell offset, k in eV/Å²)
qs = np.linspace(0, np.pi / a, 5)
print(dispersion(qs, masses, springs, a) * OMEGA_TO_THZ)  # acoustic and optical branch in THz
```

## Examples

Run from the repo root, e.g. `.venv/bin/python examples/wave_packet.py`. See [`examples/README.md`](examples/README.md) for the list.

## Layout

| Path | Content |
|---|---|
| `src/libsolidstate/lattice/` | Bravais lattices, crystal structures, reciprocal lattice, Brillouin zone, force constants |
| `src/libsolidstate/spectra/` | Normal modes, dynamical matrix, dispersions |
| `src/libsolidstate/dynamics/` | Time integration and pygame animation |
| `src/libsolidstate/units.py` | Unit system and conversions |
| `examples/` | Concrete systems (CO₂, chains, 2D crystals) |
| `tests/` | Checks against analytic results |
