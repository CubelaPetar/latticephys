# latticephys

Python library for lattice systems: crystal structures, phonon spectra (band structures along symmetry paths) and time-domain simulation of lattice oscillations.

## Setup

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pytest
.venv/bin/ruff check
```

Uses `pygame-ce` (drop-in fork of pygame, `import pygame`), since pygame has no wheels for recent Python versions.

## Units

Inputs in Å, amu and eV (spring constants in eV/Å²). Time is then in units of √(amu·Å²/eV) ≈ 10.18 fs. `latticephys.units` has the factors to convert ω to THz, meV and cm⁻¹.

## Quickstart

```python
import numpy as np
from latticephys.lattice.crystal import Crystal
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

a = 2.5                                                              # Å
chain = Crystal([[a]], [0.0, a / 2], [12.011, 15.999], ["C", "O"])   # vectors, positions, masses (amu), names
bonds = [(0, 1, 0, 10.0), (1, 0, 1, 10.0)]   # (atom mu, atom nu, cell of nu, K in eV/Å²)
qs = np.linspace(0, np.pi / a, 5)
print(dispersion(chain, bonds, qs) * OMEGA_TO_THZ)  # acoustic and optical branch in THz
```

## Examples

Run from the repo root, e.g. `.venv/bin/python examples/chains/wave_packet.py`. See [`examples/README.md`](examples/README.md) for the list.

## Layout

| Path | Content |
|---|---|
| `src/latticephys/lattice/` | Geometry, independent of the physics: Bravais lattices, crystal structures, reciprocal lattice, Brillouin zone, k-paths, bonds and supercells, plots |
| `src/latticephys/phonons/` | Springs (central or (K_L, K_T)), dynamical matrix, normal modes, dispersions, group and sound velocities, DOS, scattering intensity S(Q, ω), thermal properties, time integration, pygame animation |
| `src/latticephys/units.py` | Unit system and conversions |
| `src/latticephys/export.py` | JSON export, e.g. for web pages |
| `examples/` | Concrete systems by topic (`molecules/`, `chains/`, `lattices_2d/`, `export/`); planned ones in `IDEAS.md` |
| `tests/` | Checks against analytic results |
