"""Unit system: lengths in Å, masses in amu, energies in eV (k in eV/Å²).

With these inputs time comes out in units of sqrt(amu Å² / eV) ≈ 10.18 fs
and angular frequencies omega in the inverse of that.
"""

import math

AMU = 1.66053906892e-27  # kg (CODATA 2022)
EV = 1.602176634e-19  # J (exact)
ANGSTROM = 1e-10  # m
HBAR = 1.054571817e-34  # J s (exact)
SPEED_OF_LIGHT = 299792458.0  # m/s (exact)
K_B = 8.617333262e-5  # eV/K

TIME = math.sqrt(AMU * ANGSTROM**2 / EV)  # s
TIME_FS = TIME / 1e-15

# multiply an angular frequency omega in internal units to get ...
OMEGA_TO_THZ = 1 / (2 * math.pi * TIME) / 1e12  # ... nu = omega / 2 pi in THz
OMEGA_TO_MEV = HBAR / TIME / EV * 1e3  # ... hbar omega in meV
OMEGA_TO_INV_CM = 1 / (2 * math.pi * TIME) / SPEED_OF_LIGHT / 100  # ... nu / c in cm^-1
