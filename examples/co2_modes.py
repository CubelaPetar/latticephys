"""Longitudinal normal modes of CO2 (O=C=O) as a linear spring model.

Measured: symmetric stretch 1388 cm^-1 (Raman), asymmetric stretch 2349 cm^-1.
"""

import numpy as np

from libsolidstate.spectra.modes import normal_modes
from libsolidstate.units import OMEGA_TO_INV_CM, OMEGA_TO_THZ

SPRING_CONSTANT = 100.0  # eV/Å², C=O bond (≈ 1600 N/m)
MASS_O = 15.999  # amu
MASS_C = 12.011  # amu


def main():
    k = SPRING_CONSTANT
    K = k * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    masses = np.array([MASS_O, MASS_C, MASS_O])

    omega_squared, vectors = normal_modes(K, masses)
    displacements = vectors / np.sqrt(masses)[:, None]

    for j in range(len(masses)):
        omega = np.sqrt(max(omega_squared[j], 0.0))
        print(
            f"{omega * OMEGA_TO_INV_CM:7.1f} cm^-1  {omega * OMEGA_TO_THZ:5.1f} THz"
            f"   displacements = {np.round(displacements[:, j], 4)}"
        )


if __name__ == "__main__":
    main()
