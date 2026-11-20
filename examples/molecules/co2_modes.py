"""Longitudinal normal modes of CO2 (O=C=O) as a linear spring model.

Measured, for comparison with the printed model frequencies: symmetric stretch 1333 cm^-1
(unperturbed value), asymmetric stretch 2349 cm^-1 (Herzberg, Molecular Spectra and Molecular
Structure II, 1945, as listed in NIST CCCBDB, https://cccbdb.nist.gov/exp2x.asp?casno=124389).
"""

import numpy as np

from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.modes import normal_modes
from latticephys.units import OMEGA_TO_INV_CM, OMEGA_TO_THZ

SPRING_CONSTANT = 100.0  # eV/Å², C=O bond (≈ 1600 N/m)
MASS_O = 15.999  # amu
MASS_C = 12.011  # amu


def main():
    Phi = chain_force_constants(3, SPRING_CONSTANT, boundary="free")
    masses = np.array([MASS_O, MASS_C, MASS_O])

    omega, vectors = normal_modes(Phi, masses)
    displacements = vectors / np.sqrt(masses)[:, None]

    for j in range(len(masses)):
        print(
            f"{omega[j] * OMEGA_TO_INV_CM:7.1f} cm^-1  {omega[j] * OMEGA_TO_THZ:5.1f} THz"
            f"   displacements = {np.round(displacements[:, j], 4)}"
        )


if __name__ == "__main__":
    main()
