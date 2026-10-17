"""Longitudinal normal modes of a linear m-M-m molecule (CO2-like)."""

import numpy as np

from libsolidstate.spectra.modes import normal_modes

SPRING_CONSTANT = 1.0
MASS_SMALL = 1.0
MASS_LARGE = 2.0


def main():
    k = SPRING_CONSTANT
    K = k * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    masses = np.array([MASS_SMALL, MASS_LARGE, MASS_SMALL])

    omega_squared, vectors = normal_modes(K, masses)
    displacements = vectors / np.sqrt(masses)[:, None]

    for j in range(len(masses)):
        omega = np.sqrt(max(omega_squared[j], 0.0))
        print(f"omega = {omega:.4f}   displacements = {np.round(displacements[:, j], 4)}")


if __name__ == "__main__":
    main()
