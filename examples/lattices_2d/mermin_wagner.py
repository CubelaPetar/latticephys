"""Thermal mean-square displacement against system size in 1D, 2D and 3D (Mermin-Wagner).

Classically <u^2> = k_B T / N sum over modes of 1 / (m omega^2). Long waves dominate:
in 1D <u^2> grows linearly with the length, in 2D logarithmically, in 3D it converges.
So 1D and 2D crystals have no true long-range positional order at T > 0.
The dashed line is the 2D long-wavelength result k_B T A_cell / (2 pi m) (1/c_L^2 + 1/c_T^2) ln N.
"""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.springs import neighbour_bonds
from libsolidstate.spectra.dispersion import sound_velocities
from libsolidstate.spectra.thermal import mean_square_displacement
from libsolidstate.units import K_B

MASS = 12.011  # amu
SPRING_CONSTANT = 10.0  # eV/Å², nearest neighbours
SECOND_NEIGHBOURS = 3.0  # eV/Å², needed to stabilize the simple cubic lattice
TEMPERATURE = 300.0  # K
SIZES = {1: [8, 16, 32, 64, 128, 256], 2: [8, 16, 32, 64, 128], 3: [4, 8, 12, 16, 24]}  # cells per direction


def main():
    crystals = {
        1: (Crystal([[1.5]], [0.0], [MASS], ["C"]), [SPRING_CONSTANT]),
        2: (Crystal(primitive_vectors("hexagonal", 2.46), [[0.0, 0.0]], [MASS], ["C"]), [SPRING_CONSTANT]),
        3: (Crystal(2.5 * np.eye(3), [[0.0, 0.0, 0.0]], [MASS], ["C"]), [SPRING_CONSTANT, SECOND_NEIGHBOURS]),
    }
    labels = {1: "1D chain", 2: "2D triangular", 3: "3D simple cubic"}

    fig, axis = plt.subplots(figsize=(7, 5))
    for d, (crystal, spring_constants) in crystals.items():
        bonds = neighbour_bonds(crystal, spring_constants)
        u2 = [mean_square_displacement(crystal, bonds, n, TEMPERATURE)[0] for n in SIZES[d]]
        axis.plot(SIZES[d], u2, "o-", label=labels[d])
        if d == 2:
            c_t, c_l = sound_velocities(crystal, bonds, [[1.0, 0.0]])[0]
            area = abs(np.linalg.det(crystal.vectors))
            slope = K_B * TEMPERATURE * area / (2 * np.pi * MASS) * (1 / c_l**2 + 1 / c_t**2)
            sizes = np.array(SIZES[d])
            axis.plot(sizes, u2[-1] + slope * np.log(sizes / sizes[-1]), "k--", lw=1, label="2D: slope from c_L, c_T")

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("cells per direction N")
    axis.set_ylabel(r"$\langle u^2 \rangle$ (Å$^2$)")
    axis.set_title(f"T = {TEMPERATURE:.0f} K")
    axis.legend()
    axis.grid(alpha=0.3, which="both")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
