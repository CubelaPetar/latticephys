"""A light impurity in a monatomic chain creates a localized mode above the phonon band.

For impurity mass m' = (1 - eps) m the infinite chain has (Montroll-Potts)
omega_L = omega_max / sqrt(1 - eps^2), and the amplitude decays as
((1 - eps) / (1 + eps))^|n| with alternating sign.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.modes import normal_modes
from latticephys.units import OMEGA_TO_THZ

N_ATOMS = 41
SPRING_CONSTANT = 10.0  # eV/Å²
MASS = 12.011  # amu (12C)
IMPURITY_MASS = 6.0  # amu, placed at the centre


def main():
    masses = np.full(N_ATOMS, MASS)
    centre = N_ATOMS // 2
    masses[centre] = IMPURITY_MASS

    K = chain_force_constants(N_ATOMS, SPRING_CONSTANT, boundary="fixed")
    omega, vectors = normal_modes(K, masses)

    omega_max = 2 * np.sqrt(SPRING_CONSTANT / MASS)
    eps = 1 - IMPURITY_MASS / MASS
    print(f"band edge:       {omega_max * OMEGA_TO_THZ:.3f} THz")
    print(f"localized mode:  {omega[-1] * OMEGA_TO_THZ:.3f} THz")
    print(f"analytic:        {omega_max / np.sqrt(1 - eps**2) * OMEGA_TO_THZ:.3f} THz")

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))
    left.plot(omega * OMEGA_TO_THZ, "o", ms=4)
    left.axhline(omega_max * OMEGA_TO_THZ, color="grey", ls="--", label="band edge of the pure chain")
    left.set_xlabel("mode index")
    left.set_ylabel(r"$\nu$ (THz)")
    left.legend()

    displacement = vectors[:, -1] / np.sqrt(masses)
    sites = np.arange(N_ATOMS) - centre
    right.plot(sites, displacement / np.abs(displacement).max(), "o-")
    right.set_xlabel("site relative to impurity")
    right.set_ylabel("displacement (normalized)")
    right.set_title("localized mode")

    for axis in (left, right):
        axis.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
