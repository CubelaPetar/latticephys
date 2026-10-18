"""Phonon density of states and heat capacity of the monatomic chain, compared with Debye and Einstein.

The 1D density of states diverges at the band edge (van Hove singularity).
C_V approaches k_B per atom at high T (Dulong-Petit); at low T it is linear in T
like the Debye model, while the Einstein model is exponentially suppressed.
"""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.lattice.crystal import Crystal
from libsolidstate.spectra.dispersion import dispersion
from libsolidstate.units import K_B, OMEGA_TO_MEV, OMEGA_TO_THZ

SPRING_CONSTANT = 10.0  # eV/Å²
MASS = 12.011  # amu
LATTICE_CONSTANT = 1.5  # Å
Q_POINTS = 20000  # sampling of the Brillouin zone
TEMPERATURES = np.linspace(1, 1500, 300)  # K


def heat_capacity(omega, temperatures):
    """C_V per mode in units of k_B, averaged over the given angular frequencies."""
    energy = omega * OMEGA_TO_MEV / 1000  # eV
    x = energy[None, :] / (K_B * temperatures[:, None])
    return np.mean(x**2 * np.exp(-x) / (1 - np.exp(-x)) ** 2, axis=1)


def main():
    a, k, m = LATTICE_CONSTANT, SPRING_CONSTANT, MASS
    # midpoints of the Brillouin zone, avoiding q = 0 where omega = 0
    qs = (np.arange(Q_POINTS) + 0.5) / Q_POINTS * 2 * np.pi / a - np.pi / a
    omega = dispersion(Crystal([[a]], [0.0], [m], ["C"]), [(0, 0, 1, k)], qs)[:, 0]

    debye = a * np.sqrt(k / m) * np.abs(qs)  # linear dispersion with the same sound velocity
    einstein = np.full_like(omega, omega.mean())

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))
    left.hist(omega * OMEGA_TO_THZ, bins=100, density=True)
    left.set_xlabel(r"$\nu$ (THz)")
    left.set_ylabel("density of states (per THz)")

    for label, frequencies in (("chain", omega), ("Debye", debye), ("Einstein", einstein)):
        right.plot(TEMPERATURES, heat_capacity(frequencies, TEMPERATURES), label=label)
    right.axhline(1, color="grey", ls="--", lw=1)
    right.set_xlabel("T (K)")
    right.set_ylabel(r"$C_V$ per atom ($k_B$)")
    right.legend()

    for axis in (left, right):
        axis.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
