"""One-phonon scattering intensity S(Q, omega) of the diatomic chain over three Brillouin zones.

Neutrons or X-rays see a phonon at momentum transfer Q = q + G with intensity
|sum_mu b_mu / sqrt(m_mu) (Q e_mu) exp(-i Q r_mu)|^2 / omega. The dispersion repeats in
every zone, the intensity does not: with equal scattering lengths and atoms at 0 and a/2
the two atoms scatter in antiphase near Q = 2 pi / a, so the acoustic branch is invisible
there and the optical branch takes over. Near Q = 4 pi / a the acoustic branch is strong.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.spectra import density_of_states, structure_factor
from latticephys.units import OMEGA_TO_THZ

SPRING_CONSTANT = 10.0  # eV/Å²
MASS_1 = 12.011  # amu
MASS_2 = 15.999  # amu
SCATTERING_LENGTHS = [1.0, 1.0]  # equal, to show the extinction
LATTICE_CONSTANT = 3.0  # Å, two atoms per cell
ZONES = 3  # Q from 0 to ZONES * 2 pi / a
SMEARING = 0.4  # THz, energy resolution


def main():
    a, K = LATTICE_CONSTANT, SPRING_CONSTANT
    crystal = Crystal([[a]], [0.0, a / 2], [MASS_1, MASS_2], ["C", "O"])
    bonds = [(0, 1, 0, K), (1, 0, 1, K)]

    Qs = np.linspace(0.01, ZONES * 2 * np.pi / a, 600)
    omega, intensity = structure_factor(crystal, bonds, Qs, SCATTERING_LENGTHS)
    nu = omega * OMEGA_TO_THZ
    grid = np.linspace(0, 1.1 * nu.max(), 300)
    spectrum = np.array([density_of_states(w, grid, SMEARING, weights=i) for w, i in zip(nu, intensity)]).T

    fig, axis = plt.subplots(figsize=(10, 5))
    axis.imshow(spectrum, origin="lower", aspect="auto", cmap="magma", vmax=np.percentile(spectrum, 99),
                extent=[0, Qs[-1] * a / (2 * np.pi), 0, grid[-1]])
    axis.plot(Qs * a / (2 * np.pi), nu, color="white", lw=0.5, ls=":")
    axis.set_xlabel(r"$Q a / 2\pi$")
    axis.set_ylabel(r"$\nu$ (THz)")
    axis.set_title("S(Q, ω) of the diatomic chain; dotted: dispersion, the same in every zone")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
