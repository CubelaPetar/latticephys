"""Scattering intensity S(Q, omega) of graphene along Γ-M through three Brillouin zones.

Q runs along the reciprocal vector b1 up to 3 b1. At a reciprocal-lattice vector G = m b1 + n b2
the two atoms (at 0 and delta = (a1 + a2) / 3) scatter with the phase G . delta = 2 pi (m + n) / 3.
At G = 3 b1 they scatter in phase: the acoustic branch is strong and the optical branches, where
the atoms move against each other, vanish. Along Γ-M only longitudinal modes have Q . e != 0, so
the transverse branches (dotted) are invisible everywhere.

The extinction needs equal b / m of the two atoms; hBN (b / m of B and N differ by 8 %) is still
99.5 % extinct at 3 b1. In-plane modes only, spring constants illustrative, as in phonons_2d.py.
Scattering length: coherent neutron value of C from the NIST table (Sears 1992).
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.reciprocal import reciprocal_vectors
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.spectra import intensity_map, structure_factor
from latticephys.units import OMEGA_TO_THZ

A = 2.46  # Å
MASS = 12.011  # amu
SCATTERING_LENGTH = 6.646  # fm
SPRING_CONSTANTS = [25.0, 5.0]  # eV/Å², first and second neighbours
ZONES = 3  # Q from 0 to ZONES * b1
SMEARING = 0.5  # THz, energy resolution


def main():
    graphene = honeycomb(A, [MASS, MASS], ["C", "C"])
    b1 = reciprocal_vectors(graphene.vectors)[0]
    t = np.linspace(0.005, ZONES, 600)
    bonds = neighbour_bonds(graphene, SPRING_CONSTANTS)
    omega, intensity = structure_factor(graphene, bonds, t[:, None] * b1, [SCATTERING_LENGTH] * 2)
    nu = omega * OMEGA_TO_THZ
    grid = np.linspace(0, 1.1 * nu.max(), 300)
    spectrum = intensity_map(nu, intensity, grid, SMEARING)

    fig, axis = plt.subplots(figsize=(10, 5))
    axis.imshow(spectrum, origin="lower", aspect="auto", cmap="magma", vmax=np.percentile(spectrum, 99),
                extent=[0, ZONES, 0, grid[-1]])
    axis.plot(t, nu, color="white", lw=0.5, ls=":")
    axis.set_xlabel(r"$Q / |b_1|$ (along Γ–M)")
    axis.set_ylabel(r"$\nu$ (THz)")
    axis.set_title("S(Q, ω) of graphene; dotted: dispersion, the same in every zone")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
