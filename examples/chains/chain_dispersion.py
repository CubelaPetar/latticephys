"""Phonon dispersion of the monatomic and the diatomic chain in the first Brillouin zone."""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

SPRING_CONSTANT = 10.0  # eV/Å²
MASS_1 = 12.011  # amu (C)
MASS_2 = 15.999  # amu (O)
LATTICE_CONSTANT = 2.5  # Å


def main():
    a = LATTICE_CONSTANT
    K = SPRING_CONSTANT
    qs = np.linspace(-np.pi / a, np.pi / a, 401)

    monatomic_chain = Crystal([[a]], [0.0], [MASS_1], ["C"])
    diatomic_chain = Crystal([[a]], [0.0, a / 2], [MASS_1, MASS_2], ["C", "O"])

    # bonds (mu, nu, n, K): atom mu in cell 0 -- atom nu in cell n
    monatomic = dispersion(monatomic_chain, [(0, 0, 1, K)], qs) * OMEGA_TO_THZ
    diatomic = dispersion(diatomic_chain, [(0, 1, 0, K), (1, 0, 1, K)], qs) * OMEGA_TO_THZ

    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    left.plot(qs * a, monatomic)
    left.set_title(f"Monatomic chain, m = {MASS_1} amu")
    right.plot(qs * a, diatomic[:, 0], label="acoustic")
    right.plot(qs * a, diatomic[:, 1], label="optical")
    right.set_title(f"Diatomic chain, m1 = {MASS_1}, m2 = {MASS_2} amu")
    right.legend()

    left.set_ylabel(r"$\nu$ (THz)")
    for axis in (left, right):
        axis.set_xlabel(r"$q a$")
        axis.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
        axis.set_xlim(-np.pi, np.pi)
        axis.grid(alpha=0.3)

    fig.suptitle(f"K = {K} eV/Å², a = {a} Å")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
