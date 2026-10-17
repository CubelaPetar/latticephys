"""Phonon dispersion of the monatomic and the diatomic chain in the first Brillouin zone."""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.spectra.dispersion import dispersion

SPRING_CONSTANT = 1.0
MASS_1 = 1.0
MASS_2 = 3.0
LATTICE_CONSTANT = 1.0


def main():
    a = LATTICE_CONSTANT
    k = SPRING_CONSTANT
    qs = np.linspace(-np.pi / a, np.pi / a, 401)

    monatomic = dispersion(qs, [MASS_1], [(0, 0, 1, k)], a)
    # atom 0 -- atom 1 inside the cell, atom 1 -- atom 0 of the next cell
    diatomic = dispersion(qs, [MASS_1, MASS_2], [(0, 1, 0, k), (1, 0, 1, k)], a)

    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    left.plot(qs * a, monatomic)
    left.set_title(f"Monatomic chain, m = {MASS_1}")
    right.plot(qs * a, diatomic[:, 0], label="acoustic")
    right.plot(qs * a, diatomic[:, 1], label="optical")
    right.set_title(f"Diatomic chain, m1 = {MASS_1}, m2 = {MASS_2}")
    right.legend()

    left.set_ylabel(r"$\omega \, / \, \sqrt{k/m_1}$")
    for axis in (left, right):
        axis.set_xlabel(r"$q a$")
        axis.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
        axis.set_xlim(-np.pi, np.pi)
        axis.grid(alpha=0.3)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
