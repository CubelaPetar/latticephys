"""Normal modes of a finite chain between fixed walls lie on the infinite-chain dispersion.

The N standing waves of the finite chain pick the discrete wave vectors
q_j = j pi / ((N + 1) a), j = 1 .. N, from the continuous curve omega(q).
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.dispersion import dispersion
from latticephys.phonons.modes import normal_modes
from latticephys.units import OMEGA_TO_THZ

N_ATOMS = 10
SPRING_CONSTANT = 10.0  # eV/Å²
MASS = 12.011  # amu
LATTICE_CONSTANT = 1.5  # Å


def main():
    a = LATTICE_CONSTANT
    qs = np.linspace(0, np.pi / a, 200)
    chain = Crystal([[a]], [0.0], [MASS], ["C"])
    curve = dispersion(chain, [(0, 0, 1, SPRING_CONSTANT)], qs)[:, 0] * OMEGA_TO_THZ

    K = chain_force_constants(N_ATOMS, SPRING_CONSTANT, boundary="fixed")
    omega, vectors = normal_modes(K, np.full(N_ATOMS, MASS))
    q_finite = np.arange(1, N_ATOMS + 1) * np.pi / ((N_ATOMS + 1) * a)

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4))
    left.plot(qs * a, curve, label="infinite chain")
    left.plot(q_finite * a, omega * OMEGA_TO_THZ, "o", label=f"fixed chain, N = {N_ATOMS}")
    left.set_xlabel(r"$q a$")
    left.set_ylabel(r"$\nu$ (THz)")
    left.set_xticks([0, np.pi / 2, np.pi], ["0", r"$\pi/2$", r"$\pi$"])
    left.legend()

    sites = np.arange(N_ATOMS + 2)
    for j in (0, 1, N_ATOMS - 1):
        profile = np.concatenate([[0], vectors[:, j], [0]])  # walls at both ends
        right.plot(sites, profile, "o-", label=f"mode {j + 1}")
    right.set_xlabel("site (walls at 0 and N + 1)")
    right.set_ylabel("displacement (normalized)")
    right.legend()

    for axis in (left, right):
        axis.grid(alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
