"""Square lattice with nearest neighbours: central springs vs. springs with a transverse part (K_L, K_T).

Central springs only resist stretching, so a row of atoms sliding along itself costs no
energy: the transverse branch vanishes along Gamma-X. A transverse stiffness K_T (the
Born-von Karman model) gives shear stiffness without second neighbours. Along Gamma-X
then omega_T = 2 sqrt(K_T / m) sin(q a / 2). Transverse springs are not rotationally
invariant, so they are a model for bulk crystals, not for free molecules or ribbons.
"""

import matplotlib.pyplot as plt

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_bands
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A = 2.5  # Å
MASS = 12.011  # amu
K_L = 10.0  # eV/Å², longitudinal
K_T = 2.0  # eV/Å², transverse
LABELS = ["Γ", "X", "M", "Γ"]


def main():
    crystal = Crystal(primitive_vectors("square", A), [[0.0, 0.0]], [MASS], ["C"])
    points = special_points("square", crystal.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for axis, K in zip(axes, (K_L, (K_L, K_T))):
        frequencies = dispersion(crystal, neighbour_bonds(crystal, [K]), qs) * OMEGA_TO_THZ
        plot_bands(distance, frequencies, corners, LABELS, ax=axis)
        axis.set_title(f"K = {K} eV/Å²")
    axes[0].set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
