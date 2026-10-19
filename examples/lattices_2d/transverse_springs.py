"""Square lattice with nearest neighbours: central springs vs. springs with a transverse part (k_L, k_T).

Central springs only resist stretching, so a row of atoms sliding along itself costs no
energy: the transverse branch vanishes along Gamma-X. A transverse stiffness k_T (the
Born-von Karman model) gives shear stiffness without second neighbours. Along Gamma-X
then omega_T = 2 sqrt(k_T / m) sin(q a / 2). Transverse springs are not rotationally
invariant, so they are a model for bulk crystals, not for free molecules or ribbons.
"""

import matplotlib.pyplot as plt

from libsolidstate.lattice.bonds import neighbour_bonds
from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.crystal import Crystal
from libsolidstate.lattice.plot import plot_bands
from libsolidstate.lattice.reciprocal import k_path, special_points
from libsolidstate.phonons.dispersion import dispersion
from libsolidstate.units import OMEGA_TO_THZ

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
    for axis, k in zip(axes, (K_L, (K_L, K_T))):
        frequencies = dispersion(crystal, neighbour_bonds(crystal, [k]), qs) * OMEGA_TO_THZ
        plot_bands(distance, frequencies, corners, LABELS, ax=axis)
        axis.set_title(f"k = {k} eV/Å²")
    axes[0].set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
