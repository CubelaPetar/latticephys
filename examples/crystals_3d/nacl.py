"""Rocksalt NaCl with central springs: acoustic and optical branches of a two-atom 3D crystal.

Two atoms per primitive cell give 3 acoustic (TA, TA, LA) and 3 optical (TO, TO, LO)
branches. With nearest-neighbour springs only (Na-Cl along the cube axes), a shear along
Γ-X stretches no spring to first order: both transverse acoustic branches are 0 there and
along X-W (one of them also along Γ-K), as in the square lattice. Springs to the second neighbours (Na-Na, Cl-Cl
along the face diagonals) give it a finite frequency.

Missing here: in the real ionic crystal the long-range Coulomb forces split LO from TO at
Γ (LO-TO splitting); central short-range springs give LO = TO at Γ. Spring constants are
illustrative.
"""

import matplotlib.pyplot as plt

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.plot import plot_bands, plot_crystal
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import rocksalt
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A = 5.64  # Å, cubic lattice constant of NaCl
MASSES = [22.990, 35.453]  # amu, Na and Cl
SPRING_CONSTANTS = [1.5, 0.3]  # eV/Å², first (Na-Cl) and second (Na-Na, Cl-Cl) neighbours
LABELS = ["Γ", "X", "W", "K", "Γ", "L"]


def main():
    nacl = rocksalt(A, MASSES, ["Na", "Cl"])
    points = special_points("fcc", nacl.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for axis, spring_constants in zip(axes, (SPRING_CONSTANTS[:1], SPRING_CONSTANTS)):
        frequencies = dispersion(nacl, neighbour_bonds(nacl, spring_constants), qs) * OMEGA_TO_THZ
        plot_bands(distance, frequencies, corners, LABELS, ax=axis)
        axis.set_title(f"NaCl, K = {spring_constants} eV/Å²")
    axes[0].set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()

    ax = plot_crystal(nacl, radii=[0.4, 0.6], colors=["tab:purple", "tab:green"], cells=2)
    ax.set_title("NaCl (atoms not to scale) and its primitive cell")
    plt.show()


if __name__ == "__main__":
    main()
