"""Phonon bands of the kagome and the Lieb (CuO2) lattice, with and without second-neighbour springs.

With nearest-neighbour central springs only, both are floppy (omega = 0 modes):
- kagome: the lowest branch vanishes along the three Gamma-M lines of the zone,
- Lieb: two branches vanish everywhere (each O has two collinear Cu-O bonds and can
  move freely sideways), and a third vanishes along Gamma-X (Cu-O-Cu chains have no
  shear stiffness).
Second-neighbour springs remove these, except one Lieb mode at M: the CuO4 squares
rotate rigidly, alternating from cell to cell (a "rigid unit mode", like octahedral
rotations in perovskites). A third shell (Cu-Cu springs) makes it finite as well.
Masses and spring constants are illustrative.
"""

import matplotlib.pyplot as plt

from libsolidstate.lattice.bonds import neighbour_bonds
from libsolidstate.lattice.plot import plot_bands
from libsolidstate.lattice.reciprocal import k_path, special_points
from libsolidstate.lattice.structures import kagome, lieb
from libsolidstate.phonons.dispersion import dispersion
from libsolidstate.units import OMEGA_TO_THZ

A = 4.0  # Å
SPRING_CONSTANTS = [10.0, 3.0]  # eV/Å², first and second shell


def main():
    cases = [
        ("Kagome", kagome(A, [12.011] * 3, ["A", "B", "C"]), "hexagonal", ["Γ", "M", "K", "Γ"]),
        ("Lieb (CuO2)", lieb(A, [63.546, 15.999, 15.999], ["Cu", "O", "O"]), "square", ["Γ", "X", "M", "Γ"]),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey="row")
    for row, (name, crystal, lattice, labels) in zip(axes, cases):
        points = special_points(lattice, crystal.vectors)
        qs, distance, corners = k_path([points[label] for label in labels])
        for axis, spring_constants in zip(row, (SPRING_CONSTANTS[:1], SPRING_CONSTANTS)):
            frequencies = dispersion(crystal, neighbour_bonds(crystal, spring_constants), qs) * OMEGA_TO_THZ
            plot_bands(distance, frequencies, corners, labels, ax=axis)
            axis.set_title(f"{name}, k = {spring_constants} eV/Å²")
        row[0].set_ylabel(r"$\nu$ (THz)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
