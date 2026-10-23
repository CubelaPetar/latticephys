"""Silicon (diamond structure) with central springs to first and second neighbours, against neutron data.

Two atoms per primitive cell on the fcc lattice give TA, TA, LA and TO, TO, LO. The two
spring constants follow in closed form from two measured frequencies:
- Γ optical: the sublattices move rigidly against each other, omega^2 = 8 K1 / (3 M);
- TA at X: nearest-neighbour central springs give no transverse force there, so only the
  fcc sublattices count, omega^2 = 4 K2 / M (as for fcc Cu).

The other measured points test the model, and it does poorly away from these two: LA = LO
at X comes out 3 % high, but the optical branches rise from Γ towards the zone boundary
while the measured ones fall (TO at X 16 % too high), and at L the branches are spread too
wide. Full dispersion measurements also show the real TA branches flat over most of the
zone. Both need forces that resist bending the tetrahedral bonds (valence-force or
bond-charge models), which short-range central springs do not have.

Measured values (points): neutron scattering at room temperature, Nilsson and Nelin, Phys.
Rev. B 6, 3777 (1972), with Dolling (1963), as tabulated by Savrasov, Phys. Rev. B 54, 16470
(1996), Table I: secondary source, see docs/research/si-phonon-frequencies.md.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.plot import plot_bands, plot_crystal, plot_points
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import diamond
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ

A = 5.431  # Å, lattice constant of Si at 300 K
MASS = 28.086  # amu
OPTICAL_GAMMA = 15.53  # THz, measured optical frequency at Γ, sets K1
TA_X = 4.49  # THz, measured TA at X, sets K2
MEASURED = {  # THz, neutron data at Γ, X, L (all branches, degenerate ones once)
    "Γ": [OPTICAL_GAMMA],
    "X": [TA_X, 12.32, 13.90],
    "L": [3.43, 11.35, 12.60, 14.68],
}
LABELS = ["Γ", "X", "W", "K", "Γ", "L"]


def main():
    K1 = 3 * MASS * (OPTICAL_GAMMA / OMEGA_TO_THZ) ** 2 / 8
    K2 = MASS * (TA_X / OMEGA_TO_THZ) ** 2 / 4
    print(f"K1 = {K1:.3f} eV/Å² (from Γ), K2 = {K2:.3f} eV/Å² (from TA at X)")

    si = diamond(A, [MASS, MASS], ["Si (A)", "Si (B)"])
    bonds = neighbour_bonds(si, [K1, K2])
    points = special_points("fcc", si.vectors)
    for label, measured in MEASURED.items():
        model = dispersion(si, bonds, [points[label]])[0] * OMEGA_TO_THZ
        print(f"{label}: model {np.unique(model.round(3)).round(2).tolist()} THz, measured {measured} THz")

    qs, distance, corners = k_path([points[label] for label in LABELS])
    frequencies = dispersion(si, bonds, qs) * OMEGA_TO_THZ

    ax = plot_bands(distance, frequencies, corners, LABELS)
    plot_points(ax, corners, LABELS, MEASURED)
    ax.legend(loc="center right")
    ax.set_ylabel(r"$\nu$ (THz)")
    ax.set_title(f"Si, central springs K1 = {K1:.2f}, K2 = {K2:.2f} eV/Å²")

    plt.tight_layout()

    ax = plot_crystal(si, radii=[0.4, 0.4], colors=["tab:blue", "tab:cyan"], cells=2)
    ax.set_title("Si: two fcc sublattices A, B shifted by a/4 (1, 1, 1); primitive cell")
    plt.show()


if __name__ == "__main__":
    main()
