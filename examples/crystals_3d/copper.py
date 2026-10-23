"""Phonons of fcc copper with nearest-neighbour central springs: dispersion, DOS and C_V(T).

One atom per primitive cell, 12 nearest neighbours at a / sqrt(2), so 3 acoustic branches.
At X = 2 pi/a (1, 0, 0) the model gives omega_L^2 = 8 K / m and omega_T^2 = 4 K / m, so
nu_L / nu_T = sqrt(2) there. K is fitted to the measured longitudinal frequency at X,
K = m omega_L(X)^2 / 8. Then TA at X comes out as measured, LA at L within 1 %, but TA at L
about 6 % too high: real Cu has weaker forces out to further neighbours. Along Γ-X and Γ-L
the two transverse branches are degenerate by symmetry; along Γ-K they split.

Measured values (red points): room-temperature neutron scattering, most likely Svensson,
Brockhouse and Rowe, Phys. Rev. 155, 619 (1967), as quoted by Savrasov and Savrasov,
Phys. Rev. B 54, 16487 (1996), Table I; secondary source, see docs/research/cu-phonon-frequencies.md.

The DOS has van Hove peaks from flat parts of the branches near the zone boundary. The heat
capacity grows as T^3 at low T (3D Debye law) and approaches 3 k_B per atom (Dulong-Petit).
The finite q-grid has no modes below its smallest q != 0 (Γ is left out; lowest mode
about 0.5 THz for GRID = 24), which cuts C off at the lowest temperatures; near 8 K the
slope of log C against log T is still somewhat above 3.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.plot import plot_bands, plot_crystal, plot_points
from latticephys.lattice.reciprocal import k_path, q_grid, special_points
from latticephys.phonons.dispersion import dispersion
from latticephys.phonons.spectra import density_of_states
from latticephys.phonons.thermal import heat_capacity
from latticephys.units import OMEGA_TO_THZ

A = 3.615  # Å, lattice constant of Cu
MASS = 63.546  # amu
LA_X = 7.25  # THz, measured LA at X, sets K
MEASURED = {"X": [5.13, LA_X], "L": [3.42, 7.30]}  # THz, measured TA and LA at room temperature
LABELS = ["Γ", "X", "W", "K", "Γ", "L"]
GRID = 24  # GRID^3 wave vectors in the Brillouin zone
SMEARING = 0.1  # THz, Gaussian width of each mode in the DOS
TEMPERATURES = np.geomspace(8, 1000, 200)  # K
T3_GUIDE = (8.0, 30.0)  # K, range of the T^3 guide line, through the first temperature


def main():
    copper = Crystal(primitive_vectors("fcc", A), [[0.0, 0.0, 0.0]], [MASS], ["Cu"])
    K = MASS * (LA_X / OMEGA_TO_THZ) ** 2 / 8
    bonds = neighbour_bonds(copper, [K])
    print(f"K = {K:.3f} eV/Å²")

    points = special_points("fcc", copper.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])
    frequencies = dispersion(copper, bonds, qs) * OMEGA_TO_THZ
    for label, measured in MEASURED.items():
        model = dispersion(copper, bonds, [points[label]])[0] * OMEGA_TO_THZ
        print(f"{label}: model {np.unique(model.round(3)).round(2).tolist()} THz, measured {measured} THz")

    omega = dispersion(copper, bonds, q_grid(copper.vectors, GRID)[1:])  # without the zero modes at Γ
    per_atom = copper.dim * heat_capacity(omega, TEMPERATURES)

    fig, (left, middle, right) = plt.subplots(1, 3, figsize=(16, 4.5), width_ratios=[2, 1, 1.3])
    plot_bands(distance, frequencies, corners, LABELS, ax=left)
    plot_points(left, corners, LABELS, MEASURED)
    left.legend(loc="lower right")
    left.set_ylabel(r"$\nu$ (THz)")
    left.set_title(f"fcc Cu, nearest neighbours, K = {K:.2f} eV/Å²")

    grid = np.linspace(0, 1.05 * frequencies.max(), 400)
    middle.plot(density_of_states(omega * OMEGA_TO_THZ, grid, sigma=SMEARING), grid)
    middle.set_xlabel("DOS (per THz)")
    middle.set_ylim(left.get_ylim())
    middle.set_title(f"DOS, {GRID}³ q-points")

    right.loglog(TEMPERATURES, per_atom, label="Cu, spring model")
    guide = np.array(T3_GUIDE)
    right.loglog(guide, per_atom[0] * (guide / TEMPERATURES[0]) ** 3, "k:", label=r"$\propto T^3$")
    right.axhline(copper.dim, color="grey", ls="--", lw=0.8, label="Dulong-Petit")
    right.set_xlabel("T (K)")
    right.set_ylabel(r"$C_V$ per atom ($k_B$)")
    right.legend()

    for axis in (middle, right):
        axis.grid(alpha=0.3)
    fig.tight_layout()

    ax = plot_crystal(copper, radii=[0.4], colors=["tab:brown"], cells=2)
    ax.set_title("fcc Cu (atoms not to scale) and its primitive cell")
    plt.show()


if __name__ == "__main__":
    main()
