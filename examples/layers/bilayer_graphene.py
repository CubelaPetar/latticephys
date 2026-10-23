"""AB (Bernal) bilayer graphene: two stiff sheets coupled by weak interlayer springs.

The bilayer is one layer crystal (p = 2, d = 3) with 4 atoms per cell. Atom A2 of the
upper sheet sits directly above B1 of the lower one; one spring per stacked pair, with
K_perp along z and K_par in the plane, couples the sheets. At Γ two rigid-layer modes
appear at very low frequency: the shear mode (sheets slide against each other in the
plane, C mode) and the layer breathing mode (they move against each other along z).
For rigid sheets omega^2 = K / M per cell (one spring, two sheets of 2 M each); the
springs act on one sublattice only, so the sheets bend slightly and the example fits
K_perp and K_par numerically to the Raman values

    shear 31 cm^-1 (Tan et al., Nat. Mater. 11, 294 (2012)),
    breathing 81 cm^-1 (Lui et al., Nano Lett. 14, 4615 (2014)).

For comparison it converts the published force constants per area alpha into springs
per cell, K = alpha A_cell: alpha_par = 12.8e18 N/m^3 (Tan 2012) gives about 31 cm^-1,
but alpha_perp = 106e18 N/m^3 (Wu et al., ACS Nano 9, 7440 (2015), fitted to thicker
twisted samples) gives about 89 cm^-1, not the measured 81 (88.7 for rigid sheets, a bit
less here as the sheets bend). The interlayer distance 3.356 Å is half the c axis of
graphite (Trucano and Chen, Nature 258, 136 (1975)); no measurement for the bilayer was
found. In this model it only sets the geometry, the springs are fitted. In-plane springs as in graphene_layer.py.
Details: docs/research/bilayer-graphene-interlayer-modes.md.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.crystal import embed
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import bilayer, honeycomb
from latticephys.phonons.dispersion import dispersion
from latticephys.units import ANGSTROM, EV, OMEGA_TO_INV_CM

A = 2.46  # Å
INTERLAYER = 3.356  # Å, interlayer distance of graphite
MASS = 12.011  # amu
IN_PLANE = [(30.0, 6.0), (5.0, 0.5)]  # eV/Å², (K_L, K_T) for first and second neighbours, as graphene_layer.py
SHEAR = 31.0  # cm^-1, measured
BREATHING = 81.0  # cm^-1, measured
ALPHA_PAR = 12.8e18  # N/m^3, published shear force constant per area
ALPHA_PERP = 106e18  # N/m^3, published breathing force constant per area
Q_MAX = 0.15  # fraction of Γ-M shown


def model(K_perp, K_par):
    """AB bilayer and its bonds. Sheet 1: atoms 0 (A1), 1 (B1) at z = 0; sheet 2 shifted by B1 - A1, so
    A2 (atom 2) sits above B1, coupled by one spring with K_perp along z and K_par in the plane."""
    sheet = embed(honeycomb(A, [MASS, MASS], ["C", "C"]), 3)
    crystal, bonds = bilayer(sheet, neighbour_bonds(sheet, IN_PLANE), sheet.positions[1] + [0.0, 0.0, INTERLAYER])
    return crystal, bonds + [(1, 2, (0, 0), (K_perp, K_par))]


def rigid_layer_modes(K_perp, K_par):
    """Shear and breathing frequency at Γ in cm^-1: the lowest nonzero modes, shear doubly degenerate."""
    omega = dispersion(*model(K_perp, K_par), [np.zeros(3)])[0] * OMEGA_TO_INV_CM
    shear, shear_2, breathing = omega[omega > 0][:3]
    if not np.isclose(shear, shear_2) or np.isclose(shear_2, breathing):
        raise ValueError("lowest modes at Γ are not a shear pair and a breathing mode, check K_perp > K_par")
    return shear, breathing


def fit():
    """K_perp, K_par giving the measured shear and breathing: start from rigid sheets (K = M omega^2), rescale."""
    K_perp, K_par = MASS * (BREATHING / OMEGA_TO_INV_CM) ** 2, MASS * (SHEAR / OMEGA_TO_INV_CM) ** 2
    for _ in range(20):
        shear, breathing = rigid_layer_modes(K_perp, K_par)
        if abs(shear - SHEAR) < 1e-6 and abs(breathing - BREATHING) < 1e-6:
            return K_perp, K_par
        K_perp, K_par = K_perp * (BREATHING / breathing) ** 2, K_par * (SHEAR / shear) ** 2
    raise RuntimeError("fit of the interlayer springs did not converge")


def main():
    K_perp, K_par = fit()
    crystal, bonds = model(K_perp, K_par)
    print(f"fitted: K_perp = {K_perp:.4f} eV/Å², K_par = {K_par:.4f} eV/Å²")
    shear, breathing = rigid_layer_modes(K_perp, K_par)
    print(f"  shear {shear:.1f} cm^-1 (measured {SHEAR}), breathing {breathing:.1f} cm^-1 (measured {BREATHING})")

    cell_area = abs(np.cross(*crystal.vectors)[2]) * ANGSTROM**2  # m^2
    to_ev_per_a2 = cell_area / (EV / ANGSTROM**2)
    shear, breathing = rigid_layer_modes(ALPHA_PERP * to_ev_per_a2, ALPHA_PAR * to_ev_per_a2)
    print(f"published alpha: shear {shear:.1f} cm^-1, breathing {breathing:.1f} cm^-1")

    points = special_points("hexagonal", crystal.vectors)
    qs, distance, _ = k_path([points["Γ"], Q_MAX * points["M"]], n=100)
    frequencies = dispersion(crystal, bonds, qs) * OMEGA_TO_INV_CM

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(distance, frequencies, color="tab:blue")
    ax.plot([0, 0], [SHEAR, BREATHING], "o", color="tab:red", label="Raman at Γ: shear, breathing")
    ax.set_xlim(0, distance[-1])
    ax.set_ylim(0, 250)
    ax.set_xlabel(r"$|q|$ along Γ-M (1/Å)")
    ax.set_ylabel(r"$\tilde\nu$ (cm$^{-1}$)")
    ax.set_title("AB bilayer graphene, low-frequency branches near Γ")
    ax.legend(loc="lower right")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
