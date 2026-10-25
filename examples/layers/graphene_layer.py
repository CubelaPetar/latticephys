"""Phonons of a graphene layer in 3D space: in-plane and out-of-plane branches along Γ-M-K-Γ.

The layer is periodic in two directions (p = 2) but its atoms move in three (d = 3), so
there are 3 B = 6 branches. Their labels: polarization Z (out of plane), T (in plane,
transverse) or L (in plane, longitudinal), plus A (acoustic, the two atoms move together,
omega -> 0 at Γ) or O (optical, they move against each other): ZA, TA, LA, ZO, TO, LO.
The example finds the labels from the eigenvectors at one q near Γ.

Central springs alone give no restoring force for out-of-plane motion: ZA and ZO stay at
omega = 0. A transverse stiffness K_T acts on z as well, so for in-plane bonds z decouples
from the in-plane motion and ZA, ZO become nonzero. Near Γ this ZA is linear (for nearest
neighbours only, omega = a/2 sqrt(K_T / m) q). In real graphene ZA is quadratic,
omega ~ q^2: a rigid rotation of the sheet costs no energy, which transverse springs
violate. Getting q^2 needs bending (angle) forces, not in this model.

Spring constants are illustrative; they put LO/TO at Γ near the Raman G band
(1582 cm^-1 ≈ 47.4 THz). Branches are coloured by their out-of-plane weight.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.crystal import embed
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import atom_weights
from latticephys.phonons.modes import phonon_modes
from latticephys.units import OMEGA_TO_THZ

A = 2.46  # Å
MASS = 12.011  # amu
SPRING_CONSTANTS = [(30.0, 6.0), (5.0, 0.5)]  # eV/Å², (K_L, K_T) for first and second neighbours
LABELS = ["Γ", "M", "K", "Γ"]
LABEL_AT = 0.3  # fraction of Γ-M where the branches are named


def branch_label(e, q):
    """ZA, TA, LA, ZO, TO or LO for the eigenvector e (atom by atom, 2 atoms x 3) at q != 0."""
    e = e.reshape(2, 3)
    z = np.sum(np.abs(e[:, 2]) ** 2)
    along_q = np.sum(np.abs(e[:, :2] @ q[:2] / np.linalg.norm(q)) ** 2)
    in_phase = np.sum(np.abs(e[0] + e[1]) ** 2) / 2
    return ("Z" if z > 0.5 else "L" if along_q > 0.5 else "T") + ("A" if in_phase > 0.5 else "O")


def main():
    layer = embed(honeycomb(A, [MASS, MASS], ["C", "C"]), 3)
    bonds = neighbour_bonds(layer, SPRING_CONSTANTS)
    points = special_points("hexagonal", layer.vectors)
    qs, distance, corners = k_path([points[label] for label in LABELS])

    omega, eigenvectors = phonon_modes(layer, bonds, qs)
    frequencies = omega * OMEGA_TO_THZ
    out_of_plane = atom_weights(eigenvectors, layer.dim)[..., 2].sum(axis=-1)  # z of all atoms

    fig, ax = plt.subplots(figsize=(7, 5))
    dots = ax.scatter(np.repeat(distance, frequencies.shape[1]), frequencies.ravel(), c=out_of_plane.ravel(),
                      cmap="coolwarm", vmin=0, vmax=1, s=4)
    fig.colorbar(dots, ax=ax, label="out-of-plane weight (Z)")

    i = np.searchsorted(distance, LABEL_AT * corners[1])
    for branch, nu in enumerate(frequencies[i]):
        ax.annotate(branch_label(eigenvectors[i][:, branch], qs[i]), (distance[i], nu),
                    xytext=(16 * (branch % 2), 4), textcoords="offset points", ha="center", va="bottom")

    for corner in corners:
        ax.axvline(corner, color="grey", lw=0.5)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xticks(corners, LABELS)
    ax.set_xlim(corners[0], corners[-1])
    ax.set_ylabel(r"$\nu$ (THz)")
    ax.set_title(f"Graphene layer, (K_L, K_T) = {SPRING_CONSTANTS} eV/Å²")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
