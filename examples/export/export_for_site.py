"""Export dispersions and Bloch modes as JSON for a web page, e.g. the Phonon Playground.

Writes one file per system into the directory given as argument (default site_data/).
A page animates a mode as positions + displacements * cos(2 pi frequency t); the
schema is described in examples/README.md.
"""

import sys
from pathlib import Path

import numpy as np

from latticephys.export import save_json
from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.reciprocal import k_path, special_points
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import bloch_mode, dispersion
from latticephys.units import OMEGA_TO_THZ


def system_data(name, crystal, bonds, points, path, mode_points, cells):
    """Dispersion along path (labels of points) and all branches at mode_points in a periodic supercell."""
    qs, distance, ticks = k_path([points[label] for label in path])
    supercell, _ = periodic_supercell(crystal, bonds, cells)

    modes = []
    for label in mode_points:
        for branch in range(len(crystal.masses) * crystal.dim):
            omega, u0 = bloch_mode(crystal, bonds, points[label], branch, cells)
            modes.append({
                "point": label,
                "branch": branch,
                "frequency": omega * OMEGA_TO_THZ,
                "displacements": (u0 / np.abs(u0).max()).reshape(-1, crystal.dim),
            })

    return {
        "system": name,
        "units": {"length": "Å", "frequency": "THz", "mass": "amu", "displacements": "arbitrary, max 1"},
        "crystal": {"vectors": crystal.vectors, "positions": crystal.positions, "masses": crystal.masses,
                    "names": crystal.names},
        "dispersion": {"distance": distance, "frequencies": dispersion(crystal, bonds, qs) * OMEGA_TO_THZ,
                       "ticks": ticks, "labels": path},
        "supercell": {"cells": cells, "positions": supercell.positions, "names": supercell.names},
        "modes": modes,
    }


def main():
    output = Path(sys.argv[1] if len(sys.argv) > 1 else "site_data")
    output.mkdir(exist_ok=True)

    a = 3.0
    chain = Crystal([[a]], [0.0, a / 2], [12.011, 15.999], ["C", "O"])
    chain_bonds = [(0, 1, 0, 10.0), (1, 0, 1, 10.0)]
    chain_points = {"-π/a": [-np.pi / a], "Γ": [0.0], "π/a": [np.pi / a]}

    graphene = honeycomb(2.46, [12.011] * 2, ["C"] * 2)
    graphene_bonds = neighbour_bonds(graphene, [25.0, 5.0])
    graphene_points = special_points("hexagonal", graphene.vectors)

    systems = {
        "diatomic_chain": system_data(
            "diatomic chain", chain, chain_bonds, chain_points, ["-π/a", "Γ", "π/a"], ["Γ", "π/a"], 8
        ),
        "graphene": system_data(
            "graphene (in-plane)", graphene, graphene_bonds, graphene_points, ["Γ", "M", "K", "Γ"], ["Γ", "M", "K"], 6
        ),
    }
    for filename, data in systems.items():
        save_json(output / f"{filename}.json", data)
        print(f"wrote {output / filename}.json: {len(data['modes'])} modes")


if __name__ == "__main__":
    main()
