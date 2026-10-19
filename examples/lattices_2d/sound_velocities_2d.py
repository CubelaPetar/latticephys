"""Sound velocities of 2D crystals as a function of direction (polar plot).

The triangular lattice with nearest-neighbour springs is elastically isotropic and
obeys the Cauchy relation c_L / c_T = sqrt(3). Second-neighbour springs on the square
lattice make the velocities depend on the direction; graphene (hexagonal symmetry)
stays isotropic for long waves.
"""

import matplotlib.pyplot as plt
import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais2d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.dispersion import sound_velocities
from latticephys.units import VELOCITY_TO_KM_S

MASS = 12.011  # amu
ANGLES = np.linspace(0, 2 * np.pi, 361)


def main():
    a = 2.46
    cases = [
        ("Triangular, k = [10]", Crystal(primitive_vectors("hexagonal", a), [[0.0, 0.0]], [MASS], ["C"]), [10.0]),
        ("Square, k = [10, 3]", Crystal(primitive_vectors("square", 2.5), [[0.0, 0.0]], [MASS], ["C"]), [10.0, 3.0]),
        ("Graphene, k = [25, 5]", honeycomb(a, [MASS] * 2, ["C"] * 2),
         [25.0, 5.0]),
    ]
    directions = np.column_stack([np.cos(ANGLES), np.sin(ANGLES)])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={"projection": "polar"})
    for axis, (title, crystal, spring_constants) in zip(axes, cases):
        velocities = sound_velocities(crystal, neighbour_bonds(crystal, spring_constants), directions)
        velocities = velocities * VELOCITY_TO_KM_S
        axis.plot(ANGLES, velocities[:, 1], label="longitudinal")
        axis.plot(ANGLES, velocities[:, 0], label="transverse")
        ratio = velocities[0, 1] / velocities[0, 0]
        axis.set_title(f"{title} eV/Å²\n$c_L / c_T$ along x = {ratio:.3f}")
        axis.set_rlabel_position(90)
    axes[0].legend(loc="lower left", bbox_to_anchor=(-0.1, -0.15))
    fig.suptitle("sound velocity (km/s) against direction")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
