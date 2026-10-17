"""Animate one Bloch mode of the diatomic chain: acoustic or optical branch at wave vector q.

At q = 0 the acoustic mode is a rigid translation (starting at rest nothing moves),
in the optical mode the two
atoms of each cell move against each other. At the zone edge q = pi/a only the
heavy atoms move on the acoustic branch and only the light atoms on the optical branch.
"""

import numpy as np

from libsolidstate.dynamics.display import animate
from libsolidstate.dynamics.integrate import simulate
from libsolidstate.lattice.chain import chain_force_constants
from libsolidstate.spectra.dispersion import dynamical_matrix
from libsolidstate.units import OMEGA_TO_THZ

BRANCH = "optical"  # "acoustic" or "optical"
QA = np.pi / 2  # q times a, in [-pi, pi]; must be a multiple of 2 pi / CELLS

CELLS = 8  # periodic chain of CELLS unit cells
SPRING_CONSTANT = 10.0  # eV/Å²
MASS_1 = 12.011  # amu
MASS_2 = 15.999  # amu
LATTICE_CONSTANT = 3.0  # Å, two atoms per cell
AMPLITUDE = 0.2  # Å, exaggerated for visibility

TIMESTEP = 0.02  # internal units
PERIODS = 5
SECONDS_PER_PERIOD = 2.0

RED = (188, 39, 50)
BLUE = (100, 149, 237)


def main():
    if not np.isclose(QA * CELLS / (2 * np.pi), round(QA * CELLS / (2 * np.pi))):
        raise ValueError("QA must be a multiple of 2 pi / CELLS to fit the periodic chain")
    a = LATTICE_CONSTANT
    q = QA / a
    masses = np.array([MASS_1, MASS_2])
    springs = [(0, 1, 0, SPRING_CONSTANT), (1, 0, 1, SPRING_CONSTANT)]

    omega_squared, vectors = np.linalg.eigh(dynamical_matrix(q, masses, springs, a))
    branch = {"acoustic": 0, "optical": 1}[BRANCH]
    omega = np.sqrt(max(omega_squared[branch], 0.0))
    print(f"{BRANCH} mode at qa = {QA:.3f}: nu = {omega * OMEGA_TO_THZ:.2f} THz")

    # standing Bloch wave u_{mu,n} = Re[e_mu exp(i q n a)] / sqrt(m_mu), starting at rest
    cells = np.arange(CELLS)
    polarization = vectors[:, branch] / np.sqrt(masses)
    u0 = np.real(np.outer(np.exp(1j * q * cells * a), polarization)).ravel()
    u0 *= AMPLITUDE / np.abs(u0).max()

    all_masses = np.tile(masses, CELLS)
    K = chain_force_constants(2 * CELLS, SPRING_CONSTANT, boundary="periodic")
    # rigid translation at q = 0 has no period; show it with the optical period instead
    period = 2 * np.pi / (omega if omega > 1e-8 else np.sqrt(omega_squared[1]))
    steps = round(PERIODS * period / TIMESTEP)
    u, _ = simulate(K, all_masses, u0, np.zeros_like(u0), TIMESTEP, steps)

    x0 = a / 2 * (np.arange(2 * CELLS) - (2 * CELLS - 1) / 2)
    animate(
        x0,
        u,
        TIMESTEP,
        period / SECONDS_PER_PERIOD,
        radii=[12, 15] * CELLS,
        colors=[RED, BLUE] * CELLS,
        scale=30,
    )


if __name__ == "__main__":
    main()
