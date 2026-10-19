"""A Gaussian wave packet in a monatomic chain moves with the group velocity d omega / dq.

Space-time plot of the displacements for two carrier wave vectors q0. The dashed
line is the analytic group velocity v_g = a sqrt(k/m) cos(q0 a / 2): the packet is
fast at small q0 and nearly stops at the zone edge, and it spreads due to dispersion.
"""

import matplotlib.pyplot as plt
import numpy as np

from libsolidstate.phonons.chain import chain_force_constants
from libsolidstate.phonons.integrate import simulate
from libsolidstate.units import TIME_FS

N_ATOMS = 300
SPRING_CONSTANT = 10.0  # eV/Å²
MASS = 12.011  # amu
LATTICE_CONSTANT = 1.5  # Å
CARRIERS = [0.3 * np.pi, 0.8 * np.pi]  # q0 a
WIDTH = 10  # sites, Gaussian envelope
AMPLITUDE = 0.1  # Å

TIMESTEP = 0.05  # internal units
DURATION = 150.0  # internal units (≈ 1.5 ps)


def main():
    a, k, m = LATTICE_CONSTANT, SPRING_CONSTANT, MASS
    K = chain_force_constants(N_ATOMS, k, boundary="periodic")
    masses = np.full(N_ATOMS, m)
    steps = round(DURATION / TIMESTEP)

    sites = np.arange(N_ATOMS)
    start = N_ATOMS // 4
    envelope = AMPLITUDE * np.exp(-((sites - start) ** 2) / (2 * WIDTH**2))
    length = N_ATOMS * a
    time_fs = DURATION * TIME_FS

    fig, axes = plt.subplots(1, len(CARRIERS), figsize=(12, 5), sharey=True)
    for axis, q0a in zip(axes, CARRIERS):
        omega0 = 2 * np.sqrt(k / m) * np.sin(q0a / 2)
        # right-moving wave u = envelope cos(q0 x - omega0 t) at t = 0
        u0 = envelope * np.cos(q0a * sites)
        p0 = m * omega0 * envelope * np.sin(q0a * sites)
        u, _ = simulate(K, masses, u0, p0, TIMESTEP, steps)

        axis.imshow(u, origin="lower", aspect="auto", cmap="RdBu", extent=[0, length, 0, time_fs])
        group_velocity = a * np.sqrt(k / m) * np.cos(q0a / 2)  # Å per internal time unit
        axis.plot([start * a, start * a + group_velocity * DURATION], [0, time_fs], "k--")
        axis.set_xlim(0, length)
        axis.set_title(rf"$q_0 a = {q0a / np.pi:.1f}\,\pi$")
        axis.set_xlabel("x (Å)")
    axes[0].set_ylabel("t (fs)")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
