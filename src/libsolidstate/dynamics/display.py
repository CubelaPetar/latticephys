"""Display a precomputed 1D trajectory with pygame."""

import numpy as np
import pygame

BLACK = (0, 0, 0)


def animate(equilibrium_positions, u, dt, speed, radii, colors, scale=100, fps=60, size=(800, 800)):
    """Play the trajectory in a loop until the window is closed.

    equilibrium_positions: (N,) positions in Å
    u:                     (steps, N) displacements in Å, e.g. from simulate()
    dt:                    time step of u in internal time units
    speed:                 simulation time shown per second of wall-clock time
    scale:                 pixels per Å
    """
    pygame.init()
    window = pygame.display.set_mode(size)
    pygame.display.set_caption("libsolidstate")
    clock = pygame.time.Clock()

    x0 = np.asarray(equilibrium_positions, dtype=float)
    time = 0.0
    running = True
    while running:
        clock.tick(fps)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        frame = int(time / dt) % len(u)
        window.fill(BLACK)
        x = (x0 + u[frame]) * scale + size[0] / 2
        for xi, radius, color in zip(x, radii, colors):
            pygame.draw.circle(window, color, (xi, size[1] / 2), radius)
        pygame.display.update()

        time += speed / fps

    pygame.quit()
