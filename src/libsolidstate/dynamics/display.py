"""Display a precomputed 1D or 2D trajectory with pygame."""

import numpy as np
import pygame

BLACK = (0, 0, 0)


def animate(equilibrium_positions, u, dt, speed, radii, colors, scale=100, fps=60, size=(800, 800)):
    """Play the trajectory in a loop until the window is closed.

    equilibrium_positions: (N,) positions in 1D or (N, 2) in 2D, in Å
    u:                     (steps, N) or (steps, 2 N) displacements in Å, e.g. from
                           simulate(), ordered atom by atom [u_0x, u_0y, u_1x, ...]
    dt:                    time step of u in internal time units
    speed:                 simulation time shown per second of wall-clock time
    scale:                 pixels per Å; the origin is the window centre
    """
    x0 = np.asarray(equilibrium_positions, dtype=float)
    if x0.ndim == 1:
        x0 = np.column_stack([x0, np.zeros_like(x0)])
    u = np.asarray(u).reshape(len(u), len(x0), -1)
    if u.shape[2] == 1:
        u = np.concatenate([u, np.zeros_like(u)], axis=2)

    pygame.init()
    window = pygame.display.set_mode(size)
    pygame.display.set_caption("libsolidstate")
    clock = pygame.time.Clock()

    time = 0.0
    running = True
    while running:
        clock.tick(fps)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        frame = int(time / dt) % len(u)
        window.fill(BLACK)
        r = (x0 + u[frame]) * scale
        for (x, y), radius, color in zip(r, radii, colors):
            pygame.draw.circle(window, color, (size[0] / 2 + x, size[1] / 2 - y), radius)
        pygame.display.update()

        time += speed / fps

    pygame.quit()
