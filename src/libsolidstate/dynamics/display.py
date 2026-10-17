"""Display a precomputed 1D trajectory with pygame."""

import numpy as np
import pygame

BLACK = (0, 0, 0)


def animate(equilibrium_positions, u, radii, colors, scale=100, fps=60, size=(800, 800)):
    """Show one time step per frame and loop until the window is closed.

    equilibrium_positions: (N,) positions in model units
    u:                     (steps, N) displacements, e.g. from simulate()
    scale:                 pixels per model length unit
    """
    pygame.init()
    window = pygame.display.set_mode(size)
    pygame.display.set_caption("libsolidstate")
    clock = pygame.time.Clock()

    x0 = np.asarray(equilibrium_positions, dtype=float)
    frame = 0
    running = True
    while running:
        clock.tick(fps)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        window.fill(BLACK)
        x = (x0 + u[frame]) * scale + size[0] / 2
        for xi, radius, color in zip(x, radii, colors):
            pygame.draw.circle(window, color, (xi, size[1] / 2), radius)
        pygame.display.update()

        frame = (frame + 1) % len(u)

    pygame.quit()
