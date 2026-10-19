"""Display a precomputed 1D or 2D trajectory with pygame."""

import numpy as np
import pygame
from matplotlib import colormaps

BLACK = (0, 0, 0)


def displacement_colors(u, color_scale):
    """RGB colours (steps, N, 3) of the viridis colormap by |u_i| / color_scale, saturating at 1.

    u: (steps, N, d) displacements per atom
    """
    level = np.clip(np.linalg.norm(u, axis=2) / color_scale, 0, 1)
    return (colormaps["viridis"](level)[..., :3] * 255).astype(np.uint8)


def animate(equilibrium_positions, u, dt, speed, radii, colors, scale=100, fps=60, size=(800, 800), color_scale=None):
    """Play the trajectory in a loop until the window is closed.

    equilibrium_positions: (N,) or (N, 1) positions in 1D, (N, 2) in 2D, in Å
    u:                     (steps, N) or (steps, 2 N) displacements in Å, e.g. from
                           simulate(), ordered atom by atom [u_0x, u_0y, u_1x, ...]
    dt:                    time step of u in internal time units
    speed:                 simulation time shown per second of wall-clock time
    scale:                 pixels per Å; the origin is the window centre
    color_scale:           if given, colour each atom by its displacement |u_i| instead of
                           colors (dark at rest, yellow at |u_i| >= color_scale in Å), so
                           small waves become visible
    """
    x0 = np.asarray(equilibrium_positions, dtype=float)
    x0 = x0.reshape(len(x0), -1)
    if x0.shape[1] == 1:
        x0 = np.column_stack([x0, np.zeros_like(x0)])
    if x0.shape[1] != 2:
        raise ValueError("animate shows 1D or 2D systems")
    u = np.asarray(u).reshape(len(u), len(x0), -1)
    if u.shape[2] == 1:
        u = np.concatenate([u, np.zeros_like(u)], axis=2)
    if color_scale is None:
        frame_colors = np.broadcast_to(colors, (len(u), len(x0), 3))
    else:
        frame_colors = displacement_colors(u, color_scale)

    pygame.init()
    window = pygame.display.set_mode(size)
    pygame.display.set_caption("latticephys")
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
        for (x, y), radius, color in zip(r, radii, frame_colors[frame]):
            pygame.draw.circle(window, color, (size[0] / 2 + x, size[1] / 2 - y), radius)
        pygame.display.update()

        time += speed / fps

    pygame.quit()
