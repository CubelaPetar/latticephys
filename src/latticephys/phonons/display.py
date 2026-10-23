"""Display a precomputed 1D or 2D trajectory with pygame."""

import numpy as np
import pygame
from matplotlib import colormaps

BLACK = (0, 0, 0)
GREY = (90, 90, 90)
MARGIN = 0.1  # fraction of the window kept free on each side


def displacement_colors(u, color_scale):
    """RGB colours (steps, N, 3) of the viridis colormap by |u_i| / color_scale, saturating at 1.

    u: (steps, N, d) displacements per atom
    """
    level = np.clip(np.linalg.norm(u, axis=2) / color_scale, 0, 1)
    return (colormaps["viridis"](level)[..., :3] * 255).astype(np.uint8)


def screen_coordinates(positions, size, scale=None):
    """Pixel coordinates (..., 2) of positions (..., 2) in Å, centred in a window of size (width, height).

    The middle of the bounding box of all points goes to the window centre; y points up
    in Å and down on screen. scale is pixels per Å; if None, the largest scale (equal in
    x and y) that keeps all points inside the window with a margin of MARGIN per side.
    """
    r = np.asarray(positions, dtype=float)
    points = r.reshape(-1, 2)
    low, high = points.min(axis=0), points.max(axis=0)
    if scale is None:
        fits = [(1 - 2 * MARGIN) * s / e for s, e in zip(size, high - low) if e > 0]
        scale = min(fits, default=1.0)
    pixels = (r - (low + high) / 2) * scale * [1, -1]
    return pixels + np.array(size) / 2


def animate(equilibrium_positions, u, dt, speed, radii, colors, scale=None, fps=60, size=(800, 800), color_scale=None,
            bonds=()):
    """Play the trajectory in a loop until the window is closed.

    equilibrium_positions: (N,) or (N, 1) positions in 1D, (N, 2) in 2D, in Å
    u:                     (steps, N) or (steps, 2 N) displacements in Å, e.g. from
                           simulate(), ordered atom by atom [u_0x, u_0y, u_1x, ...]
    dt:                    time step of u in internal time units
    speed:                 simulation time shown per second of wall-clock time
    scale:                 pixels per Å; if None, the system fits the window
                           (see screen_coordinates); it is centred either way
    color_scale:           if given, colour each atom by its displacement |u_i| instead of
                           colors (dark at rest, yellow at |u_i| >= color_scale in Å), so
                           small waves become visible
    bonds:                 bonds (i, j, n, K) of the system, e.g. from periodic_supercell();
                           those with n = 0 are drawn as lines between the current
                           positions of atoms i and j, those to other cells are skipped
                           (they would cross the window)
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
    lines = [(i, j) for i, j, n, _ in bonds if not np.any(n)]
    if any(not (0 <= i < len(x0) and 0 <= j < len(x0)) for i, j in lines):
        raise ValueError(f"bond atom index outside 0 .. {len(x0) - 1}")
    pixels = screen_coordinates(x0 + u, size, scale)

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
        for i, j in lines:
            pygame.draw.line(window, GREY, pixels[frame, i], pixels[frame, j], 2)
        for centre, radius, color in zip(pixels[frame], radii, frame_colors[frame]):
            pygame.draw.circle(window, color, centre, radius)
        pygame.display.update()

        time += speed / fps

    pygame.quit()
