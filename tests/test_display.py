import numpy as np
from matplotlib import colormaps

from libsolidstate.phonons.display import displacement_colors


def test_displacement_colors_scale_and_saturate():
    """|u| = 0 gives the start of viridis, |u| >= color_scale its end; the direction does not matter."""
    u = np.array([[[0.0, 0.0], [0.3, 0.4], [-1.0, 2.0], [0.0, -0.25]]])
    colors = displacement_colors(u, color_scale=0.5)

    expected = (colormaps["viridis"]([0.0, 1.0, 1.0, 0.5])[:, :3] * 255).astype(np.uint8)
    np.testing.assert_array_equal(colors[0], expected)
