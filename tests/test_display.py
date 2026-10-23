import numpy as np
from matplotlib import colormaps

from latticephys.phonons.display import displacement_colors, screen_coordinates


def test_displacement_colors_scale_and_saturate():
    """|u| = 0 gives the start of viridis, |u| >= color_scale its end; the direction does not matter."""
    u = np.array([[[0.0, 0.0], [0.3, 0.4], [-1.0, 2.0], [0.0, -0.25]]])
    colors = displacement_colors(u, color_scale=0.5)

    expected = (colormaps["viridis"]([0.0, 1.0, 1.0, 0.5])[:, :3] * 255).astype(np.uint8)
    np.testing.assert_array_equal(colors[0], expected)


def test_screen_coordinates_centre_fit_and_keep_aspect():
    """A 4 x 1 Å rectangle off the origin in an 800 x 600 window with 10 % margins: the
    largest common scale is min(640 / 4, 480 / 1) = 160 px/Å, so it becomes 640 x 160 px
    (ratio 4 : 1 kept) centred on (400, 300), with screen y pointing down."""
    positions = np.array([[10.0, 5.0], [14.0, 5.0], [14.0, 6.0], [10.0, 6.0]])
    r = screen_coordinates(positions, (800, 600))

    np.testing.assert_allclose(r, [[80, 380], [720, 380], [720, 220], [80, 220]])


def test_screen_coordinates_explicit_scale_and_1d():
    """An explicit scale is used as given; a chain on the x axis (zero height) is fitted
    by its length alone: 3 Å at 50 px/Å is 150 px, and fitted it is 640 px."""
    chain = np.array([[0.0, 0.0], [1.0, 0.0], [3.0, 0.0]])

    np.testing.assert_allclose(screen_coordinates(chain, (800, 600), scale=50)[:, 0], [325, 375, 475])
    np.testing.assert_allclose(screen_coordinates(chain, (800, 600)), [[80, 300], [80 + 640 / 3, 300], [720, 300]])
