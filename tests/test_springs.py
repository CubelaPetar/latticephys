import numpy as np

from libsolidstate.lattice.springs import force_constants
from libsolidstate.spectra.modes import normal_modes


def test_1d_matches_scalar_springs():
    K = force_constants(3, [(0, 1, 2.0, 1.0), (1, 2, 2.0, -1.0)])
    np.testing.assert_array_equal(K, 2.0 * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]]))


def test_spring_only_resists_stretching():
    """A spring along x gives no restoring force for displacements along y."""
    K = force_constants(2, [(0, 1, 1.0, [3.0, 0.0])], dim=2)
    np.testing.assert_array_equal(K @ [0, 1, 0, 0], 0)
    np.testing.assert_array_equal(K @ [1, 0, 0, 0], [1, 0, -1, 0])


def test_triangle_has_three_zero_modes():
    """A free triangle in 2D: 2 translations + 1 rotation at omega = 0, 3 vibrations."""
    r = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    bonds = [(i, j, 1.0, r[j] - r[i]) for i, j in ((0, 1), (1, 2), (2, 0))]
    omega_squared, _ = normal_modes(force_constants(3, bonds, dim=2), np.ones(6))

    np.testing.assert_allclose(omega_squared[:3], 0, atol=1e-12)
    assert np.all(omega_squared[3:] > 0.1)
