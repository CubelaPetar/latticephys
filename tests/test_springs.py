import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.modes import normal_modes
from latticephys.phonons.springs import force_constants

UNUSED = 100.0  # lattice constant of molecules, whose bonds all stay in cell 0


def test_1d_molecule():
    molecule = Crystal([[UNUSED]], [0.0, 1.0, 2.0], [1.0] * 3, ["A"] * 3)
    K = force_constants(molecule, [(0, 1, 0, 2.0), (1, 2, 0, 2.0)])
    np.testing.assert_array_equal(K, 2.0 * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]]))


def test_spring_only_resists_stretching():
    """A spring along x gives no restoring force for displacements along y."""
    pair = Crystal(UNUSED * np.eye(2), [[0.0, 0.0], [3.0, 0.0]], [1.0] * 2, ["A"] * 2)
    K = force_constants(pair, [(0, 1, (0, 0), 1.0)])
    np.testing.assert_allclose(K @ [0, 1, 0, 0], 0)
    np.testing.assert_allclose(K @ [1, 0, 0, 0], [1, 0, -1, 0])


def test_triangle_has_three_zero_modes():
    """A free triangle in 2D: 2 translations + 1 rotation at omega = 0, 3 vibrations."""
    triangle = Crystal(UNUSED * np.eye(2), [[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]], [1.0] * 3, ["A"] * 3)
    bonds = [(i, j, (0, 0), 1.0) for i, j in ((0, 1), (1, 2), (2, 0))]
    omega, _ = normal_modes(force_constants(triangle, bonds), np.ones(6))

    np.testing.assert_array_equal(omega[:3], 0)
    assert np.all(omega[3:] > 0.1)


def test_tetrahedron_has_six_zero_modes():
    """A free tetrahedron in 3D: 3 translations + 3 rotations, 6 vibrations."""
    corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    tetrahedron = Crystal(UNUSED * np.eye(3), corners, [1.0] * 4, ["A"] * 4)
    bonds = [(i, j, (0, 0, 0), 1.0) for i in range(4) for j in range(i + 1, 4)]
    omega, _ = normal_modes(force_constants(tetrahedron, bonds), np.ones(12))

    np.testing.assert_array_equal(omega[:6], 0)
    assert np.all(omega[6:] > 0.1)


def test_central_spring_equals_pair_without_transverse_part():
    pair = Crystal(UNUSED * np.eye(2), [[0.0, 0.0], [1.0, 2.0]], [1.0] * 2, ["A"] * 2)
    np.testing.assert_array_equal(
        force_constants(pair, [(0, 1, (0, 0), 3.0)]), force_constants(pair, [(0, 1, (0, 0), (3.0, 0.0))])
    )


def test_transverse_spring_breaks_rotational_invariance():
    """With k_T a free pair keeps only the 2 translations at omega = 0: its rotation costs energy."""
    pair = Crystal(UNUSED * np.eye(2), [[0.0, 0.0], [1.0, 0.0]], [1.0] * 2, ["A"] * 2)
    omega, _ = normal_modes(force_constants(pair, [(0, 1, (0, 0), (2.0, 0.5))]), np.ones(4))
    np.testing.assert_array_equal(omega[:2], 0)
    np.testing.assert_allclose(omega[2:] ** 2, [2 * 0.5, 2 * 2.0])
