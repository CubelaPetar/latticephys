import numpy as np

from latticephys.lattice.crystal import Crystal
from latticephys.phonons.modes import normal_modes
from latticephys.phonons.springs import force_constants

UNUSED = 100.0  # lattice constant of a 2D crystal used as a molecule, all bonds in cell (0, 0)


def test_1d_molecule():
    molecule = Crystal(np.zeros((0, 1)), [0.0, 1.0, 2.0], [1.0] * 3, ["A"] * 3)
    Phi = force_constants(molecule, [(0, 1, (), 2.0), (1, 2, (), 2.0)])
    np.testing.assert_array_equal(Phi, 2.0 * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]]))


def test_spring_only_resists_stretching():
    """A spring along x gives no restoring force for displacements along y."""
    pair = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [3.0, 0.0]], [1.0] * 2, ["A"] * 2)
    Phi = force_constants(pair, [(0, 1, (), 1.0)])
    np.testing.assert_allclose(Phi @ [0, 1, 0, 0], 0)
    np.testing.assert_allclose(Phi @ [1, 0, 0, 0], [1, 0, -1, 0])


def test_triangle_has_three_zero_modes():
    """A free triangle in 2D: 2 translations + 1 rotation at omega = 0, 3 vibrations."""
    triangle = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]], [1.0] * 3, ["A"] * 3)
    bonds = [(i, j, (), 1.0) for i, j in ((0, 1), (1, 2), (2, 0))]
    omega, _ = normal_modes(force_constants(triangle, bonds), np.ones(6))

    np.testing.assert_array_equal(omega[:3], 0)
    assert np.all(omega[3:] > 0.1)


def test_tetrahedron_has_six_zero_modes():
    """A free tetrahedron in 3D: 3 translations + 3 rotations, 6 vibrations."""
    corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    tetrahedron = Crystal(np.zeros((0, 3)), corners, [1.0] * 4, ["A"] * 4)
    bonds = [(i, j, (), 1.0) for i in range(4) for j in range(i + 1, 4)]
    omega, _ = normal_modes(force_constants(tetrahedron, bonds), np.ones(12))

    np.testing.assert_array_equal(omega[:6], 0)
    assert np.all(omega[6:] > 0.1)


def test_central_spring_equals_pair_without_transverse_part():
    pair = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.0, 2.0]], [1.0] * 2, ["A"] * 2)
    np.testing.assert_array_equal(
        force_constants(pair, [(0, 1, (), 3.0)]), force_constants(pair, [(0, 1, (), (3.0, 0.0))])
    )


def test_transverse_spring_breaks_rotational_invariance():
    """With K_T a free pair keeps only the 2 translations at omega = 0: its rotation costs energy."""
    pair = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.0, 0.0]], [1.0] * 2, ["A"] * 2)
    omega, _ = normal_modes(force_constants(pair, [(0, 1, (), (2.0, 0.5))]), np.ones(4))
    np.testing.assert_array_equal(omega[:2], 0)
    np.testing.assert_allclose(omega[2:] ** 2, [2 * 0.5, 2 * 2.0])


def test_molecule_equals_crystal_with_bonds_in_home_cell():
    """A bent 3-atom molecule as p = 0 (vectors (0, 2), cell index ()) has the same force
    constants and normal modes as a 2D crystal whose bonds all stay in cell (0, 0)."""
    positions = [[0.0, 0.0], [1.0, 0.2], [1.8, 1.1]]
    springs = [(0, 1, 3.0), (1, 2, (2.0, 0.4)), (0, 2, 1.5)]
    molecule = Crystal(np.zeros((0, 2)), positions, [1.0, 2.0, 3.0], ["A", "B", "C"])
    home_cell = Crystal(UNUSED * np.eye(2), positions, [1.0, 2.0, 3.0], ["A", "B", "C"])
    assert (molecule.dim, molecule.periodicity) == (2, 0)

    Phi = force_constants(molecule, [(i, j, (), K) for i, j, K in springs])
    Phi_home = force_constants(home_cell, [(i, j, (0, 0), K) for i, j, K in springs])
    np.testing.assert_array_equal(Phi, Phi_home)

    masses = np.repeat(molecule.masses, 2)
    np.testing.assert_array_equal(normal_modes(Phi, masses)[0], normal_modes(Phi_home, masses)[0])
