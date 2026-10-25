"""Tools on a layer (periodicity 2, dimension 3) and a chain in a plane (periodicity 1, dimension 2),
checked against the supercell, the 2D crystal or analytic formulas."""

import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds, periodic_supercell
from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.plot import plot_brillouin_zone
from latticephys.lattice.reciprocal import brillouin_zone, q_grid, special_points, wigner_seitz_cell
from latticephys.lattice.structures import bilayer, honeycomb
from latticephys.phonons.dispersion import bloch_mode, dispersion, group_velocities, sound_velocities
from latticephys.phonons.display import animate
from latticephys.phonons.modes import normal_modes
from latticephys.phonons.spectra import density_of_states, structure_factor
from latticephys.phonons.springs import force_constants
from latticephys.phonons.thermal import heat_capacity, mean_square_displacement
from latticephys.units import K_B

A, M = 2.46, 12.0
SPRINGS = [(10.0, 2.0), (2.0, 0.5)]  # (K_L, K_T), first and second neighbours: no zero modes at q != 0


def graphene():
    return honeycomb(A, [M, M], ["C", "C"])


def layer():
    return embed(graphene(), 3)


def planar_chain():
    """Diatomic chain along x moving in the plane (p = 1, d = 2)."""
    return embed(Crystal([[A]], [0.0, A / 2], [M, 2 * M], ["A", "B"]), 2)


@pytest.mark.parametrize("crystal", [layer(), planar_chain()], ids=["layer", "planar chain"])
def test_supercell_modes_are_the_dispersion_on_the_q_grid(crystal):
    """A periodic supercell of n^p cells (not n^d) has exactly the modes omega_s(q) at the n^p q of q_grid."""
    n, p, d = 4, crystal.periodicity, crystal.dim
    bonds = neighbour_bonds(crystal, SPRINGS)
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, n)
    assert len(supercell.masses) == n**p * len(crystal.masses)

    omega, _ = normal_modes(force_constants(supercell, supercell_bonds), np.repeat(supercell.masses, d))
    expected = np.sort(dispersion(crystal, bonds, q_grid(crystal.vectors, n)).ravel())
    np.testing.assert_allclose(omega, expected, atol=1e-7)


def test_bloch_mode_is_a_supercell_eigenvector():
    """The Bloch mode at M of a 2 x 2 supercell solves Phi u = omega^2 m u for every branch."""
    crystal, cells = layer(), 2
    bonds = neighbour_bonds(crystal, SPRINGS)
    supercell, supercell_bonds = periodic_supercell(crystal, bonds, cells)
    Phi = force_constants(supercell, supercell_bonds)
    masses = np.repeat(supercell.masses, 3)

    q = special_points("hexagonal", crystal.vectors)["M"]
    for branch in range(6):
        omega, U = bloch_mode(crystal, bonds, q, branch, cells)
        np.testing.assert_allclose(Phi @ U, omega**2 * masses * U, atol=1e-9)
        assert np.linalg.norm(U.real) > 0.1


def test_dos_and_heat_capacity_with_central_springs():
    """Central springs: the B = 2 z branches are 0 everywhere, the other 4 are the 2D ones, so the
    normalized DOS of the 3B = 6 branches is (4 g_2D + 2 G(omega)) / 6 and C = (4 C_2D + 2) / 6."""
    qs = q_grid(layer().vectors, 12)
    bonds = neighbour_bonds(layer(), [10.0, 2.0])
    omega = dispersion(layer(), bonds, qs)
    omega_2d = dispersion(graphene(), bonds, qs[:, :2])
    assert omega.shape == (len(qs), 6)

    grid, sigma = np.linspace(-0.5, 4.0, 300), 0.05
    zero = density_of_states(np.zeros(1), grid, sigma)
    np.testing.assert_allclose(density_of_states(omega, grid, sigma),
                               (4 * density_of_states(omega_2d, grid, sigma) + 2 * zero) / 6, atol=1e-12)

    temperatures = [10.0, 300.0, 3000.0]
    np.testing.assert_allclose(heat_capacity(omega, temperatures),
                               (4 * heat_capacity(omega_2d, temperatures) + 2) / 6, rtol=1e-12)


def test_mean_square_displacement_adds_the_out_of_plane_part():
    """Nearest-neighbour (K_L, K_T): z decouples, D_z = K_T/m [[3, -f], [-f*, 3]] with
    f = 1 + e^(i q.a1) + e^(i q.a2), so <u_z^2> = k_B T / (N K_T) (sum_(q != 0) 3 / (9 - |f|^2) + 1/12)
    (diagonal of D_z^-1; at q = 0 only ZO, omega^2 = 6 K_T / m) on top of the 2D in-plane <u^2>."""
    K_T, n, temperature = 2.0, 12, 300.0
    bonds = neighbour_bonds(layer(), [(10.0, K_T)])
    a1, a2 = layer().vectors
    qs = q_grid(layer().vectors, n)[1:]
    f2 = np.abs(1 + np.exp(1j * qs @ a1) + np.exp(1j * qs @ a2)) ** 2
    u_z2 = K_B * temperature / (n**2 * K_T) * (np.sum(3 / (9 - f2)) + 1 / 12)

    u2_2d = mean_square_displacement(graphene(), bonds, n, temperature)
    np.testing.assert_allclose(mean_square_displacement(layer(), bonds, n, temperature), u2_2d + u_z2, rtol=1e-10)


def test_structure_factor_of_a_layer():
    """In-plane Q: z modes are invisible (Q . e = 0), the rest is the 2D result. Q with a z part:
    the phonon has the in-plane part of Q, and the f-sum rule sum_s omega_s I_s = |Q|^2 sum_mu b_mu^2 / m_mu
    holds with the full |Q|^2."""
    bonds, b = neighbour_bonds(layer(), SPRINGS), np.array([0.665, 0.936])
    Qs = np.random.default_rng(5).uniform(-4, 4, size=(8, 3))

    omega_2d, intensity_2d = structure_factor(graphene(), bonds, Qs[:, :2], b)
    flat = np.column_stack([Qs[:, :2], np.zeros(len(Qs))])
    _, intensity = structure_factor(layer(), bonds, flat, b)
    np.testing.assert_allclose(np.sort(intensity, axis=1),
                               np.sort(np.column_stack([intensity_2d, np.zeros((len(Qs), 2))]), axis=1), atol=1e-12)

    omega, intensity = structure_factor(layer(), bonds, Qs, b)
    expected = (b**2 / M).sum() * (Qs**2).sum(axis=1)
    np.testing.assert_allclose((omega * intensity).sum(axis=1), expected, rtol=1e-10)


def test_structure_factor_out_of_plane_part_of_Q():
    """The phonon sees only the in-plane part of Q, so omega does not depend on Q_z. Each mode is purely
    in-plane or purely z, so I(Q_z) = I_inplane + Q_z^2 I_z exactly: I(2 Q_z) - I(0) = 4 (I(Q_z) - I(0))."""
    bonds, b = neighbour_bonds(layer(), SPRINGS), np.array([0.665, 0.936])
    Q = np.random.default_rng(7).uniform(-4, 4, size=(6, 3))
    results = [structure_factor(layer(), bonds, Q * [1, 1, scale], b) for scale in (0, 1, 2)]
    (omega_0, I_0), (omega_1, I_1), (omega_2, I_2) = results

    np.testing.assert_allclose(omega_1, omega_0, atol=1e-12)
    np.testing.assert_allclose(omega_2, omega_0, atol=1e-12)
    np.testing.assert_allclose(I_2 - I_0, 4 * (I_1 - I_0), atol=1e-9)


def test_group_velocities_with_central_springs():
    """Central springs: the 2 z branches are zero modes (velocity undefined, nan); the 4 others
    have the 2D velocities (nan where degenerate) with v_z = 0."""
    bonds = neighbour_bonds(layer(), [10.0, 2.0])
    qs = q_grid(layer().vectors, 6)[1:]
    v = group_velocities(layer(), bonds, qs)
    v_2d = group_velocities(graphene(), bonds, qs[:, :2])

    assert np.all(np.isnan(v[:, :2]))
    np.testing.assert_allclose(v[:, 2:, :2], v_2d, atol=1e-10)
    np.testing.assert_allclose(v[:, 2:, 2], np.where(np.isnan(v_2d[:, :, 0]), np.nan, 0.0), atol=1e-12)


@pytest.mark.parametrize("crystal", [layer(), planar_chain()], ids=["layer", "planar chain"])
def test_sound_velocities_only_along_the_lattice(crystal):
    """A sound velocity needs a direction in which the crystal is periodic; the last axis (z or y) raises."""
    bonds = neighbour_bonds(crystal, SPRINGS)
    along, across = np.eye(crystal.dim)[0], np.eye(crystal.dim)[-1]
    assert sound_velocities(crystal, bonds, [along]).shape == (1, crystal.dim)
    with pytest.raises(ValueError, match="span"):
        sound_velocities(crystal, bonds, [along + across])


def test_sound_velocities_raise_for_a_molecule():
    molecule = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.0, 0.0]], [M, M], ["A", "B"])
    with pytest.raises(ValueError, match="molecule"):
        sound_velocities(molecule, [(0, 1, (), 1.0)], [[1.0, 0.0]])


@pytest.mark.parametrize("vectors", [layer().vectors, [[A, 0.0]]], ids=["layer", "planar chain"])
def test_2d_zone_tools_raise_for_a_layer_or_a_ribbon(vectors):
    """Wigner-Seitz cell and Brillouin zone exist here only for p = d = 2."""
    for tool in (wigner_seitz_cell, brillouin_zone, plot_brillouin_zone):
        with pytest.raises(ValueError, match="p = d = 2"):
            tool(vectors)


def test_animate_raises_for_3d_motion():
    """animate shows d <= 2; a layer moves in 3D."""
    positions = layer().positions
    with pytest.raises(ValueError, match="1D or 2D"):
        animate(positions, np.zeros((2, positions.size)), 0.1, 1.0, radii=[5, 5], colors=[(0, 0, 0)] * 2)


def test_bilayer_stacks_two_copies():
    """bilayer: the copy is shifted, its atoms come after the originals, and each sheet keeps its own bonds."""
    sheet = layer()
    bonds = neighbour_bonds(sheet, [(10.0, 2.0)])
    crystal, both = bilayer(sheet, bonds, [0.5, 0.0, 3.0])

    np.testing.assert_array_equal(crystal.vectors, sheet.vectors)
    np.testing.assert_allclose(crystal.positions, np.vstack([sheet.positions, sheet.positions + [0.5, 0.0, 3.0]]))
    assert both == bonds + [(mu + 2, nu + 2, n, K) for mu, nu, n, K in bonds]


def test_aa_bilayer_breathing_mode_and_its_scattering():
    """AA bilayer, in-plane (K_L, K_T) springs and one vertical spring K_perp per stacked pair: the breathing mode
    (rigid sheets against each other along z, no in-plane spring stretched) has omega^2 = 2 K_perp / M; K_T
    separates it from A and B breathing on their own. At Q = (0, 0, Q_z) its eigenvector e_z = (1, 1, -1, -1) / 2
    gives I = |sum_mu Q_z e_z,mu e^(-i Q_z z_mu)|^2 / (M omega) = 4 Q_z^2 sin^2(Q_z c / 2) / (M omega),
    which vanishes at Q_z c = 2 pi: the phase between the sheets matters."""
    K_perp, c = 0.3, 3.356
    aa, bonds = bilayer(layer(), neighbour_bonds(layer(), [(10.0, 2.0)]), [0.0, 0.0, c])
    bonds += [(mu, mu + 2, (0, 0), K_perp) for mu in (0, 1)]

    breathing = np.sqrt(2 * K_perp / M)
    omega = dispersion(aa, bonds, [np.zeros(3)])[0]
    assert np.sum(np.isclose(omega, breathing, rtol=1e-12)) == 1

    Qz = np.array([0.5, 1.0, 2 * np.pi / c])
    omega, intensity = structure_factor(aa, bonds, np.column_stack([np.zeros((3, 2)), Qz]))
    branch = np.argmin(np.abs(omega[0] - breathing))
    np.testing.assert_allclose(intensity[:, branch], 4 * Qz**2 * np.sin(Qz * c / 2) ** 2 / (M * breathing), atol=1e-12)
