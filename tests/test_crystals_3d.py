"""3D crystals (fcc, bcc, rocksalt, diamond, hcp) against closed forms and symmetries."""

import numpy as np

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import Crystal
from latticephys.lattice.reciprocal import q_grid, special_points
from latticephys.lattice.structures import diamond, hcp, rocksalt
from latticephys.phonons.dispersion import dispersion
from latticephys.units import OMEGA_TO_THZ


def two_by_two(A, B, m1, m2):
    """Frequencies of D = [[A/m1, B/sqrt(m1 m2)], [B/sqrt(m1 m2), A/m2]], ascending."""
    D = np.array([[A / m1, B / np.sqrt(m1 * m2)], [B / np.sqrt(m1 * m2), A / m2]])
    return np.sqrt(np.clip(np.linalg.eigvalsh(D), 0, None))


def test_rocksalt_along_gamma_x():
    """Rocksalt, central K1 to the 6 nearest neighbours (a/2 along the axes), K2 to the 12 second ones
    (same species, a/2 (+-1, +-1, 0)). For q = (q, 0, 0), c = cos(q a / 2), each polarization is a 2 x 2 problem
    [[A/m1, B/sqrt(m1 m2)], [B/sqrt(m1 m2), A/m2]]:
    longitudinal A = 2 K1 + 4 K2 (1 - c), B = -2 K1 c; transverse (twice) A = 2 K1 + 2 K2 (1 - c), B = -2 K1."""
    a, K1, K2, m1, m2 = 5.64, 1.5, 0.3, 22.99, 35.45
    crystal = rocksalt(a, [m1, m2], ["Na", "Cl"])
    bonds = neighbour_bonds(crystal, [K1, K2])
    assert len(bonds) == 6 + 6 + 6  # per cell: 6 Na-Cl, and 12 neighbours each for Na-Na and Cl-Cl listed once

    for t in np.linspace(0, 1, 7):
        q = 2 * np.pi / a * t
        c = np.cos(q * a / 2)
        longitudinal = two_by_two(2 * K1 + 4 * K2 * (1 - c), -2 * K1 * c, m1, m2)
        transverse = two_by_two(2 * K1 + 2 * K2 * (1 - c), -2 * K1, m1, m2)
        expected = np.sort(np.concatenate([longitudinal, transverse, transverse]))
        np.testing.assert_allclose(dispersion(crystal, bonds, [[q, 0.0, 0.0]])[0], expected, atol=1e-9)


def test_diamond_optical_mode_at_gamma_and_stability():
    """Diamond, central K1 to the 4 nearest neighbours along <111>, K2 to the 12 second ones (same sublattice).
    At Γ the sublattices move rigidly against each other: sum e e^T = 4/3 over the 4 bonds, reduced mass M/2,
    so omega^2 = 8 K1 / (3 M), triply degenerate; K2 does not enter. With K2 > 0 no mode is unstable."""
    a, K1, K2, M = 5.431, 7.0, 1.0, 28.086
    si = diamond(a, [M, M], ["Si", "Si"])
    bonds = neighbour_bonds(si, [K1, K2])

    omega = dispersion(si, bonds, [np.zeros(3)])[0]
    np.testing.assert_allclose(omega, [0, 0, 0] + [np.sqrt(8 * K1 / (3 * M))] * 3, atol=1e-12)
    assert np.all(dispersion(si, bonds, q_grid(si.vectors, 8)) >= 0)


def test_diamond_ta_at_x_depends_only_on_second_neighbours():
    """At X nearest-neighbour central springs exert no transverse force, so TA is the fcc sublattice
    value omega^2 = 4 K2 / M (as for fcc with nearest neighbours), for any K1."""
    a, K2, M = 5.431, 0.6, 28.086
    si = diamond(a, [M, M], ["Si", "Si"])
    x = special_points("fcc", si.vectors)["X"]
    for K1 in (1.0, 10.0, 30.0):
        omega = dispersion(si, neighbour_bonds(si, [K1, K2]), [x])[0]
        np.testing.assert_allclose(omega[:2], np.sqrt(4 * K2 / M), rtol=1e-12)


def test_hcp_optical_modes_at_gamma():
    """hcp, central springs K to the 6 neighbours in the other sublattice (distance d^2 = a^2/3 + c^2/4) and K2 in
    the plane. At Γ the sublattices move rigidly against each other, only the 6 interlayer springs act (reduced
    mass M/2): in plane omega^2 = 2 K a^2 / (M d^2) (twice), out of plane omega^2 = 3 K c^2 / (M d^2); K2 drops out."""
    a, c, K, K2, M = 3.209, 5.211, 1.2, 0.7, 24.305
    mg = hcp(a, c, [M, M], ["Mg", "Mg"])
    bonds = neighbour_bonds(mg, [K, K2])  # out-of-plane neighbours (3.198 Å) come before in-plane (3.209 Å)
    d2 = a**2 / 3 + c**2 / 4

    omega = dispersion(mg, bonds, [np.zeros(3)])[0]
    expected = np.sqrt([0, 0, 0, 2 * K * a**2 / (M * d2), 2 * K * a**2 / (M * d2), 3 * K * c**2 / (M * d2)])
    np.testing.assert_allclose(omega, np.sort(expected), atol=1e-12)


def fcc(K, a=3.615, m=63.546):
    crystal = Crystal(primitive_vectors("fcc", a), [[0.0, 0.0, 0.0]], [m], ["Cu"])
    return crystal, neighbour_bonds(crystal, [K])


def test_fcc_nearest_neighbours_at_x():
    """fcc, central K to the 12 neighbours a/2 (+-1, +-1, 0) etc.: D = K/m sum_r (1 - cos q.r) e e^T. At
    X = 2 pi/a (1, 0, 0) the 8 neighbours with x != 0 have q.r = +-pi, so D = diag(8, 4, 4) K / m:
    omega_L^2 = 8 K / m, omega_T^2 = 4 K / m (twice)."""
    K, m = 1.7, 63.546
    crystal, bonds = fcc(K, m=m)
    omega = dispersion(crystal, bonds, [special_points("fcc", crystal.vectors)["X"]])[0]
    np.testing.assert_allclose(omega**2, [4 * K / m, 4 * K / m, 8 * K / m], rtol=1e-12)


def test_fcc_transverse_branches_degenerate_on_symmetry_axes():
    """Along the fourfold axis Γ-X and the threefold axis Γ-L the two transverse branches are degenerate by symmetry;
    along the twofold axis Γ-K they split."""
    crystal, bonds = fcc(1.0)
    points = special_points("fcc", crystal.vectors)
    for t in (0.25, 0.5, 0.75):
        for label in ("X", "L"):
            omega = dispersion(crystal, bonds, [t * points[label]])[0]
            assert omega[1] - omega[0] < 1e-12 * omega[2]
        omega = dispersion(crystal, bonds, [t * points["K"]])[0]
        assert omega[1] - omega[0] > 0.01 * omega[2]


def test_bcc_nearest_neighbours_floppy_along_gamma_n():
    """bcc, central K1 to the 8 neighbours a/2 (+-1, +-1, +-1). For q = t (1, 1, 0) and polarization
    (1, -1, 0) / sqrt(2) every bond has q.r = 0 (signs of x and y opposite) or is perpendicular to the polarization
    (same signs), so this transverse branch is 0 on the whole line Γ-N; K2 to the 6 second neighbours stiffens it."""
    crystal = Crystal(primitive_vectors("bcc", 2.87), [[0.0, 0.0, 0.0]], [55.845], ["Fe"])
    n = special_points("bcc", crystal.vectors)["N"]
    line = np.linspace(0.1, 1, 10)[:, None] * n
    np.testing.assert_allclose(dispersion(crystal, neighbour_bonds(crystal, [1.0]), line)[:, 0], 0, atol=1e-12)
    assert np.all(dispersion(crystal, neighbour_bonds(crystal, [1.0, 0.5]), line)[:, 0] > 0.01)


def test_hcp_branches_pair_up_at_a():
    """hcp is nonsymmorphic: its screw axis and glide plane map one layer onto the other with a shift c/2. At
    A = (0, 0, pi / c), the zone face along c, these operations together with time reversal make every branch
    twofold degenerate."""
    mg = hcp(3.209, 5.211, [24.305] * 2, ["Mg"] * 2)
    omega = dispersion(mg, neighbour_bonds(mg, [1.2, 0.7]), [special_points("hexagonal", mg.vectors)["A"]])[0]
    np.testing.assert_allclose(omega[0::2], omega[1::2], rtol=1e-12)


def test_rocksalt_nearest_neighbours_decouple_into_chains():
    """Rocksalt with central K1 along the cube axes only: each Cartesian direction alpha is a diatomic chain with
    [[2 K1 / m1, -2 K1 c / sqrt(m1 m2)], [., 2 K1 / m2]], c = cos(q_alpha a / 2). Its acoustic mode vanishes when
    c = +-1: both TA along Γ-X and X-W, one along Γ-K; at Γ the three optical modes are degenerate (LO = TO)."""
    a, K1, m1, m2 = 5.64, 1.5, 22.99, 35.45
    crystal = rocksalt(a, [m1, m2], ["Na", "Cl"])
    bonds = neighbour_bonds(crystal, [K1])
    points = special_points("fcc", crystal.vectors)

    lines = [(points["Γ"], points["X"], 2), (points["X"], points["W"], 2), (points["Γ"], points["K"], 1)]
    for start, end, zeros in lines:
        for t in (0.3, 0.6):
            q = start + t * (end - start)
            chains = [two_by_two(2 * K1, -2 * K1 * np.cos(q_alpha * a / 2), m1, m2) for q_alpha in q]
            omega = dispersion(crystal, bonds, [q])[0]
            np.testing.assert_allclose(omega, np.sort(np.concatenate(chains)), atol=1e-9)
            assert np.sum(omega < 1e-9) == zeros

    optical = dispersion(crystal, bonds, [np.zeros(3)])[0][3:]
    np.testing.assert_allclose(optical, np.sqrt(2 * K1 * (1 / m1 + 1 / m2)), rtol=1e-12)


def test_silicon_central_springs_raise_the_to_branches():
    """Si with K1, K2 fitted as in examples/crystals_3d/silicon.py (8 K1 / (3 M) = omega_Γ^2 with 15.53 THz,
    4 K2 / M = omega_TA(X)^2 with 4.49 THz): the two TO branches lie higher at the zone boundary points X, L, K, W
    than at Γ, while the next branch (LO, or LA = LO at X) lies lower."""
    a, M = 5.431, 28.086
    gamma, ta_x = 15.53 / OMEGA_TO_THZ, 4.49 / OMEGA_TO_THZ
    si = diamond(a, [M, M], ["Si", "Si"])
    bonds = neighbour_bonds(si, [3 * M * gamma**2 / 8, M * ta_x**2 / 4])
    points = special_points("fcc", si.vectors)

    omega = dispersion(si, bonds, [points[label] for label in "XLKW"])
    assert np.all(omega[:, -2:] > gamma)
    assert np.all(omega[:, 3] < gamma)
