import numpy as np
import pytest

from latticephys.lattice.bonds import neighbour_bonds
from latticephys.lattice.structures import honeycomb
from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.dispersion import dispersion
from latticephys.phonons.modes import diagonalize, dynamical_matrix, normal_modes, phonon_modes


@pytest.mark.parametrize("K, m, M", [(1.0, 1.0, 2.0), (10.0, 16.0, 12.0)])
def test_co2_frequencies(K, m, M):
    """Linear m-M-m molecule: omega^2 = 0, K/m, K/m (1 + 2m/M)."""
    Phi = K * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    omega, _ = normal_modes(Phi, [m, M, m])

    np.testing.assert_allclose(omega**2, [0, K / m, K / m * (1 + 2 * m / M)], atol=1e-12)
    assert omega[0] == 0


def test_light_impurity_localized_mode():
    """Montroll-Potts: omega_L = omega_max / sqrt(1 - eps^2) for m' = (1 - eps) m."""
    K, m, eps, n = 1.0, 1.0, 0.5, 41
    masses = np.full(n, m)
    masses[n // 2] = (1 - eps) * m
    omega, _ = normal_modes(chain_force_constants(n, K, boundary="fixed"), masses)

    assert omega[-1] == pytest.approx(2 * np.sqrt(K / m) / np.sqrt(1 - eps**2), rel=1e-9)


def test_diagonalize_signs_and_zero_modes():
    """omega = sign(omega^2) sqrt(|omega^2|); round-off of a zero mode becomes exactly 0."""
    omega, vectors = diagonalize(np.diag([9.0, -4.0, 1e-15]))
    np.testing.assert_array_equal(omega, [-2.0, 0.0, 3.0])
    np.testing.assert_array_equal(np.abs(vectors), np.eye(3)[:, [1, 2, 0]])


def test_phonon_modes_are_orthonormal_eigenvectors_of_the_dynamical_matrix():
    """D(q) e_s = omega_s^2 e_s and e^H e = 1 at each q, and omega equals dispersion().

    omega_s^2 means sign(omega_s) omega_s^2, so unstable modes would be checked too."""
    graphene = honeycomb(2.46, [12.0, 14.0], ["B", "N"])
    bonds = neighbour_bonds(graphene, [10.0, (2.0, 0.5)])
    qs = np.array([[0.0, 0.0], [0.4, -0.3], [1.1, 0.7]])
    omega, vectors = phonon_modes(graphene, bonds, qs)

    assert omega.shape == (3, 4) and vectors.shape == (3, 4, 4)
    np.testing.assert_array_equal(omega, dispersion(graphene, bonds, qs))
    for q, w, e in zip(qs, omega, vectors):
        np.testing.assert_allclose(e.conj().T @ e, np.eye(4), atol=1e-12)
        np.testing.assert_allclose(dynamical_matrix(graphene, bonds, q) @ e, e * np.sign(w) * w**2, atol=1e-10)


def test_light_impurity_mode_decays_with_alternating_sign():
    """Montroll-Potts: around the impurity at n = 0 the localized mode is u_n = u_0 (-(1 - eps) / (1 + eps))^|n|."""
    K, m, eps, n = 1.0, 1.0, 0.5, 41
    masses = np.full(n, m)
    masses[n // 2] = (1 - eps) * m
    _, vectors = normal_modes(chain_force_constants(n, K, boundary="fixed"), masses)
    u, c = vectors[:, -1] / np.sqrt(masses), n // 2

    ratio = -(1 - eps) / (1 + eps)
    np.testing.assert_allclose(u[c + 1 : c + 6] / u[c : c + 5], ratio, rtol=1e-6)
    np.testing.assert_allclose(u[c - 5 : c] / u[c - 4 : c + 1], ratio, rtol=1e-6)
