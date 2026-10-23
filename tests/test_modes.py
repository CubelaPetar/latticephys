import numpy as np
import pytest

from latticephys.phonons.chain import chain_force_constants
from latticephys.phonons.modes import diagonalize, normal_modes


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
