import numpy as np
import pytest

from libsolidstate.lattice.chain import chain_force_constants
from libsolidstate.spectra.modes import normal_modes


@pytest.mark.parametrize("k, m, M", [(1.0, 1.0, 2.0), (10.0, 16.0, 12.0)])
def test_co2_frequencies(k, m, M):
    """Linear m-M-m molecule: omega^2 = 0, k/m, k/m (1 + 2m/M)."""
    K = k * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    omega_squared, _ = normal_modes(K, [m, M, m])

    np.testing.assert_allclose(omega_squared, [0, k / m, k / m * (1 + 2 * m / M)], atol=1e-12)


def test_light_impurity_localized_mode():
    """Montroll-Potts: omega_L = omega_max / sqrt(1 - eps^2) for m' = (1 - eps) m."""
    k, m, eps, n = 1.0, 1.0, 0.5, 41
    masses = np.full(n, m)
    masses[n // 2] = (1 - eps) * m
    omega_squared, _ = normal_modes(chain_force_constants(n, k, boundary="fixed"), masses)

    assert np.sqrt(omega_squared[-1]) == pytest.approx(2 * np.sqrt(k / m) / np.sqrt(1 - eps**2), rel=1e-9)
