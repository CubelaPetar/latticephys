import numpy as np
import pytest

from libsolidstate.spectra.modes import normal_modes


@pytest.mark.parametrize("k, m, M", [(1.0, 1.0, 2.0), (10.0, 16.0, 12.0)])
def test_co2_frequencies(k, m, M):
    """Linear m-M-m molecule: omega^2 = 0, k/m, k/m (1 + 2m/M)."""
    K = k * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]])
    omega_squared, _ = normal_modes(K, [m, M, m])

    np.testing.assert_allclose(omega_squared, [0, k / m, k / m * (1 + 2 * m / M)], atol=1e-12)
