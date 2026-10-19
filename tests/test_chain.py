import numpy as np
import pytest

from libsolidstate.phonons.chain import chain_force_constants


def test_free_chain():
    K = chain_force_constants(3, 2.0, boundary="free")
    np.testing.assert_array_equal(K, 2.0 * np.array([[1, -1, 0], [-1, 2, -1], [0, -1, 1]]))


def test_fixed_chain():
    K = chain_force_constants(3, 1.0, boundary="fixed")
    np.testing.assert_array_equal(K, [[2, -1, 0], [-1, 2, -1], [0, -1, 2]])


def test_periodic_chain_dispersion():
    """Eigenvalues of K/m match omega^2 = 4k/m sin^2(q a / 2) with q a = 2 pi j / n."""
    n, k, m = 8, 1.5, 2.0
    omega_squared = np.linalg.eigvalsh(chain_force_constants(n, k, boundary="periodic") / m)

    qa = 2 * np.pi * np.arange(n) / n
    expected = np.sort(4 * k / m * np.sin(qa / 2) ** 2)
    np.testing.assert_allclose(omega_squared, expected, atol=1e-12)


def test_unknown_boundary():
    with pytest.raises(ValueError):
        chain_force_constants(3, 1.0, boundary="closed")
