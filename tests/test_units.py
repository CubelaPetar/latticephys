import pytest

from latticephys import units


def test_time_unit():
    assert units.TIME_FS == pytest.approx(10.1805, rel=1e-5)


def test_frequency_conversions_consistent():
    """1 THz = 4.135667 meV = 33.35641 cm^-1."""
    assert units.OMEGA_TO_MEV / units.OMEGA_TO_THZ == pytest.approx(4.135667, rel=1e-6)
    assert units.OMEGA_TO_INV_CM / units.OMEGA_TO_THZ == pytest.approx(33.35641, rel=1e-6)
