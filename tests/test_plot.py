import numpy as np
import pytest

from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.plot import plot_crystal


@pytest.mark.parametrize(
    "crystal",
    [
        Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.0, 0.0]], [1.0] * 2, ["A"] * 2),  # molecule, p = 0
        embed(Crystal(np.eye(2), [[0.0, 0.0]], [1.0], ["A"]), 3),  # layer, d = 3
    ],
)
def test_plot_crystal_raises_for_what_it_cannot_draw(crystal):
    """plot_crystal draws dimension 2 with periodicity 1 or 2 only; anything else raises instead of a wrong picture."""
    with pytest.raises(ValueError, match="dimension 2"):
        plot_crystal(crystal, radii=[0.3] * len(crystal.masses), colors=["red"] * len(crystal.masses))
