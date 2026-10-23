import numpy as np
import pytest

from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.plot import plot_crystal


@pytest.mark.parametrize(
    "crystal",
    [
        embed(Crystal(np.eye(2), [[0.0, 0.0]], [1.0], ["A"]), 3),  # layer, d = 3
    ],
)
def test_plot_crystal_raises_for_what_it_cannot_draw(crystal):
    """plot_crystal draws dimension 2 only; a layer (d = 3) raises instead of giving a wrong picture."""
    with pytest.raises(ValueError, match="dimension 2"):
        plot_crystal(crystal, radii=[0.3] * len(crystal.masses), colors=["red"] * len(crystal.masses))


def test_plot_crystal_draws_a_molecule():
    """p = 0: one circle per atom and no cell outline."""
    molecule = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.2, 0.0], [0.6, 1.0]], [1.0] * 3, ["A"] * 3)
    ax = plot_crystal(molecule, radii=[0.3] * 3, colors=["red"] * 3)
    assert len(ax.patches) == 3
