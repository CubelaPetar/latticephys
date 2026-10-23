import numpy as np
import pytest

from latticephys.lattice.bravais3d import primitive_vectors
from latticephys.lattice.crystal import Crystal, embed
from latticephys.lattice.plot import plot_crystal
from latticephys.lattice.structures import honeycomb


def test_plot_crystal_raises_for_dimension_1():
    """plot_crystal draws d = 2 or 3; a chain in 1D raises instead of giving a wrong picture."""
    chain = Crystal([[1.0]], [0.0], [1.0], ["A"])
    with pytest.raises(ValueError, match="dimension 2 or 3"):
        plot_crystal(chain, radii=[0.3], colors=["red"])


def test_plot_crystal_draws_a_molecule():
    """p = 0: one circle per atom and no cell outline."""
    molecule = Crystal(np.zeros((0, 2)), [[0.0, 0.0], [1.2, 0.0], [0.6, 1.0]], [1.0] * 3, ["A"] * 3)
    ax = plot_crystal(molecule, radii=[0.3] * 3, colors=["red"] * 3)
    assert len(ax.patches) == 3


def drawn_atoms(ax):
    return sum(len(c.get_offsets()) for c in ax.collections)


def test_3d_fcc_crystal_atoms_and_cell_edges():
    """fcc with cells = 3 shows the box |x|, |y|, |z| <= 1.5 a: the lattice points (n1, n2, n3) a/2 with
    -3 <= n_i <= 3 and n1 + n2 + n3 even (3 even, 4 odd values: 3^3 + 3 * 4^2 * 3 = 171 atoms); the primitive
    cell has 12 edges."""
    a = 3.6
    fcc = Crystal(primitive_vectors("fcc", a), [[0.0, 0.0, 0.0]], [1.0], ["Cu"])
    ax = plot_crystal(fcc, radii=[0.5], colors=["tab:orange"])

    assert ax.name == "3d"
    assert drawn_atoms(ax) == 171
    assert len(ax.lines) == 12


def test_3d_layer_shows_one_sheet_and_its_cell():
    """A graphene layer (p = 2, d = 3): atoms of both sublattices, the cell outline is a parallelogram (4 edges)."""
    layer = embed(honeycomb(2.46, [12.0] * 2, ["C"] * 2), 3)
    ax = plot_crystal(layer, radii=[0.3] * 2, colors=["red", "blue"])

    assert len(ax.collections) == 2 and drawn_atoms(ax) > 10
    assert len(ax.lines) == 4


def test_3d_molecule_has_atoms_and_no_cell():
    """p = 0 in 3D: the 3 atoms, no cell edges."""
    molecule = Crystal(np.zeros((0, 3)), [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.6, 1.0, 0.5]], [1.0] * 3, ["A"] * 3)
    ax = plot_crystal(molecule, radii=[0.3] * 3, colors=["red"] * 3)
    assert drawn_atoms(ax) == 3 and len(ax.lines) == 0


@pytest.mark.parametrize(
    "radii, colors, cells, ax_2d",
    [([0.3], ["red"] * 2, 3, False), ([0.3] * 2, ["red"], 3, False), ([0.3] * 2, ["red"] * 2, 0, False),
     ([0.3] * 2, ["red"] * 2, 3, True)],
    ids=["too few radii", "too few colors", "no cells", "2D axis for 3D"],
)
def test_plot_crystal_rejects_inconsistent_input(radii, colors, cells, ax_2d):
    """One radius and colour per basis atom, at least one cell, and a 3D axis for a 3D crystal."""
    import matplotlib.pyplot as plt

    layer = embed(honeycomb(2.46, [12.0] * 2, ["C"] * 2), 3)
    with pytest.raises(ValueError):
        plot_crystal(layer, radii=radii, colors=colors, cells=cells, ax=plt.subplots()[1] if ax_2d else None)
