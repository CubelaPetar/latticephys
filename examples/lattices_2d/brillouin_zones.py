"""First Brillouin zones of the five 2D Bravais lattices, with reciprocal lattice points and b1, b2."""

import matplotlib.pyplot as plt

from libsolidstate.lattice.bravais2d import primitive_vectors
from libsolidstate.lattice.plot import plot_brillouin_zone

LATTICES = [
    ("Oblique (a = 3, b = 4 Å, 75°)", primitive_vectors("oblique", 3.0, 4.0, 75.0)),
    ("Rectangular (a = 3, b = 4 Å)", primitive_vectors("rectangular", 3.0, 4.0)),
    ("Centered rectangular (a = 3, b = 4 Å)", primitive_vectors("centered rectangular", 3.0, 4.0)),
    ("Square (a = 3 Å)", primitive_vectors("square", 3.0)),
    ("Hexagonal (a = 3 Å)", primitive_vectors("hexagonal", 3.0)),
]


def main():
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for axis, (title, vectors) in zip(axes.flat, LATTICES):
        plot_brillouin_zone(vectors, ax=axis)
        axis.set_title(title)
    axes.flat[-1].axis("off")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
