"""Regenerate the README images: run an example headless and save one of its figures.

Run from the repo root: .venv/bin/python docs/images/make_images.py
"""

import runpy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).parent.parent.parent
IMAGES = {  # image name: (example, number of the figure it creates)
    "nacl_crystal.png": ("crystals_3d/nacl.py", 2),
    "graphene_layer.png": ("layers/graphene_layer.py", 1),
    "graphene_structure_factor.png": ("lattices_2d/structure_factor_2d.py", 1),
    "copper.png": ("crystals_3d/copper.py", 1),
}


def main():
    plt.show = lambda: None
    for name, (example, figure) in IMAGES.items():
        runpy.run_path(str(ROOT / "examples" / example), run_name="__main__")
        plt.figure(figure).savefig(Path(__file__).parent / name, dpi=100)
        plt.close("all")


if __name__ == "__main__":
    main()
