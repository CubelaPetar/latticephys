"""Every example runs headless: plots go to the Agg backend, animate() is skipped, exports go to a temp dir."""

import runpy
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pytest

import latticephys.phonons.display as display

matplotlib.use("Agg")
EXAMPLES = sorted((Path(__file__).parent.parent / "examples").glob("*/*.py"))


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda path: f"{path.parent.name}/{path.name}")
def test_example_runs(path, tmp_path, monkeypatch):
    monkeypatch.setattr(display, "animate", lambda *args, **kwargs: None)
    monkeypatch.setattr(plt, "show", lambda: None)
    monkeypatch.setattr(sys, "argv", [str(path), str(tmp_path)])
    runpy.run_path(str(path), run_name="__main__")
    plt.close("all")
