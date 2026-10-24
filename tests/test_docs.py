"""The Python code blocks in the README and the guide run (plots to the Agg backend)."""

import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pytest

matplotlib.use("Agg")
ROOT = Path(__file__).parent.parent
DOCS = [ROOT / "README.md", *sorted((ROOT / "docs" / "guide").glob("*.md"))]


@pytest.mark.parametrize("path", DOCS, ids=lambda path: path.name)
def test_code_blocks_run(path, monkeypatch):
    monkeypatch.setattr(plt, "show", lambda: None)
    blocks = re.findall(r"```python\n(.*?)```", path.read_text(), re.S)
    exec("\n".join(blocks), {})  # blocks of one page build on each other
    plt.close("all")
