import json

import numpy as np
import pytest

from latticephys.export import save_json


def test_round_trip_with_numpy(tmp_path):
    data = {"labels": ["Γ", "K"], "frequencies": np.array([[0.0, 1.5]]), "omega": np.float64(2.5), "cells": np.int64(6)}
    save_json(tmp_path / "out.json", data)
    loaded = json.loads((tmp_path / "out.json").read_text(encoding="utf-8"))
    assert loaded == {"labels": ["Γ", "K"], "frequencies": [[0.0, 1.5]], "omega": 2.5, "cells": 6}


def test_unsupported_type_raises(tmp_path):
    with pytest.raises(TypeError):
        save_json(tmp_path / "out.json", {"value": object()})
