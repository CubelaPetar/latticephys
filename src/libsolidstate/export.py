"""Write results to JSON, e.g. for web pages that plot or animate them."""

import json

import numpy as np


def save_json(path, data):
    """Write a dict to path as JSON; numpy arrays become nested lists, numpy numbers plain numbers."""

    def convert(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"cannot write {type(value).__name__} to JSON")

    with open(path, "w", encoding="utf-8") as file:
        json.dump(data, file, default=convert, ensure_ascii=False)
