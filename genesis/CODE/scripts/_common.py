"""Shared paths for the reproduction scripts."""
import os
import sys
import json

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
ROOT = os.path.dirname(CODE)
FIG = os.path.join(ROOT, "FIGURES")
DATA = os.path.join(ROOT, "DATA")
sys.path.insert(0, CODE)
os.makedirs(FIG, exist_ok=True)
os.makedirs(DATA, exist_ok=True)


def _default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def dump(name, obj):
    path = os.path.join(DATA, name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_default)
    print("wrote", os.path.relpath(path, ROOT))
    return path
