"""Test helpers."""
import dataclasses
import math

import numpy as np

REL = 1e-11
ABS = 1e-11


def assert_same(a, b, path="model", ignore=()):
    """Deep comparison of two resolved models (floats to a relative 1e-11)."""
    if dataclasses.is_dataclass(a) and not isinstance(a, type):
        assert type(a) is type(b), path
        for f in dataclasses.fields(a):
            if f.name in ignore:
                continue
            assert_same(getattr(a, f.name), getattr(b, f.name), f"{path}.{f.name}", ignore)
    elif isinstance(a, dict):
        assert list(a) == list(b), path
        for k in a:
            assert_same(a[k], b[k], f"{path}[{k!r}]", ignore)
    elif isinstance(a, (list, tuple)):
        assert isinstance(b, (list, tuple)) and len(a) == len(b), path
        for k, (x, y) in enumerate(zip(a, b)):
            assert_same(x, y, f"{path}[{k}]", ignore)
    elif isinstance(a, np.ndarray):
        assert a.shape == np.shape(b), path
        for k, (x, y) in enumerate(zip(a.ravel(), np.asarray(b).ravel())):
            assert_same(float(x), float(y), f"{path}[{k}]", ignore)
    elif isinstance(a, float) or isinstance(b, float):
        assert a is not None and b is not None, path
        if math.isinf(a) or math.isinf(b):
            assert a == b, path
        else:
            assert abs(a - b) <= max(ABS, REL * max(abs(a), abs(b))), f"{path}: {a!r} != {b!r}"
    else:
        assert a == b, f"{path}: {a!r} != {b!r}"
