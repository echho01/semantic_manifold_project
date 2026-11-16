"""Fixed-length mixed-radix code for a separated local diagram chart.

The ordered reference features and K are shared. Each endpoint offset is
in [-eps0, eps0]. Their disjoint matching neighborhoods identify features;
the payload encodes all 2K endpoint digits as ONE base-m integer.
"""
from __future__ import annotations

import math
import numpy as np


def bit_length(K: int, m: int) -> int:
    if K < 1 or m < 1:
        raise ValueError("K and m must be positive")
    return (m ** (2 * K) - 1).bit_length()


def encode(reference, diagram, eps0: float, eps1: float):
    reference = np.asarray(reference, float).reshape(-1, 2)
    diagram = np.asarray(diagram, float).reshape(-1, 2)
    if reference.shape != diagram.shape or len(reference) == 0:
        raise ValueError("known feature correspondence and positive K are required")
    if not 0 < eps1 <= eps0:
        raise ValueError("require 0 < eps1 <= eps0")
    offsets = diagram - reference
    if np.max(np.abs(offsets)) > eps0 * (1 + 1e-12):
        raise ValueError("source is outside the local chart")
    m = math.ceil(eps0 / eps1)
    digits = np.clip(np.floor((offsets + eps0) * m / (2 * eps0)), 0, m - 1).astype(int)
    value = 0
    for digit in digits.ravel():
        value = value * m + int(digit)
    length = bit_length(len(reference), m)
    message = format(value, f"0{length}b") if length else ""
    return message, m


def decode(reference, message: str, eps0: float, m: int):
    reference = np.asarray(reference, float).reshape(-1, 2)
    K = len(reference)
    if len(message) != bit_length(K, m) or set(message) - {"0", "1"}:
        raise ValueError("invalid fixed-length payload")
    value = int(message, 2) if message else 0
    if value >= m ** (2 * K):
        raise ValueError("unused binary codeword")
    digits = np.empty(2 * K, int)
    for j in range(2 * K - 1, -1, -1):
        value, digits[j] = divmod(value, m)
    return reference - eps0 + (digits.reshape(K, 2) + 0.5) * (2 * eps0 / m)
