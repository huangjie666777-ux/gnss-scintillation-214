"""Satellite position (Lagrange 8-node) and clock (linear) interpolation."""

from __future__ import annotations

import bisect

import numpy as np

from .sp3 import SatRecord

N_LAGRANGE = 8


def _window(times: list[float], t: float, n: int) -> list[int] | None:
    """Indices of n consecutive nodes bracketing t, or None if out of range."""
    idx = bisect.bisect_right(times, t)
    lo = idx - n // 2
    lo = max(0, min(lo, len(times) - n))
    if lo < 0 or lo + n > len(times):
        return None
    if t < times[lo] or t > times[lo + n - 1]:
        return None
    return list(range(lo, lo + n))


def satellite_position(rec: SatRecord, t: float) -> tuple[np.ndarray | None, str | None]:
    """8-node Lagrange interpolation of position at time t (s since GPS epoch).

    Returns (position_m, None) or (None, reason) without extrapolating or
    crossing nodes with invalid (zero) coordinates.
    """
    win = _window(rec.times, t, N_LAGRANGE)
    if win is None:
        return None, "outside SP3 coverage (needs 8-node window, no extrapolation)"
    if not all(rec.pos_valid[k] for k in win):
        return None, "zero/invalid coordinates inside interpolation window"
    # never interpolate across a gap in the ephemeris time nodes
    nominal = min(b - a for a, b in zip(rec.times, rec.times[1:]))
    for a, b in zip(win, win[1:]):
        if rec.times[b] - rec.times[a] > 1.5 * nominal:
            return None, "gap in SP3 time nodes inside interpolation window"
    ts = [rec.times[k] for k in win]
    pos = np.zeros(3)
    for j, k in enumerate(win):
        w = 1.0
        for m, tm in enumerate(ts):
            if m != j:
                w *= (t - tm) / (ts[j] - tm)
        pos += w * np.asarray(rec.pos[k])
    return pos, None


def satellite_clock(rec: SatRecord, t: float) -> tuple[float | None, str | None]:
    """Linear interpolation of satellite clock (seconds) between adjacent nodes.

    Does not cross missing clock values or extrapolate.
    """
    times = rec.times
    idx = bisect.bisect_right(times, t)
    if idx == 0 or idx >= len(times) + 1:
        return None, "outside SP3 clock coverage (no extrapolation)"
    if idx == len(times):
        return None, "outside SP3 clock coverage (no extrapolation)"
    a, b = idx - 1, idx
    if not (rec.clk_valid[a] and rec.clk_valid[b]):
        return None, "missing satellite clock at adjacent node"
    t0, t1 = times[a], times[b]
    if t1 == t0:
        return None, "degenerate SP3 time nodes"
    f = (t - t0) / (t1 - t0)
    return rec.clk[a] * (1.0 - f) + rec.clk[b] * f, None
