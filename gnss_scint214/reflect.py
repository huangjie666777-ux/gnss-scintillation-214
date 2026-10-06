"""Single-frequency GNSS-IR reflectometry from L1 C/A SNR (S1C).

Model: a planar, static reflector (water surface) at vertical distance H
below the antenna phase center produces an interference oscillation in the
SNR series with phase 4*pi*H*sin(elevation)/lambda1. The reflector height H
is found by a least-squares grid search over H (unevenly sampled sin(elev)
data cannot be FFT'd directly); water level = Z - H where Z is the antenna
phase-center height above the gauge zero.

Assumptions and limitations: static planar reflector, no tropospheric
refraction correction, specular point assumed directly below the antenna.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geodesy import C_LIGHT, elevation_deg, epoch_to_seconds
from .interp import satellite_position
from .rinex import RinexData
from .sp3 import Sp3Data

F1 = 1575.42e6          # GPS L1 frequency, Hz
LAM1 = C_LIGHT / F1     # L1 wavelength, m

MIN_ELEVATION_DEG = 5.0
MAX_ELEVATION_DEG = 25.0
MAX_GAP_S = 120.0
MIN_ARC_POINTS = 12
MIN_ARC_SPAN_DEG = 5.0
MAX_GRID_POINTS = 5001
SNR_OBS = "S1C"


@dataclass
class _Sample:
    time_iso: str
    t_sec: float
    s1c: float | None = None          # dB-Hz
    elevation: float | None = None    # deg
    problem: str | None = None        # break/filter reason


def _collect_samples(obs: RinexData, sp3: Sp3Data, prn: str,
                     station_ecef: np.ndarray,
                     lat_deg: float, lon_deg: float) -> list[_Sample]:
    samples: list[_Sample] = []
    rec = sp3.sats.get(prn)
    for ep in obs.epochs:
        t_sec = epoch_to_seconds(ep.time)
        iso = ep.time.isoformat()
        entry = ep.obs.get(prn)
        if entry is None:
            samples.append(_Sample(iso, t_sec,
                                   problem="satellite not observed at epoch"))
            continue
        if SNR_OBS not in entry:
            samples.append(_Sample(iso, t_sec,
                                   problem=f"missing observation(s): {SNR_OBS}"))
            continue
        s1c = entry[SNR_OBS][0]
        if rec is None:
            samples.append(_Sample(iso, t_sec, s1c=s1c,
                                   problem="satellite not present in SP3 file"))
            continue
        pos, reason = satellite_position(rec, t_sec)
        if pos is None:
            samples.append(_Sample(iso, t_sec, s1c=s1c, problem=reason))
            continue
        elev = elevation_deg(station_ecef, lat_deg, lon_deg, pos)
        if not (MIN_ELEVATION_DEG <= elev <= MAX_ELEVATION_DEG):
            samples.append(_Sample(
                iso, t_sec, s1c=s1c, elevation=elev,
                problem=(f"elevation {elev:.2f} deg outside "
                         f"[{MIN_ELEVATION_DEG:.0f}, {MAX_ELEVATION_DEG:.0f}] "
                         "deg window")))
            continue
        samples.append(_Sample(iso, t_sec, s1c=s1c, elevation=elev))
    return samples


def _split_arcs(samples: list[_Sample]) -> list[list[_Sample]]:
    """Split ok samples into elevation-monotonic arcs.

    Breaks on filtered/missing samples (never crosses a filter gap), on time
    gaps > MAX_GAP_S, and on rising/setting reversals.
    """
    arcs: list[list[_Sample]] = []
    cur: list[_Sample] = []
    prev: _Sample | None = None
    direction = 0  # +1 rising, -1 setting, 0 not yet established
    for s in samples:
        if s.problem is not None:
            if cur:
                arcs.append(cur)
                cur = []
            prev = None
            direction = 0
            continue
        if prev is not None:
            if s.t_sec - prev.t_sec > MAX_GAP_S:
                arcs.append(cur)
                cur = []
                direction = 0
            else:
                d = (s.elevation > prev.elevation) - (s.elevation < prev.elevation)
                if d != 0:
                    if direction == 0:
                        direction = d
                    elif d != direction:
                        arcs.append(cur)
                        cur = []
                        direction = 0
        cur.append(s)
        prev = s
    if cur:
        arcs.append(cur)
    return arcs


def _invert_arc(arc: list[_Sample], h_grid: np.ndarray) -> dict:
    """Grid-search reflector height for one arc; returns arc result dict."""
    result: dict = {
        "status": "failed",
        "reason": None,
        "warning": None,
        "reflector_height_m": None,
        "amplitude": None,
        "residual_rms": None,
        "height_curve": None,
    }
    x = np.array([math.sin(math.radians(s.elevation)) for s in arc])
    snr = np.array([s.s1c for s in arc])
    amp = np.power(10.0, snr / 20.0)  # dB-Hz -> linear amplitude

    # remove a quadratic trend in sin(elevation)
    design_trend = np.column_stack([x * x, x, np.ones_like(x)])
    coeffs, _, rank, _ = np.linalg.lstsq(design_trend, amp, rcond=None)
    if rank < 3:
        result["reason"] = "detrend design matrix rank deficient"
        return result
    y = amp - design_trend @ coeffs
    if float(np.std(y)) <= 1e-9:
        result["reason"] = "no residual oscillation after quadratic detrend"
        return result

    best_idx = -1
    best_rss = math.inf
    best_coef: np.ndarray | None = None
    curve: list[dict] = []
    n_rank_ok = 0
    for i, h in enumerate(h_grid):
        phi = 4.0 * math.pi * float(h) * x / LAM1
        design = np.column_stack([np.cos(phi), np.sin(phi), np.ones_like(x)])
        coef, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
        if rank < 3:
            curve.append({"height_m": float(h), "rss": None})
            continue
        n_rank_ok += 1
        rss = float(np.sum((y - design @ coef) ** 2))
        curve.append({"height_m": float(h), "rss": rss})
        if rss < best_rss:  # strict: ties keep the smaller H
            best_rss = rss
            best_idx = i
            best_coef = coef
    result["height_curve"] = curve
    if n_rank_ok == 0 or best_coef is None:
        result["reason"] = ("rank-deficient fit at every grid height "
                            "(insufficient sin(elevation) spread)")
        return result

    result["status"] = "ok"
    result["reflector_height_m"] = float(h_grid[best_idx])
    result["amplitude"] = float(math.hypot(best_coef[0], best_coef[1]))
    result["residual_rms"] = math.sqrt(best_rss / len(arc))
    if best_idx == 0 or best_idx == len(h_grid) - 1:
        result["warning"] = ("optimum at search-grid boundary; "
                             "true reflector height may lie outside "
                             f"[{h_grid[0]:.3f}, {h_grid[-1]:.3f}] m")
    return result


def build_grid(h_min: float, h_max: float, h_step: float) -> np.ndarray:
    # tolerate floating-point rounding so a legal upper bound that is an
    # exact multiple of the step (e.g. (2.5-0.5)/0.1 = 19.999...999) is
    # not silently dropped from the search grid
    quotient = (h_max - h_min) / h_step
    n = int(math.floor(quotient + 1e-9)) + 1
    grid = h_min + h_step * np.arange(n)
    if abs(grid[-1] - h_max) <= 1e-9 * max(1.0, abs(h_max)):
        grid[-1] = h_max
    return grid


def compute_reflectometry(obs: RinexData, sp3: Sp3Data,
                          station_ecef: np.ndarray, lat_deg: float,
                          lon_deg: float, antenna_height_m: float,
                          h_min: float, h_max: float, h_step: float) -> dict:
    """Per-satellite records plus per-arc reflector-height inversions."""
    h_grid = build_grid(h_min, h_max, h_step)
    all_prns = sorted({prn for ep in obs.epochs for prn in ep.obs})
    out_sats: dict[str, dict] = {}
    for prn in all_prns:
        samples = _collect_samples(obs, sp3, prn, station_ecef, lat_deg, lon_deg)
        arcs = _split_arcs(samples)
        arc_of: dict[int, int] = {}  # id(sample) -> arc number
        arc_results: list[dict] = []
        for arc_no, arc in enumerate(arcs, start=1):
            for s in arc:
                arc_of[id(s)] = arc_no
            direction = ("rising" if arc[-1].elevation > arc[0].elevation
                         else "setting")
            entry: dict = {
                "arc_id": f"{prn}-A{arc_no}",
                "start_time": arc[0].time_iso,
                "end_time": arc[-1].time_iso,
                "direction": direction,
                "n_samples": len(arc),
                "elevation_span_deg": abs(arc[-1].elevation - arc[0].elevation),
                "sample_source": SNR_OBS,
                "samples": [
                    {"time": s.time_iso,
                     "elevation_deg": s.elevation,
                     "s1c_dbhz": s.s1c}
                    for s in arc
                ],
                "water_level_m": None,
            }
            if len(arc) < MIN_ARC_POINTS:
                entry.update(status="failed",
                             reason=(f"arc has {len(arc)} points, fewer than "
                                     f"the required {MIN_ARC_POINTS}"),
                             warning=None, reflector_height_m=None,
                             amplitude=None, residual_rms=None,
                             height_curve=None)
            elif abs(arc[-1].elevation - arc[0].elevation) < MIN_ARC_SPAN_DEG:
                entry.update(status="failed",
                             reason=(f"elevation span "
                                     f"{abs(arc[-1].elevation - arc[0].elevation):.2f} "
                                     f"deg below the required "
                                     f"{MIN_ARC_SPAN_DEG:.0f} deg"),
                             warning=None, reflector_height_m=None,
                             amplitude=None, residual_rms=None,
                             height_curve=None)
            else:
                entry.update(_invert_arc(arc, h_grid))
                if entry["status"] == "ok":
                    entry["water_level_m"] = (antenna_height_m
                                              - entry["reflector_height_m"])
            arc_results.append(entry)

        records: list[dict] = []
        for s in samples:
            record = {
                "time": s.time_iso,
                "status": "ok" if s.problem is None else "excluded",
                "reason": s.problem,
                "elevation_deg": s.elevation,
                "s1c_dbhz": s.s1c,
                "arc_id": (f"{prn}-A{arc_of[id(s)]}"
                           if s.problem is None and id(s) in arc_of else None),
            }
            records.append(record)

        out_sats[prn] = {
            "n_arcs": len(arcs),
            "records": records,
            "arcs": arc_results,
        }
    return {
        "antenna_height_m": antenna_height_m,
        "grid": {
            "h_min_m": float(h_grid[0]),
            "h_max_m": float(h_grid[-1]),
            "h_step_m": h_step,
            "n_points": len(h_grid),
        },
        "elevation_window_deg": [MIN_ELEVATION_DEG, MAX_ELEVATION_DEG],
        "min_arc_points": MIN_ARC_POINTS,
        "min_arc_span_deg": MIN_ARC_SPAN_DEG,
        "lambda1_m": LAM1,
        "satellites": out_sats,
    }
