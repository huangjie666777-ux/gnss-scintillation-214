"""Dual-frequency GPS ionosphere TEC monitoring.

Geometry-free combinations per satellite:
  code  GF = (P2 - P1) - c * DCB(P1-P2)   [m]
  phase GF = L1*lam1 - L2*lam2            [m]
Both equal 40.3 * STEC * (1/f2^2 - 1/f1^2) plus (for phase) an arc-constant
ambiguity offset. Each continuous arc is leveled by the arc mean of
(code GF - phase GF); arcs never share offsets. Slant TEC in TECU keeps its
sign (negative values are preserved).

Vertical TEC uses the thin-shell approximation on a sphere of radius
SHELL_RADIUS_M: VTEC = STEC * cos(theta), theta = angle between the
station->satellite ray and the shell normal (geocentric radial) at the
ionosphere pierce point.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .geodesy import C_LIGHT, epoch_to_seconds
from .interp import satellite_position
from .rinex import RinexData
from .sp3 import Sp3Data

F1 = 1575.42e6   # GPS L1 frequency, Hz
F2 = 1227.60e6   # GPS L2 frequency, Hz
LAM1 = C_LIGHT / F1
LAM2 = C_LIGHT / F2
# meters of geometry-free combination per (electron/m^2) of slant TEC
K_TEC_M = 40.3 * (1.0 / F2**2 - 1.0 / F1**2)
TECU = 1e16

SHELL_RADIUS_M = 6821e3
MIN_ARC_SAMPLES = 3
MAX_GAP_S = 120.0
MAX_PHASE_JUMP_M = 1.0
MIN_ELEVATION_DEG = 10.0
REQUIRED_OBS = ("C1C", "C2W", "L1C", "L2W")


@dataclass
class _Sample:
    time_iso: str
    t_sec: float
    gf_code: float | None = None
    gf_phase: float | None = None
    problem: str | None = None  # data-level break reason


def _collect_samples(obs: RinexData, prn: str, bias_ns: float | None) -> list[_Sample]:
    samples: list[_Sample] = []
    for ep in obs.epochs:
        t_sec = epoch_to_seconds(ep.time)
        iso = ep.time.isoformat()
        entry = ep.obs.get(prn)
        if entry is None:
            samples.append(_Sample(iso, t_sec, problem="satellite not observed at epoch"))
            continue
        missing = [o for o in REQUIRED_OBS if o not in entry]
        if missing:
            samples.append(_Sample(
                iso, t_sec, problem=f"missing observation(s): {', '.join(missing)}"))
            continue
        lli_bad = [o for o in ("L1C", "L2W") if entry[o][1] != 0]
        if lli_bad:
            samples.append(_Sample(
                iso, t_sec, problem=f"loss-of-lock indicator set on {', '.join(lli_bad)}"))
            continue
        if bias_ns is None:
            samples.append(_Sample(
                iso, t_sec, problem="no merged code bias (DCB) provided for satellite"))
            continue
        p1 = entry["C1C"][0]
        p2 = entry["C2W"][0]
        l1_m = entry["L1C"][0] * LAM1
        l2_m = entry["L2W"][0] * LAM2
        gf_code = (p2 - p1) - C_LIGHT * bias_ns * 1e-9
        gf_phase = l1_m - l2_m
        samples.append(_Sample(iso, t_sec, gf_code=gf_code, gf_phase=gf_phase))
    return samples


def _split_arcs(samples: list[_Sample]) -> list[list[_Sample]]:
    """Split into arcs; samples with data problems are returned separately."""
    arcs: list[list[_Sample]] = []
    cur: list[_Sample] = []
    prev_ok: _Sample | None = None
    for s in samples:
        if s.problem is not None:
            if cur:
                arcs.append(cur)
                cur = []
            prev_ok = None
            continue
        if prev_ok is not None:
            gap = s.t_sec - prev_ok.t_sec
            jump = abs(s.gf_phase - prev_ok.gf_phase)
            if gap > MAX_GAP_S or jump > MAX_PHASE_JUMP_M:
                arcs.append(cur)
                cur = []
        cur.append(s)
        prev_ok = s
    if cur:
        arcs.append(cur)
    return arcs


def _geometry(station_ecef: np.ndarray, sat_ecef: np.ndarray):
    """Elevation (deg), pierce point geocentric lat/lon (deg), mapping cos."""
    up = station_ecef / np.linalg.norm(station_ecef)  # station radial = zenith
    ray = sat_ecef - station_ecef
    dist = float(np.linalg.norm(ray))
    d = ray / dist
    sin_el = float(np.dot(d, up))
    elev = math.degrees(math.asin(max(-1.0, min(1.0, sin_el))))
    if elev < MIN_ELEVATION_DEG:
        return elev, None, None
    # intersect ray station + t*d with sphere of radius SHELL_RADIUS_M
    b = float(np.dot(station_ecef, d))
    c = float(np.dot(station_ecef, station_ecef) - SHELL_RADIUS_M**2)
    disc = b * b - c
    if disc <= 0.0:
        return elev, None, None
    t = -b + math.sqrt(disc)
    if t <= 0.0:
        return elev, None, None
    pierce = station_ecef + t * d
    normal = pierce / np.linalg.norm(pierce)
    cos_theta = float(np.dot(d, normal))
    lat = math.degrees(math.asin(max(-1.0, min(1.0, pierce[2] / SHELL_RADIUS_M))))
    lon = math.degrees(math.atan2(float(pierce[1]), float(pierce[0])))
    return elev, {"lat_deg": lat, "lon_deg": lon}, cos_theta


def compute_tec(obs: RinexData, sp3: Sp3Data, station_ecef: np.ndarray,
                biases_ns: dict[str, float]) -> dict:
    """Per-satellite, per-epoch slant/vertical TEC with arc and pierce info."""
    all_prns = sorted({prn for ep in obs.epochs for prn in ep.obs})
    out_sats: dict[str, dict] = {}
    for prn in all_prns:
        bias = biases_ns.get(prn)
        samples = _collect_samples(obs, prn, bias)
        arcs = _split_arcs(samples)
        arc_of: dict[int, tuple[int, float]] = {}  # id(sample) -> (arc_no, level offset)
        short_arc_ids: set[int] = set()
        for arc_no, arc in enumerate(arcs, start=1):
            if len(arc) < MIN_ARC_SAMPLES:
                for s in arc:
                    short_arc_ids.add(id(s))
                continue
            offset = float(np.mean([s.gf_code - s.gf_phase for s in arc]))
            for s in arc:
                arc_of[id(s)] = (arc_no, offset)

        rec = sp3.sats.get(prn)
        records: list[dict] = []
        for s in samples:
            record = {
                "time": s.time_iso,
                "arc_id": None,
                "status": "ok",
                "reason": None,
                "stec_tecu": None,
                "vtec_tecu": None,
                "elevation_deg": None,
                "pierce_point": None,
            }
            if s.problem is not None:
                record["status"] = "invalid"
                record["reason"] = s.problem
                records.append(record)
                continue
            if id(s) in short_arc_ids:
                record["status"] = "invalid"
                record["reason"] = (f"arc has fewer than {MIN_ARC_SAMPLES} "
                                    "complete samples; no level offset")
                records.append(record)
                continue
            arc_no, offset = arc_of[id(s)]
            arc_id = f"{prn}-A{arc_no}"
            record["arc_id"] = arc_id
            leveled_phase = s.gf_phase + offset
            stec_tecu = leveled_phase / (K_TEC_M * TECU)
            record["stec_tecu"] = stec_tecu

            if rec is None:
                record["status"] = "excluded"
                record["reason"] = "satellite not present in SP3 file"
                records.append(record)
                continue
            pos, reason = satellite_position(rec, s.t_sec)
            if pos is None:
                record["status"] = "excluded"
                record["reason"] = reason
                records.append(record)
                continue
            elev, pierce, cos_theta = _geometry(station_ecef, pos)
            record["elevation_deg"] = elev
            if pierce is None:
                record["status"] = "excluded"
                record["reason"] = (f"elevation {elev:.2f} deg below "
                                    f"{MIN_ELEVATION_DEG:.0f} deg cutoff")
                records.append(record)
                continue
            record["pierce_point"] = pierce
            record["vtec_tecu"] = stec_tecu * cos_theta
            records.append(record)

        out_sats[prn] = {
            "bias_ns": bias,
            "n_arcs": len(arcs),
            "records": records,
        }
    return {
        "shell_radius_m": SHELL_RADIUS_M,
        "min_elevation_deg": MIN_ELEVATION_DEG,
        "satellites": out_sats,
    }
