"""GPS L1 scintillation monitoring: S4 index and phase scintillation sigma_phi.

Pipeline (50 Hz, INTERVAL 0.02 s):
  1. Epoch-grid validation: times strictly increasing, every spacing an
     integer multiple of 0.02 s within 1 us.
  2. Per-satellite quality control: epochs missing S1C or L1C, or with a
     nonzero L1C LLI, are excluded (kept in the output with their reason)
     and break the arc; time gaps larger than one interval also break arcs.
     Excluded samples are never interpolated back.
  3. Per arc, L1C (continuous carrier phase in cycles, never modded or
     unwrapped) is referenced to the arc's first phase and scaled by 2*pi
     to radians, then detrended by a 6th-order 0.1 Hz Butterworth
     high-pass applied as a causal SOS filter with zero initial state.
     The filter state resets between arcs but runs continuously across
     minute boundaries inside an arc.
  4. Metrics are delivered in 60 s half-open windows aligned to GPS whole
     minutes. A window is valid only with a complete 3000 samples, all
     from a single arc, and a window start at least 60 s after the arc
     start (filter warm-up). Incomplete, tail, cross-arc and warm-up
     windows are returned as failed with their reason.

S4: intensity I = 10^(S1C/10); S4 = std(I)/mean(I) (population std).
sigma_phi: population standard deviation of the detrended phase, radians.
No noise correction is applied to either metric; a high value is not by
itself proof of ionospheric scintillation or receiver fault.
"""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass

import numpy as np
from scipy import signal

from .errors import RejectedContentError
from .geodesy import GPS_EPOCH, epoch_to_seconds
from .rinex import RinexData

FS_HZ = 50.0
INTERVAL_S = 1.0 / FS_HZ          # 0.02 s
TIME_TOL_S = 1e-6                 # 1 us tolerance on the epoch grid
WINDOW_S = 60.0
WINDOW_N = int(WINDOW_S * FS_HZ)  # 3000 samples per window
WARMUP_S = 60.0
MAX_EPOCHS = 18000
MAX_SATELLITES = 4
FILTER_ORDER = 6
CUTOFF_HZ = 0.1
REQUIRED_OBS = ("S1C", "L1C")

_SOS = signal.butter(FILTER_ORDER, CUTOFF_HZ, btype="highpass",
                     fs=FS_HZ, output="sos")


@dataclass
class _Sample:
    time_iso: str
    t_sec: float
    s1c: float | None = None        # dB-Hz
    l1c: float | None = None        # carrier phase, cycles
    problem: str | None = None      # exclusion reason
    arc_no: int = 0                 # 1-based arc number; 0 = not in an arc
    phi_detrended: float | None = None  # high-pass filtered phase, rad


def _validate_epoch_grid(obs: RinexData) -> None:
    prev: float | None = None
    for ep in obs.epochs:
        t = epoch_to_seconds(ep.time)
        if prev is not None:
            step = t - prev
            if step <= 0.0:
                raise RejectedContentError(
                    f"epoch times not strictly increasing at {ep.time.isoformat()}",
                    "rinex", 0)
            k = round(step / INTERVAL_S)
            if abs(step - k * INTERVAL_S) > TIME_TOL_S:
                raise RejectedContentError(
                    f"epoch spacing {step:.6f} s at {ep.time.isoformat()} is not "
                    f"an integer multiple of {INTERVAL_S} s (tolerance 1 us)",
                    "rinex", 0)
        prev = t


def _collect_samples(obs: RinexData, prn: str) -> list[_Sample]:
    samples: list[_Sample] = []
    for ep in obs.epochs:
        t_sec = epoch_to_seconds(ep.time)
        iso = ep.time.isoformat()
        entry = ep.obs.get(prn)
        if entry is None:
            samples.append(_Sample(iso, t_sec,
                                   problem="satellite not observed at epoch"))
            continue
        s1c = entry["S1C"][0] if "S1C" in entry else None
        l1c = entry["L1C"] if "L1C" in entry else None
        missing = [t for t in REQUIRED_OBS if t not in entry]
        if missing:
            samples.append(_Sample(iso, t_sec, s1c=s1c,
                                   l1c=l1c[0] if l1c else None,
                                   problem="missing observation(s): "
                                           + ", ".join(missing)))
            continue
        if l1c[1] != 0:
            samples.append(_Sample(iso, t_sec, s1c=s1c, l1c=l1c[0],
                                   problem=f"L1C loss-of-lock indicator "
                                           f"{l1c[1]} nonzero"))
            continue
        samples.append(_Sample(iso, t_sec, s1c=s1c, l1c=l1c[0]))
    return samples


def _split_arcs(samples: list[_Sample]) -> list[list[_Sample]]:
    """Split ok samples into arcs; breaks on excluded samples and time gaps."""
    arcs: list[list[_Sample]] = []
    cur: list[_Sample] = []
    prev_ok_t: float | None = None
    for s in samples:
        if s.problem is not None:
            if cur:
                arcs.append(cur)
                cur = []
            prev_ok_t = None
            continue
        if prev_ok_t is not None and s.t_sec - prev_ok_t > INTERVAL_S + TIME_TOL_S:
            arcs.append(cur)
            cur = []
        cur.append(s)
        prev_ok_t = s.t_sec
    if cur:
        arcs.append(cur)
    return arcs


def _filter_arc(arc: list[_Sample]) -> None:
    """Detrend one arc: reference phase to arc start, scale to radians,
    causal SOS high-pass with zero initial state (state resets per arc)."""
    phase0 = arc[0].l1c
    phase = np.array([(s.l1c - phase0) * 2.0 * math.pi for s in arc])
    zi = np.zeros((_SOS.shape[0], 2))
    detrended, _ = signal.sosfilt(_SOS, phase, zi=zi)
    for s, value in zip(arc, detrended):
        s.phi_detrended = float(value)


def _minute_floor(t_sec: float) -> float:
    return math.floor(t_sec / WINDOW_S + 1e-9) * WINDOW_S


def _minute_iso(m_sec: float) -> str:
    return (GPS_EPOCH + dt.timedelta(seconds=m_sec)).isoformat()


def _build_windows(ok_samples: list[_Sample],
                   arc_start: dict[int, float]) -> list[dict]:
    """Emit one record per GPS whole-minute 60 s half-open window covering
    the satellite's data span; invalid windows keep their failure reason."""
    windows: list[dict] = []
    if not ok_samples:
        return windows
    by_minute: dict[float, list[_Sample]] = {}
    for s in ok_samples:
        by_minute.setdefault(_minute_floor(s.t_sec), []).append(s)
    m = _minute_floor(ok_samples[0].t_sec)
    last_m = _minute_floor(ok_samples[-1].t_sec)
    while m <= last_m + 1e-9:
        group = by_minute.get(m, [])
        entry: dict = {
            "window_start": _minute_iso(m),
            "status": "failed",
            "reason": None,
            "arc_id": None,
            "n_samples": len(group),
            "s4": None,
            "sigma_phi_rad": None,
        }
        arc_nos = {s.arc_no for s in group}
        if not group:
            entry["reason"] = "no valid samples in window"
        elif len(group) != WINDOW_N:
            entry["reason"] = (f"window has {len(group)} of the required "
                               f"{WINDOW_N} samples (data gap or tail window)")
        elif len(arc_nos) > 1:
            entry["reason"] = "window samples span multiple arcs"
        elif m - arc_start[group[0].arc_no] < WARMUP_S - TIME_TOL_S:
            entry["arc_id"] = group[0].arc_no
            entry["reason"] = (f"filter warm-up: window start less than "
                               f"{WARMUP_S:.0f} s after arc start")
        else:
            s1c = np.array([s.s1c for s in group])
            intensity = np.power(10.0, s1c / 10.0)
            mean_i = float(np.mean(intensity))
            phi = np.array([s.phi_detrended for s in group])
            entry.update(
                status="ok",
                arc_id=group[0].arc_no,
                s4=float(np.std(intensity) / mean_i),
                sigma_phi_rad=float(np.std(phi)),
            )
        windows.append(entry)
        m += WINDOW_S
    return windows


def compute_scintillation(obs: RinexData) -> dict:
    """Full scintillation pipeline; returns the JSON-ready result dict."""
    if obs.interval_s is None:
        raise RejectedContentError(
            "INTERVAL header missing; scintillation monitoring requires "
            "a nominal interval of 0.02 s", "rinex", 0)
    if abs(obs.interval_s - INTERVAL_S) > TIME_TOL_S:
        raise RejectedContentError(
            f"INTERVAL must be {INTERVAL_S} s (50 Hz), got "
            f"{obs.interval_s} s", "rinex", 0)
    _validate_epoch_grid(obs)

    all_prns = sorted({prn for ep in obs.epochs for prn in ep.obs})
    if len(all_prns) > MAX_SATELLITES:
        raise RejectedContentError(
            f"{len(all_prns)} GPS satellites in file, more than the "
            f"{MAX_SATELLITES} supported", "rinex", 0)

    out_sats: dict[str, dict] = {}
    for prn in all_prns:
        samples = _collect_samples(obs, prn)
        arcs = _split_arcs(samples)
        arc_start: dict[int, float] = {}
        arc_infos: list[dict] = []
        for arc_no, arc in enumerate(arcs, start=1):
            arc_id = f"{prn}-A{arc_no}"
            arc_start[arc_no] = arc[0].t_sec
            _filter_arc(arc)
            for s in arc:
                s.arc_no = arc_no
            arc_infos.append({
                "arc_id": arc_id,
                "start_time": arc[0].time_iso,
                "end_time": arc[-1].time_iso,
                "n_samples": len(arc),
            })
        ok_samples = [s for s in samples if s.problem is None]
        windows = _build_windows(ok_samples, arc_start)
        for w in windows:
            if isinstance(w["arc_id"], int):
                w["arc_id"] = f"{prn}-A{w['arc_id']}"
        records = [{
            "time": s.time_iso,
            "status": "ok" if s.problem is None else "excluded",
            "reason": s.problem,
            "s1c_dbhz": s.s1c,
            "l1c_cycles": s.l1c,
            "arc_id": f"{prn}-A{s.arc_no}" if s.arc_no else None,
        } for s in samples]
        out_sats[prn] = {
            "n_arcs": len(arcs),
            "arcs": arc_infos,
            "records": records,
            "windows": windows,
        }

    return {
        "n_epochs": len(obs.epochs),
        "interval_s": INTERVAL_S,
        "sample_rate_hz": FS_HZ,
        "filter": {
            "type": "Butterworth high-pass",
            "order": FILTER_ORDER,
            "cutoff_hz": CUTOFF_HZ,
            "realization": ("causal SOS, zero initial state, state reset "
                            "per arc, continuous across minute boundaries"),
        },
        "window": {
            "length_s": WINDOW_S,
            "alignment": "GPS whole minute, half-open [t, t+60 s)",
            "n_samples_required": WINDOW_N,
            "warmup_s": WARMUP_S,
        },
        "satellites": out_sats,
        "notes": [
            "S4 = std(I)/mean(I) with intensity I = 10^(S1C/10); "
            "population standard deviation.",
            "sigma_phi = population standard deviation of the detrended "
            "carrier phase in radians; L1C cycles referenced to arc start, "
            "scaled by 2*pi, 6th-order 0.1 Hz Butterworth high-pass.",
            "No noise correction applied to S4 or sigma_phi; the metrics "
            "alone do not establish ionospheric scintillation or a "
            "receiver fault.",
        ],
    }

