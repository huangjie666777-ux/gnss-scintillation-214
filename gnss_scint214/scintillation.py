"""GPS L1 scintillation monitoring from S1C and continuous L1C phase."""

from __future__ import annotations

import datetime as dt
import math
from dataclasses import dataclass

import numpy as np
from scipy.signal import butter, sosfilt

from .geodesy import epoch_to_seconds
from .rinex import RinexData

SAMPLE_RATE_HZ = 50.0
WINDOW_S = 60.0
WINDOW_SAMPLES = 3000
WARMUP_S = 60.0
S1C_OBS = "S1C"
L1C_OBS = "L1C"
REQUIRED_OBS = (S1C_OBS, L1C_OBS)


@dataclass
class _Sample:
    epoch_index: int
    time: dt.datetime
    t_sec: float
    s1c_dbhz: float
    phase_cycles: float
    filtered_phase_rad: float | None = None


def _minute_floor(t: dt.datetime) -> dt.datetime:
    return t.replace(second=0, microsecond=0)


def _build_arcs(obs: RinexData, prn: str) -> tuple[list[list[_Sample]], list[dict]]:
    quality: list[dict] = []
    raw: list[_Sample | None] = []
    breaks_before: set[int] = set()
    prev_good_time: dt.datetime | None = None

    for idx, epoch in enumerate(obs.epochs):
        entry = epoch.obs.get(prn)
        reason = None
        sample: _Sample | None = None
        observed_s1c = entry[S1C_OBS][0] if entry and S1C_OBS in entry else None
        observed_l1c = entry[L1C_OBS][0] if entry and L1C_OBS in entry else None
        if entry is None:
            reason = "satellite record missing at epoch"
        elif S1C_OBS not in entry:
            reason = f"missing {S1C_OBS}"
        elif L1C_OBS not in entry:
            reason = f"missing {L1C_OBS}"
        elif entry[L1C_OBS][1] != 0:
            reason = f"{L1C_OBS} LLI={entry[L1C_OBS][1]} (loss-of-lock indicator nonzero)"
        else:
            s1c, _ = entry[S1C_OBS]
            phase, _ = entry[L1C_OBS]
            sample = _Sample(idx, epoch.time, epoch_to_seconds(epoch.time), s1c, phase)
            if prev_good_time is not None:
                gap_s = (epoch.time - prev_good_time).total_seconds()
                if abs(gap_s - 1.0 / SAMPLE_RATE_HZ) > 1e-6:
                    breaks_before.add(idx)
            prev_good_time = epoch.time

        quality.append({
            "time": epoch.time.isoformat(),
            "arc_id": None,
            "status": "ok" if reason is None else "excluded",
            "reason": reason,
            "s1c_dbhz": observed_s1c,
            "l1c_cycles": observed_l1c,
        })
        raw.append(sample)

    arcs: list[list[_Sample]] = []
    current: list[_Sample] = []
    for sample in raw:
        if sample is None:
            if current:
                arcs.append(current)
                current = []
            continue
        if current:
            missing_epochs = sample.epoch_index != current[-1].epoch_index + 1
            if missing_epochs or sample.epoch_index in breaks_before:
                arcs.append(current)
                current = []
        current.append(sample)
    if current:
        arcs.append(current)

    sos = butter(6, 0.1, btype="highpass", fs=SAMPLE_RATE_HZ, output="sos")
    for arc_no, arc in enumerate(arcs, start=1):
        phase = np.array([s.phase_cycles for s in arc], dtype=float)
        phase_rad = (phase - phase[0]) * 2.0 * math.pi
        filtered = sosfilt(sos, phase_rad, axis=0, zi=None)
        for sample, value in zip(arc, filtered):
            sample.filtered_phase_rad = float(value)
        for sample in arc:
            quality[sample.epoch_index]["arc_id"] = f"{prn}-A{arc_no}"

    return arcs, quality


def _window_starts(obs: RinexData) -> list[dt.datetime]:
    first = _minute_floor(obs.epochs[0].time)
    last_minute = _minute_floor(obs.epochs[-1].time)
    n = int((last_minute - first).total_seconds() // 60) + 1
    return [first + dt.timedelta(seconds=60 * i) for i in range(n)]


def _evaluate_window(prn: str, arc_no: int, arc: list[_Sample],
                     start: dt.datetime, end: dt.datetime,
                     invalid_in_window: list[str],
                     file_end: dt.datetime) -> dict:
    arc_id = f"{prn}-A{arc_no}"
    points = [s for s in arc if start <= s.time < end]
    n_samples = len(points)
    record = {
        "satellite": prn,
        "time": start.isoformat(),
        "arc_id": arc_id if n_samples else None,
        "status": "failed",
        "reason": None,
        "n_samples": n_samples,
        "s4": None,
        "sigma_phi_rad": None,
    }
    reasons = list(dict.fromkeys(invalid_in_window))
    if n_samples == WINDOW_SAMPLES:
        if reasons:
            record["reason"] = "; ".join(reasons)
            return record
        warmup_s = (start - arc[0].time).total_seconds()
        if warmup_s + 1e-6 < WARMUP_S:
            reasons.append(f"insufficient filter warm-up: window start is "
                           f"{warmup_s:.3f} s after arc start, requires >=60 s")
        else:
            intensity = np.power(10.0, np.array([s.s1c_dbhz for s in points]) / 10.0)
            mean_i = float(np.mean(intensity))
            phase = np.array([s.filtered_phase_rad for s in points], dtype=float)
            if mean_i <= 0.0:
                reasons.append("non-positive mean S1C-derived intensity")
            elif not np.all(np.isfinite(phase)):
                reasons.append("non-finite filtered phase")
            else:
                record.update({
                    "status": "ok",
                    "reason": None,
                    "s4": float(np.std(intensity, ddof=0) / mean_i),
                    "sigma_phi_rad": float(np.std(phase, ddof=0)),
                })
                return record
    else:
        reasons.append(f"incomplete window: {n_samples}/{WINDOW_SAMPLES} samples")
        if end > file_end + dt.timedelta(microseconds=20000):
            reasons.append("trailing window extends beyond available observations")
    record["reason"] = "; ".join(reasons)
    return record


def compute_scintillation(obs: RinexData) -> dict:
    all_prns = sorted({prn for epoch in obs.epochs for prn in epoch.obs})
    starts = _window_starts(obs)
    file_end = obs.epochs[-1].time + dt.timedelta(milliseconds=20)
    satellites: dict[str, dict] = {}

    for prn in all_prns:
        arcs, quality = _build_arcs(obs, prn)
        windows: list[dict] = []
        for start in starts:
            end = start + dt.timedelta(seconds=WINDOW_S)
            invalid_counts: dict[str, int] = {}
            for q in quality:
                qtime = dt.datetime.fromisoformat(q["time"])
                if start <= qtime < end and q["status"] != "ok":
                    invalid_counts[q["reason"]] = invalid_counts.get(q["reason"], 0) + 1
            invalid = [f"{count} excluded sample(s): {reason}"
                       for reason, count in invalid_counts.items()]
            covering = [(no, arc) for no, arc in enumerate(arcs, start=1)
                        if any(start <= s.time < end for s in arc)]
            if len(covering) == 1:
                arc_no, arc = covering[0]
            elif len(covering) > 1:
                arc_no, arc = covering[0]
                invalid.append("window crosses an arc boundary")
            else:
                arc_no, arc = 0, []
                if not invalid:
                    invalid.append(f"no continuous same-arc samples in window")
            windows.append(_evaluate_window(prn, arc_no, arc, start, end,
                                            invalid, file_end))
        satellites[prn] = {
            "n_arcs": len(arcs),
            "quality": quality,
            "windows": windows,
        }

    return {
        "sample_rate_hz": SAMPLE_RATE_HZ,
        "interval_s": 1.0 / SAMPLE_RATE_HZ,
        "window_seconds": WINDOW_S,
        "window_samples": WINDOW_SAMPLES,
        "warmup_seconds": WARMUP_S,
        "highpass": {"filter": "Butterworth", "order": 6,
                     "cutoff_hz": 0.1, "form": "SOS", "causal": True,
                     "initial_state": "zero", "reset": "per arc"},
        "satellites": satellites,
    }
