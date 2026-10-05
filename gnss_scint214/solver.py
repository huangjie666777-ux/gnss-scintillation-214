"""Per-epoch iterative least-squares receiver position + clock solution."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .geodesy import C_LIGHT, OMEGA_E, ecef_to_geodetic, epoch_to_seconds
from .interp import satellite_clock, satellite_position
from .rinex import EpochObs
from .sp3 import Sp3Data

MAX_ITER = 20
CONV_TOL = 1e-4  # meters


@dataclass
class SatSolution:
    prn: str
    residual_m: float


@dataclass
class EpochResult:
    time: str
    status: str  # "ok" | "failed"
    reason: str | None = None
    ecef_m: list[float] | None = None
    geodetic: dict | None = None
    receiver_clock_s: float | None = None
    used_satellites: list[str] = field(default_factory=list)
    excluded_satellites: list[dict] = field(default_factory=list)
    residuals_m: dict[str, float] = field(default_factory=dict)
    rms_m: float | None = None


def solve_epoch(obs: EpochObs, sp3: Sp3Data,
                initial_ecef: np.ndarray) -> EpochResult:
    t_rx = epoch_to_seconds(obs.time)
    result = EpochResult(time=obs.time.isoformat(), status="failed")

    # --- gather usable satellites: transmit-time position + clock ---
    sats: list[tuple[str, float, np.ndarray, float]] = []  # prn, P, pos, clk
    for prn, pr_obs in sorted(obs.pseudoranges.items()):
        rec = sp3.sats.get(prn)
        if rec is None:
            result.excluded_satellites.append(
                {"prn": prn, "reason": "satellite not present in SP3 file"})
            continue
        # transmit time: t_rx - (P/c) - satellite clock; iterate twice
        t_tx = t_rx - pr_obs / C_LIGHT
        clk, reason = satellite_clock(rec, t_tx)
        if clk is None:
            result.excluded_satellites.append({"prn": prn, "reason": reason})
            continue
        t_tx = t_rx - (pr_obs / C_LIGHT + clk)
        clk, reason = satellite_clock(rec, t_tx)
        if clk is None:
            result.excluded_satellites.append({"prn": prn, "reason": reason})
            continue
        pos, reason = satellite_position(rec, t_tx)
        if pos is None:
            result.excluded_satellites.append({"prn": prn, "reason": reason})
            continue
        sats.append((prn, pr_obs, pos, clk))

    if len(sats) < 4:
        result.reason = (f"only {len(sats)} usable satellites (need >= 4); "
                         f"excluded: {[e['prn'] for e in result.excluded_satellites]}")
        return result

    # --- iterative least squares ---
    x = np.array([*initial_ecef, 0.0], dtype=float)  # x,y,z, c*dt (m)
    converged = False
    for _ in range(MAX_ITER):
        H = np.zeros((len(sats), 4))
        v = np.zeros(len(sats))
        for i, (prn, pr_obs, pos_tx, clk) in enumerate(sats):
            # Earth rotation during signal flight
            tau = np.linalg.norm(pos_tx - x[:3]) / C_LIGHT
            theta = OMEGA_E * tau
            rot = np.array([
                [np.cos(theta), np.sin(theta), 0.0],
                [-np.sin(theta), np.cos(theta), 0.0],
                [0.0, 0.0, 1.0],
            ])
            s = rot @ pos_tx
            r_vec = s - x[:3]
            rho = np.linalg.norm(r_vec)
            predicted = rho + x[3] - C_LIGHT * clk
            v[i] = pr_obs - predicted
            H[i, :3] = -r_vec / rho
            H[i, 3] = 1.0
        # rank check
        if np.linalg.matrix_rank(H) < 4:
            result.reason = "geometry matrix rank deficient"
            return result
        try:
            dx, *_ = np.linalg.lstsq(H, v, rcond=None)
        except np.linalg.LinAlgError:
            result.reason = "least-squares solve failed (singular geometry)"
            return result
        x += dx
        if np.linalg.norm(dx[:3]) < CONV_TOL:
            converged = True
            break
    if not converged:
        result.reason = f"not converged after {MAX_ITER} iterations"
        return result

    # --- final residuals ---
    residuals: dict[str, float] = {}
    ss = 0.0
    for i, (prn, pr_obs, pos_tx, clk) in enumerate(sats):
        tau = np.linalg.norm(pos_tx - x[:3]) / C_LIGHT
        theta = OMEGA_E * tau
        rot = np.array([
            [np.cos(theta), np.sin(theta), 0.0],
            [-np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ])
        s = rot @ pos_tx
        rho = np.linalg.norm(s - x[:3])
        predicted = rho + x[3] - C_LIGHT * clk
        res = pr_obs - predicted
        residuals[prn] = res
        ss += res * res
    rms = float(np.sqrt(ss / len(sats)))

    lat, lon, h = ecef_to_geodetic(x[:3])
    result.status = "ok"
    result.ecef_m = [float(v) for v in x[:3]]
    result.geodetic = {"lat_deg": lat, "lon_deg": lon, "height_m": h}
    result.receiver_clock_s = float(x[3] / C_LIGHT)
    result.used_satellites = [s[0] for s in sats]
    result.residuals_m = residuals
    result.rms_m = rms
    return result


def solve_all(epochs: list[EpochObs], sp3: Sp3Data,
              approx_position: tuple[float, float, float] | None) -> list[EpochResult]:
    """Solve each epoch independently; approx position is only an initial guess."""
    init = np.array(approx_position) if approx_position else np.zeros(3)
    results = []
    for obs in epochs:
        results.append(solve_epoch(obs, sp3, init))
    return results
