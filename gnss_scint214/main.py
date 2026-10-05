"""FastAPI entry point: RINEX 3.04 + SP3-c positioning, GPS dual-band TEC
and single-frequency (S1C) GNSS-IR water-level reflectometry."""

from __future__ import annotations

import json
import math

from fastapi import FastAPI, File, Form, HTTPException, UploadFile

from .errors import GnssError
from .geodesy import geodetic_to_ecef
from .reflect import MAX_GRID_POINTS, SNR_OBS, build_grid, compute_reflectometry
from .rinex import RinexData, parse_rinex
from .solver import solve_all
from .sp3 import Sp3Data, parse_sp3
from .tec import REQUIRED_OBS, compute_tec

app = FastAPI(title="GNSS Positioning + Ionosphere TEC + Reflectometry",
              version="0.3.0")


async def _read_inputs(rinex: UploadFile, sp3: UploadFile) -> tuple[RinexData, Sp3Data]:
    if rinex.filename and rinex.filename.endswith((".gz", ".zip", ".Z")):
        raise HTTPException(400, "compressed RINEX not accepted; upload uncompressed file")
    try:
        rinex_text = (await rinex.read()).decode("ascii")
    except UnicodeDecodeError:
        raise HTTPException(400, "rinex: file is not plain ASCII text")
    try:
        sp3_text = (await sp3.read()).decode("ascii")
    except UnicodeDecodeError:
        raise HTTPException(400, "sp3: file is not plain ASCII text")
    try:
        obs = parse_rinex(rinex_text)
    except GnssError as e:
        raise HTTPException(422, f"RINEX rejected: {e}")
    try:
        eph = parse_sp3(sp3_text)
    except GnssError as e:
        raise HTTPException(422, f"SP3 rejected: {e}")
    return obs, eph


@app.post("/position")
async def position(rinex: UploadFile = File(...), sp3: UploadFile = File(...)):
    obs, eph = await _read_inputs(rinex, sp3)
    if "C1C" not in obs.obs_types_gps:
        raise HTTPException(422, "RINEX rejected: C1C pseudorange not declared "
                                 "in header (required for positioning)")
    results = solve_all(obs.epochs, eph, obs.approx_position)
    return {
        "n_epochs": len(results),
        "n_ok": sum(1 for r in results if r.status == "ok"),
        "notes": [
            "GPS time, GPS satellites, C1C only; no atmosphere/carrier corrections.",
            "Single-epoch independent solutions; approx position used as initial guess only.",
            "Accuracy boundary: standalone C1C pseudorange positioning, "
            "meter-level at best without ionosphere/troposphere modeling.",
        ],
        "epochs": [
            {
                "time": r.time,
                "status": r.status,
                "reason": r.reason,
                "ecef_m": r.ecef_m,
                "geodetic": r.geodetic,
                "receiver_clock_s": r.receiver_clock_s,
                "used_satellites": r.used_satellites,
                "excluded_satellites": r.excluded_satellites,
                "residuals_m": r.residuals_m,
                "rms_m": r.rms_m,
            }
            for r in results
        ],
    }


@app.post("/tec")
async def tec(
    rinex: UploadFile = File(...),
    sp3: UploadFile = File(...),
    station_lat_deg: float = Form(...),
    station_lon_deg: float = Form(...),
    station_height_m: float = Form(...),
    biases: str = Form(...),
):
    """Dual-frequency GPS ionosphere TEC from C1C/C2W/L1C/L2W.

    biases is a JSON object mapping PRN (e.g. G01) to the merged P1-P2
    code bias in nanoseconds; satellites without a bias are invalid.
    """
    obs, eph = await _read_inputs(rinex, sp3)

    for name, val, lo, hi in (
        ("station_lat_deg", station_lat_deg, -90.0, 90.0),
        ("station_lon_deg", station_lon_deg, -180.0, 180.0),
        ("station_height_m", station_height_m, -1000.0, 10000.0),
    ):
        if not math.isfinite(val) or not (lo <= val <= hi):
            raise HTTPException(422, f"{name}: value {val} not finite or out of "
                                     f"[{lo}, {hi}]")

    try:
        raw = json.loads(biases)
    except json.JSONDecodeError as e:
        raise HTTPException(422, f"biases: invalid JSON: {e}")
    if not isinstance(raw, dict):
        raise HTTPException(422, "biases: must be a JSON object {PRN: bias_ns}")
    biases_ns: dict[str, float] = {}
    for key, val in raw.items():
        if not isinstance(key, str) or not key.startswith("G"):
            raise HTTPException(422, f"biases: bad satellite id {key!r} (expect e.g. 'G01')")
        if not isinstance(val, (int, float)) or not math.isfinite(val):
            raise HTTPException(422, f"biases: non-finite bias for {key}")
        biases_ns[key.strip()] = float(val)

    missing_types = [t for t in REQUIRED_OBS if t not in obs.obs_types_gps]
    if missing_types:
        raise HTTPException(
            422, f"RINEX rejected: obs type(s) {missing_types} required for TEC "
                 "not declared in header")

    station_ecef = geodetic_to_ecef(station_lat_deg, station_lon_deg,
                                    station_height_m)
    result = compute_tec(obs, eph, station_ecef, biases_ns)
    result["station"] = {
        "lat_deg": station_lat_deg,
        "lon_deg": station_lon_deg,
        "height_m": station_height_m,
        "ecef_m": [float(v) for v in station_ecef],
    }
    result["notes"] = [
        "Code GF = (P2 - P1) - c*DCB(P1-P2); phase GF = L1*lam1 - L2*lam2 (m).",
        "Bias sign: subtracted term is the merged P1-P2 code bias in ns.",
        "Arcs break on missing obs, LLI!=0, gaps >120 s or phase jumps >1 m; "
        "each arc is leveled by its own mean(code - phase); arcs < 3 samples invalid.",
        "Slant TEC = GF / (40.3*(1/f2^2 - 1/f1^2)) in TECU; sign preserved.",
        "Thin-shell mapping at 6821 km geocentric radius; VTEC = STEC*cos(theta); "
        "pierce point reported in geocentric lat/lon.",
    ]
    return result


@app.post("/reflectometry")
async def reflectometry(
    rinex: UploadFile = File(...),
    sp3: UploadFile = File(...),
    station_lat_deg: float = Form(...),
    station_lon_deg: float = Form(...),
    station_height_m: float = Form(...),
    antenna_height_m: float = Form(...),
    h_min_m: float = Form(...),
    h_max_m: float = Form(...),
    h_step_m: float = Form(...),
):
    """GNSS-IR water level from L1 C/A SNR (S1C) alone.

    antenna_height_m is the antenna phase-center height Z above the gauge
    zero; h_min_m/h_max_m/h_step_m define the reflector-height search grid
    (positive, increasing bounds, positive step, at most 5001 grid points).
    """
    obs, eph = await _read_inputs(rinex, sp3)

    for name, val, lo, hi in (
        ("station_lat_deg", station_lat_deg, -90.0, 90.0),
        ("station_lon_deg", station_lon_deg, -180.0, 180.0),
        ("station_height_m", station_height_m, -1000.0, 10000.0),
        ("antenna_height_m", antenna_height_m, -1000.0, 10000.0),
    ):
        if not math.isfinite(val) or not (lo <= val <= hi):
            raise HTTPException(422, f"{name}: value {val} not finite or out of "
                                     f"[{lo}, {hi}]")
    for name, val in (("h_min_m", h_min_m), ("h_max_m", h_max_m),
                      ("h_step_m", h_step_m)):
        if not math.isfinite(val) or val <= 0.0:
            raise HTTPException(422, f"{name}: must be a positive finite value")
    if not h_min_m < h_max_m:
        raise HTTPException(422, "h_min_m must be strictly less than h_max_m")
    n_grid = len(build_grid(h_min_m, h_max_m, h_step_m))
    if n_grid > MAX_GRID_POINTS:
        raise HTTPException(422, f"height grid would have {n_grid} points, "
                                 f"more than the {MAX_GRID_POINTS} allowed; "
                                 "increase h_step_m or narrow the bounds")

    if SNR_OBS not in obs.obs_types_gps:
        raise HTTPException(422, f"RINEX rejected: {SNR_OBS} (L1 C/A SNR) not "
                                 "declared in header")

    station_ecef = geodetic_to_ecef(station_lat_deg, station_lon_deg,
                                    station_height_m)
    result = compute_reflectometry(obs, eph, station_ecef, station_lat_deg,
                                   station_lon_deg, antenna_height_m,
                                   h_min_m, h_max_m, h_step_m)
    result["station"] = {
        "lat_deg": station_lat_deg,
        "lon_deg": station_lon_deg,
        "height_m": station_height_m,
        "ecef_m": [float(v) for v in station_ecef],
    }
    result["notes"] = [
        "Static planar reflector; no tropospheric refraction correction.",
        "S1C (dB-Hz) -> linear amplitude 10^(S1C/20); quadratic trend in "
        "sin(elevation) removed per arc.",
        "Phase model 4*pi*H*sin(elevation)/lambda1; least-squares grid "
        "search, ties resolve to the smaller H.",
        "Arcs break on missing S1C, missing ephemeris, gaps >120 s, "
        "rising/setting reversals and any filtered epoch; arcs need "
        ">=12 points and >=5 deg elevation span.",
        "water_level_m = antenna_height_m - reflector_height_m.",
    ]
    return result


@app.get("/health")
async def health():
    return {"status": "ok"}
