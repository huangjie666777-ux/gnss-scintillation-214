"""Generate synthetic RINEX 3.04 obs, SP3-c ephemeris and code biases.

Simulates 8 GPS satellites on circular orbits, a fixed receiver with a clock
offset, a smooth slant-TEC ionosphere, per-satellite P1-P2 code biases and
one LLI-flagged epoch (G03) to demonstrate arc splitting.
Writes examples/obs.rnx, examples/eph.sp3, examples/biases.json.
"""

import datetime as dt
import json
import math
import os

import numpy as np

GPS_EPOCH = dt.datetime(1980, 1, 6, tzinfo=dt.timezone.utc)
C = 299792458.0
OMEGA = 7.2921151467e-5

# receiver: lat 30N, lon 114E, h 50 m
LAT, LON, H = math.radians(30.0), math.radians(114.0), 50.0
A = 6378137.0
F = 1 / 298.257223563
E2 = F * (2 - F)
N = A / math.sqrt(1 - E2 * math.sin(LAT) ** 2)
RX = np.array([(N + H) * math.cos(LAT) * math.cos(LON),
               (N + H) * math.cos(LAT) * math.sin(LON),
               (N * (1 - E2) + H) * math.sin(LAT)])
RX_CLK = 1.5e-4  # seconds

ORB_R = 26560e3
N_SATS = 8
F1 = 1575.42e6
F2 = 1227.60e6
LAM1 = C / F1
LAM2 = C / F2


def stec_tecu(prn_idx: int, t: float) -> float:
    """Synthetic slant TEC in TECU (smooth in time, varies per satellite)."""
    return 15.0 + 8.0 * math.sin(2 * math.pi * t / 3600.0 + prn_idx * 1.3)


def dcb_ns(prn_idx: int) -> float:
    """Synthetic merged P1-P2 code bias in ns."""
    return 0.5 * (prn_idx + 1)


def sat_pos(prn_idx: int, t: float) -> np.ndarray:
    """Circular orbit position at t seconds since GPS epoch (inertial-ish frame
    rotated to ECEF at time t)."""
    raan = math.radians(prn_idx * 45.0)
    inc = math.radians(55.0)
    u0 = math.radians(prn_idx * 37.0 + 160.0)  # phase chosen so 4+ sats are visible
    n = 2 * math.pi / 43082.0  # ~half sidereal day
    u = u0 + n * t
    x_o = ORB_R * math.cos(u)
    y_o = ORB_R * math.sin(u)
    # rotate by inclination then RAAN, then Earth rotation to ECEF
    x1 = x_o
    y1 = y_o * math.cos(inc)
    z1 = y_o * math.sin(inc)
    x2 = x1 * math.cos(raan) - y1 * math.sin(raan)
    y2 = x1 * math.sin(raan) + y1 * math.cos(raan)
    theta = OMEGA * t
    x3 = x2 * math.cos(theta) + y2 * math.sin(theta)
    y3 = -x2 * math.sin(theta) + y2 * math.cos(theta)
    return np.array([x3, y3, z1])


def sat_clk(prn_idx: int, t: float) -> float:
    return 2e-5 * math.sin(2 * math.pi * t / 86400.0 + prn_idx) + prn_idx * 1e-7


def fmt_dt(t: dt.datetime) -> str:
    return (f"{t.year:5d}{t.month:3d}{t.day:3d}{t.hour:3d}{t.minute:3d}"
            f"{t.second + t.microsecond / 1e6:11.7f}")


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    t0 = dt.datetime(2024, 6, 1, 0, 10, 0, tzinfo=dt.timezone.utc)
    n_epochs = 5
    step = 30  # seconds between obs epochs

    # --- SP3: 900 s spacing, wide margin for 8-node windows ---
    sp3_start = t0 - dt.timedelta(seconds=3600)
    sp3_n = 12
    lines = []
    first = sp3_start
    lines.append(f"#cP{first.year:4d} {first.month:2d} {first.day:2d} "
                 f"{first.hour:2d} {first.minute:2d} {first.second:11.8f} "
                 f"{sp3_n:7d}    u+U  IGS14 FIT  AIUB")
    lines.append("## 0000      0.00000000    900.00000000   00000   0.0000000000000")
    lines.append("+   10   G01G02G03G04G05G06G07G08G09G10  0  0  0  0  0  0  0")
    lines.append("++          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0")
    lines.append("%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc")
    lines.append("%f  0.0000000  0.000000000  0.00000000000  0.000000000000000")
    lines.append("%i    0    0    0    0         0         0         0         0")
    lines.append("/* SYNTHETIC TEST EPHEMERIS - NOT REAL NAVIGATION DATA")
    for k in range(sp3_n):
        t = sp3_start + dt.timedelta(seconds=900 * k)
        ts = (t - GPS_EPOCH).total_seconds()
        lines.append(f"*  {t.year:4d} {t.month:2d} {t.day:2d} {t.hour:2d} "
                     f"{t.minute:2d} {t.second:11.8f}")
        for s in range(N_SATS):
            p = sat_pos(s, ts) / 1000.0
            c_us = sat_clk(s, ts) * 1e6
            lines.append(f"PG{s + 1:02d}{p[0]:14.6f}{p[1]:14.6f}{p[2]:14.6f}"
                         f"{c_us:14.6f}")
    lines.append("EOF")
    with open(os.path.join(outdir, "eph.sp3"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # --- code biases ---
    biases = {f"G{s + 1:02d}": dcb_ns(s) for s in range(N_SATS)}
    with open(os.path.join(outdir, "biases.json"), "w") as f:
        json.dump(biases, f, indent=2)

    # --- RINEX obs ---
    rng = np.random.default_rng(42)
    lines = []
    lines.append(f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE")
    lines.append(f"{'SYNTH':<20}{'TEST':<20}{'20240601 000000 UTC':<20}PGM / RUN BY / DATE")
    lines.append(f"{RX[0]:14.4f}{RX[1]:14.4f}{RX[2]:14.4f}{'':18}APPROX POSITION XYZ")
    lines.append(f"{0:6d}{'':54}RCV CLOCK OFFS APPL")
    lines.append(f"{'G':<1}{4:5d} {'C1C':<4}{'C2W':<4}{'L1C':<4}{'L2W':<4}"
                 f"{'':37}SYS / # / OBS TYPES")
    lines.append(f"{t0.year:6d}{t0.month:6d}{t0.day:6d}{t0.hour:6d}{t0.minute:6d}"
                 f"{t0.second:12.7f}GPS{'':9}TIME OF FIRST OBS")
    lines.append(f"{'':60}END OF HEADER")
    for e in range(n_epochs):
        t = t0 + dt.timedelta(seconds=step * e)
        ts = (t - GPS_EPOCH).total_seconds()
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{t.second:11.7f}  0{N_SATS:3d}{'':6}")
        for s in range(N_SATS):
            # transmit time iteration matching the solver
            pr_guess = 22000e3
            for _ in range(4):
                t_tx = ts - (pr_guess / C + sat_clk(s, ts - pr_guess / C))
                p = sat_pos(s, t_tx)
                # Earth rotation during flight
                tau = np.linalg.norm(p - RX) / C
                th = OMEGA * tau
                rot = np.array([[math.cos(th), math.sin(th), 0],
                                [-math.sin(th), math.cos(th), 0],
                                [0, 0, 1]])
                pr_guess = (np.linalg.norm(rot @ p - RX)
                            + C * RX_CLK - C * sat_clk(s, t_tx))
            # ionosphere: code delayed, phase advanced, by 40.3*STEC/f^2
            stec_si = stec_tecu(s, ts) * 1e16
            i1 = 40.3 * stec_si / F1**2
            i2 = 40.3 * stec_si / F2**2
            c1c = pr_guess + i1 + rng.normal(0, 0.3)
            c2w = pr_guess + i2 + C * dcb_ns(s) * 1e-9 + rng.normal(0, 0.05)
            n1 = 100000 + s * 137  # constant integer ambiguities (cycles)
            n2 = 80000 + s * 91
            l1c = (pr_guess - i1) / LAM1 + n1 + rng.normal(0, 0.002)
            l2w = (pr_guess - i2) / LAM2 + n2 + rng.normal(0, 0.002)
            # G03 loses lock on L1 at epoch index 3 -> arc break demo
            lli1 = 1 if (s == 2 and e == 3) else 0
            lines.append(f"G{s + 1:02d}"
                         f"{c1c:14.3f}  "
                         f"{c2w:14.3f}  "
                         f"{l1c:14.3f}{lli1} "
                         f"{l2w:14.3f}  ")
    with open(os.path.join(outdir, "obs.rnx"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote examples/obs.rnx, examples/eph.sp3, examples/biases.json")


if __name__ == "__main__":
    main()
