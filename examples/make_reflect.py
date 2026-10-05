"""Generate a synthetic S1C-only RINEX 3.04 file for reflectometry.

Simulates the same 8 GPS satellites as make_synthetic.py (reuses its orbit
model and examples/eph.sp3) plus a planar water surface at a known reflector
height below the antenna. SNR (S1C, dB-Hz) carries the interference
oscillation cos(4*pi*H*sin(elev)/lambda1). G04 has a 3-epoch S1C outage to
demonstrate arc splitting.

Ground truth: reflector height H_TRUE_M below the antenna phase center;
with ANTENNA_HEIGHT_M (Z, above gauge zero) the known water level is
WATER_LEVEL_M = ANTENNA_HEIGHT_M - H_TRUE_M.

Writes examples/obs_reflect.rnx. Run examples/make_synthetic.py first (or
afterwards; both write their own files).
"""

import datetime as dt
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from make_synthetic import GPS_EPOCH, N_SATS, sat_pos  # noqa: E402

from gnss_scint214.geodesy import elevation_deg, geodetic_to_ecef  # noqa: E402

C = 299792458.0
F1 = 1575.42e6
LAM1 = C / F1

STATION_LAT = 30.0
STATION_LON = 114.0
STATION_H = 50.0
ANTENNA_HEIGHT_M = 12.0   # Z: phase center above gauge zero
H_TRUE_M = 3.6            # reflector (water surface) height below antenna
WATER_LEVEL_M = ANTENNA_HEIGHT_M - H_TRUE_M  # 8.4 m

T0 = dt.datetime(2024, 6, 1, 0, 10, 0, tzinfo=dt.timezone.utc)
N_EPOCHS = 100
STEP_S = 60


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    station = geodetic_to_ecef(STATION_LAT, STATION_LON, STATION_H)
    rng = np.random.default_rng(7)

    lines = []
    lines.append(f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE")
    lines.append(f"{'SYNTH-REFL':<20}{'TEST':<20}{'20240601 000000 UTC':<20}"
                 "PGM / RUN BY / DATE")
    lines.append(f"{station[0]:14.4f}{station[1]:14.4f}{station[2]:14.4f}"
                 f"{'':18}APPROX POSITION XYZ")
    lines.append(f"{0:6d}{'':54}RCV CLOCK OFFS APPL")
    lines.append(f"{'G':<1}{1:5d} {'S1C':<4}{'':49}SYS / # / OBS TYPES")
    lines.append(f"{T0.year:6d}{T0.month:6d}{T0.day:6d}{T0.hour:6d}"
                 f"{T0.minute:6d}{T0.second:12.7f}GPS{'':9}TIME OF FIRST OBS")
    lines.append(f"{'':60}END OF HEADER")

    for e in range(N_EPOCHS):
        t = T0 + dt.timedelta(seconds=STEP_S * e)
        ts = (t - GPS_EPOCH).total_seconds()
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{t.second:11.7f}  0{N_SATS:3d}{'':6}")
        for s in range(N_SATS):
            pos = sat_pos(s, ts)
            el = elevation_deg(station, STATION_LAT, STATION_LON, pos)
            # interference oscillation + slow gain trend + noise
            phase = 4.0 * math.pi * H_TRUE_M * math.sin(math.radians(el)) / LAM1
            amp = 180.0 + 1.2 * el + 70.0 * math.cos(phase)
            amp += rng.normal(0.0, 4.0)
            s1c = 20.0 * math.log10(max(amp, 1.0))
            # G03: 3-epoch S1C outage (epochs 60-62) -> arc split demo
            if s == 2 and 60 <= e <= 62:
                lines.append(f"G{s + 1:02d}" + " " * 16)
            else:
                lines.append(f"G{s + 1:02d}{s1c:14.3f}  ")
    with open(os.path.join(outdir, "obs_reflect.rnx"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote examples/obs_reflect.rnx "
          f"(H={H_TRUE_M} m, Z={ANTENNA_HEIGHT_M} m, "
          f"water level={WATER_LEVEL_M} m)")


if __name__ == "__main__":
    main()
