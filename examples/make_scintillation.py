"""Generate a 50 Hz GPS L1 S1C/L1C RINEX 3.04 scintillation example."""

from __future__ import annotations

import datetime as dt
import math
import os

T0 = dt.datetime(2024, 6, 1, 0, 0, 0, tzinfo=dt.timezone.utc)
N_EPOCHS = 18000
INTERVAL_S = 0.02
N_SATS = 4


def main() -> None:
    outdir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(outdir, "obs_scint.rnx")
    lines = [
        f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE",
        "SYNTH-SCINT         TEST                20240601 000000 UTC" + " " * 20 +
        "PGM / RUN BY / DATE",
        f"{0.0:14.4f}{0.0:14.4f}{0.0:14.4f}{'':18}APPROX POSITION XYZ",
        f"{0:6d}{'':54}RCV CLOCK OFFS APPL",
        f"{'G':<1}{2:5d} {'S1C':<4}{'L1C':<4}{'':45}SYS / # / OBS TYPES",
        f"{INTERVAL_S:10.3f}{'':50}INTERVAL",
        f"{T0.year:6d}{T0.month:6d}{T0.day:6d}{T0.hour:6d}" +
        f"{T0.minute:6d}{T0.second:12.7f}{'':6}GPS{'':9}TIME OF FIRST OBS",
        " " * 60 + "END OF HEADER",
    ]

    for epoch in range(N_EPOCHS):
        t = T0 + dt.timedelta(seconds=INTERVAL_S * epoch)
        fractional_second = t.second + t.microsecond / 1_000_000.0
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{fractional_second:11.7f}  0{N_SATS:3d}{'':6}")
        time_s = epoch * INTERVAL_S
        for sat in range(1, N_SATS + 1):
            prn = f"G{sat:02d}"
            s1c = 42.0 + 0.8 * math.sin(2.0 * math.pi * 1.1 * time_s + sat)
            phase = 80.0 * time_s + 3.0 * math.sin(
                2.0 * math.pi * 0.8 * time_s + sat * 0.7)

            # G04 loses L1C at 130 s: the point is retained as an anomaly and
            # the phase filter restarts on the following arc.
            if prn == "G04" and epoch == 6500:
                lines.append(f"{prn}{s1c:14.3f}{' ':18}")
            else:
                lines.append(f"{prn}{s1c:14.3f}  {phase:14.6f}  ")

    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
