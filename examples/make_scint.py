"""Generate a synthetic 50 Hz S1C+L1C RINEX 3.04 file for scintillation.

Two GPS satellites, 300 s at 0.02 s interval (15000 epochs), GPS time,
normal epochs only. G01 is clean; G02 has a 2-epoch S1C/L1C outage at
t=90.00-90.02 s and one L1C LLI=1 epoch at t=150.00 s to demonstrate arc
splitting, warm-up rejection and window failure reasons.

S1C carries a 0.3 Hz scintillation-like intensity modulation; L1C is a
continuous carrier phase in cycles with a 0.5 Hz phase oscillation.

Writes examples/obs_scint.rnx.
"""

import datetime as dt
import math
import os

T0 = dt.datetime(2024, 6, 1, 0, 10, 0, tzinfo=dt.timezone.utc)
DURATION_S = 300.0
STEP_S = 0.02
N_EPOCHS = round(DURATION_S / STEP_S)  # 15000
SATS = ("G01", "G02")


def s1c_dbhz(t: float) -> float:
    return 45.0 + 3.0 * math.sin(2.0 * math.pi * 0.3 * t)


def l1c_cycles(t: float) -> float:
    return 1.0e6 + 8.4e4 * t + 100.0 * math.sin(2.0 * math.pi * 0.5 * t)


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    lines = []
    lines.append(f"{'3.04':<9}{'':11}{'O':<20}{'G':<20}RINEX VERSION / TYPE")
    lines.append(f"{'SYNTH-SCINT':<20}{'TEST':<20}{'20240601 000000 UTC':<20}"
                 "PGM / RUN BY / DATE")
    lines.append(f"{0:6d}{'':54}RCV CLOCK OFFS APPL")
    lines.append(f"{'G':<1}{2:5d} {'S1C':<4}{'L1C':<4}{'':45}"
                 "SYS / # / OBS TYPES")
    lines.append(f"{STEP_S:10.4f}{'':50}INTERVAL")
    lines.append(f"{T0.year:6d}{T0.month:6d}{T0.day:6d}{T0.hour:6d}"
                 f"{T0.minute:6d}{T0.second:12.7f}GPS{'':9}TIME OF FIRST OBS")
    lines.append(f"{'':60}END OF HEADER")

    for k in range(N_EPOCHS):
        t = T0 + dt.timedelta(milliseconds=20 * k)
        rel = 20 * k / 1000.0  # seconds since T0
        sec = t.second + t.microsecond / 1e6
        lines.append(f"> {t.year:4d} {t.month:02d} {t.day:02d} {t.hour:02d} "
                     f"{t.minute:02d}{sec:11.7f}  0{len(SATS):3d}{'':6}")
        for prn in SATS:
            # G02: 2-epoch outage at t=90.00/90.02, LLI=1 at t=150.00
            outage = prn == "G02" and 90.0 <= rel < 90.04
            lli = 1 if (prn == "G02" and abs(rel - 150.0) < 1e-9) else 0
            if outage:
                lines.append(f"{prn}" + " " * 32)
            else:
                lines.append(f"{prn}{s1c_dbhz(rel):14.3f}  "
                             f"{l1c_cycles(rel):14.5f}{lli} ")
    path = os.path.join(outdir, "obs_scint.rnx")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote examples/obs_scint.rnx "
          f"({N_EPOCHS} epochs, {len(SATS)} satellites, 0.02 s interval)")


if __name__ == "__main__":
    main()
