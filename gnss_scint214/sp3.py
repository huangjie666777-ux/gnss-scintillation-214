"""SP3-c precise ephemeris parser (GPS satellites, position/clock records)."""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field

from .errors import RejectedContentError, Sp3ParseError
from .geodesy import epoch_to_seconds

FILE = "sp3"


@dataclass
class SatRecord:
    # parallel arrays, sorted by time
    times: list[float] = field(default_factory=list)          # s since GPS epoch
    pos: list[tuple[float, float, float]] = field(default_factory=list)  # meters
    clk: list[float] = field(default_factory=list)            # seconds
    pos_valid: list[bool] = field(default_factory=list)
    clk_valid: list[bool] = field(default_factory=list)


@dataclass
class Sp3Data:
    sats: dict[str, SatRecord]


def parse_sp3(text: str) -> Sp3Data:
    lines = text.splitlines()
    if not lines:
        raise Sp3ParseError("empty file", FILE, 0)
    if not lines[0].startswith("#c"):
        raise RejectedContentError("only SP3-c format supported (first line must start '#c')",
                                   FILE, 1)
    try:
        n_epochs = int(lines[0][32:39])
    except ValueError:
        raise Sp3ParseError("cannot parse epoch count in header", FILE, 1)

    sats: dict[str, SatRecord] = {}
    epoch_count = 0
    cur_time: float | None = None
    cur_ln = 0
    i = 1
    # skip header until first epoch record
    while i < len(lines) and not lines[i].startswith("* "):
        i += 1
    if i >= len(lines):
        raise Sp3ParseError("no epoch records found (file truncated?)", FILE, len(lines))

    while i < len(lines):
        line = lines[i]
        ln = i + 1
        if line.startswith("EOF"):
            break
        if line.startswith("* "):
            try:
                parts = line.split()
                y, mo, d, h, mi = (int(parts[k]) for k in range(1, 6))
                sec = float(parts[6])
                whole = int(sec)
                t = dt.datetime(y, mo, d, h, mi, whole,
                                int(round((sec - whole) * 1e6)), tzinfo=dt.timezone.utc)
            except (ValueError, IndexError) as e:
                raise Sp3ParseError(f"invalid epoch line: {e}", FILE, ln)
            cur_time = epoch_to_seconds(t)
            cur_ln = ln
            epoch_count += 1
            i += 1
            continue
        if line.startswith("P"):
            if cur_time is None:
                raise Sp3ParseError("position record before any epoch", FILE, ln)
            if len(line) < 60:
                raise Sp3ParseError("position record truncated", FILE, ln)
            prn = line[1:4].strip()
            try:
                x = float(line[4:18]); y = float(line[18:32]); z = float(line[32:46])
                clk = float(line[46:60])
            except ValueError:
                raise Sp3ParseError("cannot parse position/clock values", FILE, ln)
            for v in (x, y, z, clk):
                if v != v or v in (float("inf"), float("-inf")):
                    raise Sp3ParseError("non-finite value in record", FILE, ln)
            if not prn.startswith("G"):
                i += 1
                continue  # GPS only
            rec = sats.setdefault(prn, SatRecord())
            # km -> m; zero coordinates mark invalid position
            pos_ok = not (x == 0.0 and y == 0.0 and z == 0.0)
            # us -> s; 999999.999999 marks missing clock
            clk_ok = abs(clk) < 999999.0
            rec.times.append(cur_time)
            rec.pos.append((x * 1000.0, y * 1000.0, z * 1000.0))
            rec.clk.append(clk * 1e-6)
            rec.pos_valid.append(pos_ok)
            rec.clk_valid.append(clk_ok)
            i += 1
            continue
        if line.startswith(("V", "EP", "EV", "E")):
            # velocity / correlation records not used; skip silently
            i += 1
            continue
        raise Sp3ParseError(f"unexpected record {line[:2]!r}", FILE, ln)

    if epoch_count != n_epochs:
        raise Sp3ParseError(
            f"epoch count mismatch: header says {n_epochs}, found {epoch_count}",
            FILE, cur_ln)
    if not sats:
        raise Sp3ParseError("no GPS satellite records found", FILE, len(lines))
    return Sp3Data(sats=sats)
