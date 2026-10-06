"""RINEX 3.04 observation file parser (GPS, normal epochs only).

Parses every declared GPS observation type; each value carries its LLI
(loss-of-lock indicator). C1C pseudoranges feed positioning; C1C/C2W/L1C/L2W
feed the dual-frequency TEC pipeline; S1C (L1 C/A SNR) feeds single-frequency
reflectometry. Any subset of these may be declared (e.g. an S1C-only header);
each endpoint checks for the observation types it needs.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field

from .errors import RejectedContentError, RinexParseError

FILE = "rinex"
MAX_EPOCHS = 100
SCINTILLATION_MAX_EPOCHS = 18000
SCINTILLATION_MAX_SATELLITES = 4
SCINTILLATION_INTERVAL_S = 0.02
TIME_TOLERANCE_S = 1e-6


@dataclass
class EpochObs:
    time: dt.datetime
    # prn -> pseudorange in meters (C1C); missing sats absent
    pseudoranges: dict[str, float] = field(default_factory=dict)
    # prn -> obs type -> (value, lli); missing values absent
    obs: dict[str, dict[str, tuple[float, int]]] = field(default_factory=dict)


@dataclass
class RinexData:
    epochs: list[EpochObs]
    approx_position: tuple[float, float, float] | None
    obs_types_gps: list[str]


def _parse_float_field(text: str, file: str, line: int) -> float | None:
    """Parse a fixed-width RINEX float field; blank or zero -> None (missing)."""
    s = text.strip()
    if not s:
        return None
    try:
        v = float(s.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise RinexParseError(f"invalid numeric field {s!r}", file, line)
    if not (v == v) or v in (float("inf"), float("-inf")):
        raise RinexParseError("non-finite observation value", file, line)
    if v == 0.0:
        return None
    return v


def _parse_header(lines: list[str]) -> tuple[dict, int]:
    header: dict[str, list[tuple[str, int]]] = {}
    for i, raw in enumerate(lines, start=1):
        line = raw.rstrip("\n")
        if len(line) < 60:
            raise RinexParseError("header line shorter than 60 chars (truncated?)", FILE, i)
        label = line[60:].strip()
        header.setdefault(label, []).append((line[:60], i))
        if label == "END OF HEADER":
            return header, i
    raise RinexParseError("END OF HEADER not found (file truncated)", FILE, len(lines))


def parse_rinex(text: str, strict_scintillation: bool = False) -> RinexData:
    lines = text.splitlines()
    if not lines:
        raise RinexParseError("empty file", FILE, 0)
    header, n_header = _parse_header(lines)

    # --- RINEX VERSION / TYPE ---
    ver_line, ver_ln = header.get("RINEX VERSION / TYPE", [("", 0)])[0]
    try:
        version = float(ver_line[0:9])
    except ValueError:
        raise RinexParseError("cannot parse RINEX version", FILE, ver_ln)
    if abs(version - 3.04) > 1e-9:
        raise RejectedContentError(f"only RINEX 3.04 supported, got {version}", FILE, ver_ln)
    if ver_line[20:21] != "O":
        raise RejectedContentError("only observation files (type O) supported", FILE, ver_ln)
    if ver_line[40:41] not in ("G", " "):
        raise RejectedContentError("only GPS satellite system files supported", FILE, ver_ln)

    # --- TIME OF FIRST OBS must be GPS ---
    tobs = header.get("TIME OF FIRST OBS")
    if tobs is not None:
        tline, tln = tobs[0]
        if tline[48:60].strip() != "GPS":
            raise RejectedContentError("only GPS time system supported", FILE, tln)
    elif strict_scintillation:
        raise RejectedContentError("TIME OF FIRST OBS GPS header is required", FILE, 0)

    interval_s = None
    if "INTERVAL" in header:
        iline, iln = header["INTERVAL"][0]
        try:
            interval_s = float(iline[0:10])
        except ValueError:
            raise RinexParseError("bad INTERVAL value", FILE, iln)
    if strict_scintillation:
        if interval_s is None:
            raise RejectedContentError("INTERVAL header is required", FILE, 0)
        if abs(interval_s - SCINTILLATION_INTERVAL_S) > TIME_TOLERANCE_S:
            raise RejectedContentError(
                f"INTERVAL must be {SCINTILLATION_INTERVAL_S:.2f} s for "
                f"scintillation monitoring, got {interval_s}", FILE, iln)

    # --- reject receiver clock offset pre-applied ---
    for key in ("RCV CLOCK OFFS APPL", "LEAP SECONDS"):
        if key in header and key == "RCV CLOCK OFFS APPL":
            vline, vln = header[key][0]
            try:
                if int(vline[0:6]) != 0:
                    raise RejectedContentError(
                        "receiver clock offsets pre-applied not supported", FILE, vln)
            except ValueError:
                raise RinexParseError("bad RCV CLOCK OFFS APPL value", FILE, vln)

    # --- observation types (with continuation lines) ---
    sys_obs = header.get("SYS / # / OBS TYPES")
    if not sys_obs:
        raise RinexParseError("SYS / # / OBS TYPES header missing", FILE, 0)
    obs_types: list[str] | None = None
    for content, ln in sys_obs:
        sys_id = content[0:1]
        if sys_id != "G":
            continue
        try:
            n_types = int(content[1:6])
        except ValueError:
            raise RinexParseError("bad GPS obs type count", FILE, ln)
        tokens = content[7:60].split()
        need = n_types - len(tokens)
        obs_types = tokens
        # continuation lines: 13 types per line, 4 chars each
        cont_idx = sys_obs.index((content, ln)) + 1
        while need > 0:
            if cont_idx >= len(sys_obs):
                raise RinexParseError("obs types continuation truncated", FILE, ln)
            ccontent, cln = sys_obs[cont_idx]
            if ccontent[0:1] not in (" ", "G"):
                raise RinexParseError("obs types continuation truncated", FILE, cln)
            ctokens = ccontent[7:60].split()
            obs_types.extend(ctokens)
            need -= len(ctokens)
            cont_idx += 1
        obs_types = obs_types[:n_types]
        break
    if obs_types is None:
        raise RejectedContentError("no GPS (G) observation types in header", FILE, 0)

    # --- approximate position (initial guess only) ---
    approx = None
    if "APPROX POSITION XYZ" in header:
        aline, aln = header["APPROX POSITION XYZ"][0]
        try:
            approx = tuple(float(aline[i:i + 14]) for i in (0, 14, 28))
        except ValueError:
            raise RinexParseError("bad APPROX POSITION XYZ", FILE, aln)

    # --- epoch records ---
    epochs: list[EpochObs] = []
    i = n_header  # 0-based index of first data line
    while i < len(lines):
        line = lines[i]
        ln = i + 1
        if not line.strip():
            i += 1
            continue
        if line[0] != ">":
            raise RinexParseError("expected epoch record starting with '>'", FILE, ln)
        if len(line) < 41:
            raise RinexParseError("epoch header line truncated", FILE, ln)
        try:
            year = int(line[2:6]); month = int(line[7:9]); day = int(line[10:12])
            hour = int(line[13:15]); minute = int(line[16:18]); sec = float(line[18:29])
            flag = int(line[29:32]); n_sat = int(line[32:35])
        except ValueError:
            raise RinexParseError("cannot parse epoch header fields", FILE, ln)
        if flag != 0:
            raise RejectedContentError(
                f"epoch flag {flag} (event/special record) not supported", FILE, ln)
        if not (0 <= sec < 61):
            raise RinexParseError("invalid epoch seconds", FILE, ln)
        try:
            whole = int(sec)
            t = dt.datetime(year, month, day, hour, minute, whole,
                            int(round((sec - whole) * 1e6)), tzinfo=dt.timezone.utc)
        except ValueError as e:
            raise RinexParseError(f"invalid epoch date: {e}", FILE, ln)
        max_sats = (SCINTILLATION_MAX_SATELLITES
                    if strict_scintillation else 64)
        if n_sat < 0 or n_sat > max_sats:
            raise RinexParseError(f"implausible satellite count {n_sat}", FILE, ln)
        i += 1
        obs = EpochObs(time=t)
        for _ in range(n_sat):
            if i >= len(lines):
                raise RinexParseError("file truncated inside epoch observations", FILE, ln)
            oline = lines[i]
            oln = i + 1
            if len(oline) < 3:
                raise RinexParseError("observation line truncated (no satellite id)", FILE, oln)
            prn = oline[0:3]
            if prn[0] != "G":
                i += 1
                continue  # non-GPS satellites ignored per spec
            need_len = 3 + 16 * len(obs_types)
            if len(oline) < need_len:
                raise RinexParseError(
                    f"observation line truncated: need {need_len} chars for "
                    f"{len(obs_types)} obs types", FILE, oln)
            sat_obs: dict[str, tuple[float, int]] = {}
            for idx, otype in enumerate(obs_types):
                field_start = 3 + 16 * idx
                field = oline[field_start:field_start + 16]
                val = _parse_float_field(field[:14], FILE, oln)
                if val is None:
                    continue
                lli_raw = field[14:15].strip()
                if lli_raw and lli_raw not in "01234567":
                    raise RinexParseError(
                        f"invalid LLI flag {lli_raw!r} for {otype}", FILE, oln)
                sat_obs[otype] = (val, int(lli_raw) if lli_raw else 0)
            # keep the record even when every field is blank/zero so the
            # satellite does not silently vanish from downstream reports
            obs.obs[prn.strip()] = sat_obs
            if "C1C" in sat_obs:
                obs.pseudoranges[prn.strip()] = sat_obs["C1C"][0]
            i += 1
        epochs.append(obs)
        epoch_limit = (SCINTILLATION_MAX_EPOCHS
                       if strict_scintillation else MAX_EPOCHS)
        if len(epochs) > epoch_limit:
            raise RejectedContentError(
                f"more than {epoch_limit} epochs not supported", FILE, ln)
    if not epochs:
        raise RinexParseError("no epochs found", FILE, len(lines))
    if strict_scintillation:
        all_prns = {prn for epoch in epochs for prn in epoch.obs}
        if len(all_prns) > SCINTILLATION_MAX_SATELLITES:
            raise RejectedContentError(
                f"more than {SCINTILLATION_MAX_SATELLITES} GPS satellites "
                f"not supported ({len(all_prns)} found)", FILE, n_header + 1)
        nominal_us = int(round(SCINTILLATION_INTERVAL_S * 1e6))
        tolerance_us = int(round(TIME_TOLERANCE_S * 1e6))
        for prev, cur in zip(epochs, epochs[1:]):
            delta = cur.time - prev.time
            delta_us = ((delta.days * 86400 + delta.seconds) * 1_000_000
                        + delta.microseconds)
            if delta_us <= 0:
                raise RejectedContentError("epoch times must be strictly increasing",
                                           FILE, 0)
            steps = int(round(delta_us / nominal_us))
            if steps < 1 or abs(delta_us - steps * nominal_us) > tolerance_us:
                raise RejectedContentError(
                    f"epoch time difference {(cur.time - prev.time).total_seconds():.6f} s "
                    f"is not an integer multiple of {SCINTILLATION_INTERVAL_S:.2f} s "
                    "within 1 us", FILE, 0)
    return RinexData(epochs=epochs, approx_position=approx, obs_types_gps=obs_types)
