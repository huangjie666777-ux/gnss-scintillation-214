import datetime as dt
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gnss_scint214.errors import RejectedContentError, RinexParseError, Sp3ParseError
from gnss_scint214.geodesy import ecef_to_geodetic, epoch_to_seconds, geodetic_to_ecef
from gnss_scint214.rinex import parse_rinex
from gnss_scint214.solver import solve_all
from gnss_scint214.sp3 import parse_sp3
from gnss_scint214.tec import compute_tec
import json


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_synthetic.py")],
                   check=True, cwd=ROOT)


@pytest.fixture(scope="module")
def solved():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    return solve_all(obs.epochs, eph, obs.approx_position)


def test_position_accuracy(solved):
    true = np.array([-2248562.1597, 5050353.2992, 3170398.7354])
    for r in solved:
        assert r.status == "ok", r.reason
        assert np.linalg.norm(np.array(r.ecef_m) - true) < 5.0  # iono in C1C
        assert abs(r.receiver_clock_s - 1.5e-4) < 1e-6
        assert r.rms_m < 2.0
        assert len(r.used_satellites) == 8


def test_geodetic(solved):
    lat, lon, h = (solved[0].geodetic[k] for k in ("lat_deg", "lon_deg", "height_m"))
    assert abs(lat - 30.0) < 1e-4
    assert abs(lon - 114.0) < 1e-4
    assert abs(h - 50.0) < 2.0


def test_ecef_to_geodetic_roundtrip():
    lat, lon, h = ecef_to_geodetic(np.array([6378137.0, 0.0, 0.0]))
    assert abs(lat) < 1e-9 and abs(lon) < 1e-9 and abs(h) < 1e-6


def test_reject_event_flag():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    bad = text.replace("  0  8", "  2  8", 1)
    with pytest.raises(RejectedContentError):
        parse_rinex(bad)


@pytest.fixture(scope="module")
def tec_result():
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    station = geodetic_to_ecef(30.0, 114.0, 50.0)
    return compute_tec(obs, eph, station, biases)


def test_tec_recovers_synthetic_model(tec_result):
    import examples.make_synthetic as synth
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    t0 = epoch_to_seconds(obs.epochs[0].time)
    for prn, sat in tec_result["satellites"].items():
        idx = int(prn[1:]) - 1
        for rec, ep in zip(sat["records"], obs.epochs):
            if rec["status"] != "ok":
                continue
            t = epoch_to_seconds(ep.time)
            expected = synth.stec_tecu(idx, t)
            assert abs(rec["stec_tecu"] - expected) < 2.5, (prn, rec["time"])
            assert rec["vtec_tecu"] <= rec["stec_tecu"] + 1e-9


def test_tec_arc_rules(tec_result):
    g03 = tec_result["satellites"]["G03"]
    assert g03["n_arcs"] == 2  # LLI at epoch 3 splits the arc
    statuses = [r["status"] for r in g03["records"]]
    assert statuses[:3] == ["ok", "ok", "ok"]
    assert statuses[3:] == ["invalid", "invalid"]


def test_geodetic_to_ecef_roundtrip():
    xyz = geodetic_to_ecef(30.0, 114.0, 50.0)
    lat, lon, h = ecef_to_geodetic(xyz)
    assert abs(lat - 30.0) < 1e-9
    assert abs(lon - 114.0) < 1e-9
    assert abs(h - 50.0) < 1e-6


def test_reject_truncated_epoch():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    lines = text.splitlines()
    cut = "\n".join(lines[:9]) + "\n"  # header + epoch line, no obs
    with pytest.raises(RinexParseError) as e:
        parse_rinex(cut)
    assert "truncated" in str(e.value)


def test_reject_rcv_clock_applied():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    bad = text.replace("     0 ", "     1 ", 1)
    with pytest.raises(RejectedContentError):
        parse_rinex(bad)


def test_reject_too_many_epochs():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    header, _, body = text.partition("END OF HEADER")
    one = body.strip().split("> ")[1]
    epoch_block = "> " + one
    big = header + "END OF HEADER" + "\n" + epoch_block * 101
    with pytest.raises(RejectedContentError):
        parse_rinex(big)


def test_sp3_epoch_count_mismatch():
    text = (ROOT / "examples" / "eph.sp3").read_text()
    bad = text.replace("      12 ", "      13 ", 1)
    with pytest.raises(Sp3ParseError):
        parse_sp3(bad)


def test_sp3_zero_coordinates_excluded():
    text = (ROOT / "examples" / "eph.sp3").read_text()
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("PG01"):
            lines[i] = "PG01" + f"{0.0:14.6f}" * 3 + f"{1.0:14.6f}"
            break
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3("\n".join(lines) + "\n")
    res = solve_all(obs.epochs, eph, obs.approx_position)
    g01 = [e for e in res[0].excluded_satellites if e["prn"] == "G01"]
    assert g01 and "coordinates" in g01[0]["reason"]
    assert res[0].status == "ok"  # 7 satellites remain


def test_sp3_missing_clock_excluded():
    text = (ROOT / "examples" / "eph.sp3").read_text()
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("PG02"):
            lines[i] = ln[:46] + f"{999999.999999:14.6f}" + ln[60:]
    obs = parse_rinex((ROOT / "examples" / "obs.rnx").read_text())
    eph = parse_sp3("\n".join(lines) + "\n")
    res = solve_all(obs.epochs, eph, obs.approx_position)
    g02 = [e for e in res[0].excluded_satellites if e["prn"] == "G02"]
    assert g02 and "clock" in g02[0]["reason"]


def test_insufficient_satellites():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    lines = text.splitlines()
    out = []
    skip = 0
    for ln in lines:
        if skip and (ln.startswith("G") or ln.startswith(">")):
            skip -= 1
            continue
        if ln.startswith("> "):
            ln = ln[:32] + "  3" + ln[35:]
            skip = 5
        out.append(ln)
    obs = parse_rinex("\n".join(out) + "\n")
    eph = parse_sp3((ROOT / "examples" / "eph.sp3").read_text())
    res = solve_all(obs.epochs, eph, obs.approx_position)
    assert all(r.status == "failed" and "usable satellites" in r.reason for r in res)
