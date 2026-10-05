"""Reflectometry endpoint, arc rules and parser/interp regression tests."""

import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import app
from gnss_scint214.interp import satellite_position
from gnss_scint214.rinex import parse_rinex
from gnss_scint214.sp3 import SatRecord


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_synthetic.py")],
                   check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_reflect.py")],
                   check=True, cwd=ROOT)


def _post(client, **over):
    f1 = open(ROOT / "examples" / "obs_reflect.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    data = {"station_lat_deg": "30.0", "station_lon_deg": "114.0",
            "station_height_m": "50.0", "antenna_height_m": "12.0",
            "h_min_m": "1.0", "h_max_m": "8.0", "h_step_m": "0.005"}
    data.update(over)
    r = client.post("/reflectometry",
                    files={"rinex": ("obs_reflect.rnx", f1),
                           "sp3": ("eph.sp3", f2)},
                    data=data)
    f1.close()
    f2.close()
    return r


def test_recovers_known_water_level():
    r = _post(TestClient(app))
    assert r.status_code == 200, r.json()
    body = r.json()
    ok_arcs = [a for sat in body["satellites"].values()
               for a in sat["arcs"] if a["status"] == "ok"]
    assert len(ok_arcs) >= 3
    for a in ok_arcs:
        assert abs(a["reflector_height_m"] - 3.6) < 0.05
        assert abs(a["water_level_m"] - 8.4) < 0.05
        assert a["amplitude"] > 0.0
        assert a["residual_rms"] > 0.0
        assert a["warning"] is None
        assert len(a["height_curve"]) == body["grid"]["n_points"]
        assert a["sample_source"] == "S1C"
        assert a["direction"] in ("rising", "setting")
        assert len(a["samples"]) == a["n_samples"]


def test_arc_splitting_and_filtering():
    body = _post(TestClient(app)).json()
    sats = body["satellites"]
    # G03: 3-epoch S1C outage splits the arc; leading 11-point arc fails
    g03 = sats["G03"]["arcs"]
    assert sats["G03"]["n_arcs"] == 2
    assert g03[0]["status"] == "failed" and "fewer than" in g03[0]["reason"]
    assert g03[1]["status"] == "ok"
    # outage epochs are reported, not dropped
    reasons = [r["reason"] for r in sats["G03"]["records"] if r["status"] != "ok"]
    assert any("S1C" in x for x in reasons)
    # G08: rising/setting reversal splits the arc
    assert sats["G08"]["n_arcs"] == 2
    dirs = {a["direction"] for a in sats["G08"]["arcs"]}
    assert dirs == {"rising", "setting"}
    # G01 is never in the 5-25 deg window: excluded with reason, still listed
    g01 = sats["G01"]["records"]
    assert len(g01) == 100
    assert all(r["status"] == "excluded" and "elevation" in r["reason"]
               for r in g01)


def test_boundary_optimum_warns():
    r = _post(TestClient(app), h_max_m="3.0")  # true H=3.6 m lies outside
    assert r.status_code == 200
    body = r.json()
    g03a2 = body["satellites"]["G03"]["arcs"][1]
    assert g03a2["status"] == "ok"
    assert g03a2["reflector_height_m"] == 3.0
    assert g03a2["warning"] and "boundary" in g03a2["warning"]


def test_grid_validation():
    client = TestClient(app)
    assert _post(client, h_min_m="8.0", h_max_m="1.0").status_code == 422
    assert _post(client, h_step_m="0.0").status_code == 422
    assert _post(client, h_step_m="-1.0").status_code == 422
    r = _post(client, h_min_m="0.1", h_max_m="10.0", h_step_m="0.001")
    assert r.status_code == 422 and "5001" in r.json()["detail"]
    assert _post(client, antenna_height_m="nan").status_code == 422


def test_s1c_only_header_accepted_but_position_rejected():
    client = TestClient(app)
    r = _post(client)
    assert r.status_code == 200  # S1C-only header is fine for reflectometry
    f1 = open(ROOT / "examples" / "obs_reflect.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/position",
                    files={"rinex": ("obs_reflect.rnx", f1),
                           "sp3": ("eph.sp3", f2)})
    f1.close()
    f2.close()
    assert r.status_code == 422 and "C1C" in r.json()["detail"]


def test_all_empty_satellite_record_kept():
    text = (ROOT / "examples" / "obs_reflect.rnx").read_text()
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("G05"):
            lines[i] = "G05" + " " * 16  # all-blank S1C field
            break
    obs = parse_rinex("\n".join(lines) + "\n")
    assert "G05" in obs.epochs[0].obs  # record survives, does not vanish
    assert obs.epochs[0].obs["G05"] == {}


def test_interp_does_not_cross_node_gap():
    rec = SatRecord()
    # 15 nodes at 900 s spacing with one missing node (k=11) mid-record
    rec.times = [900.0 * k for k in range(16) if k != 11]
    rec.pos = [(20000e3, 10000e3, 15000e3)] * len(rec.times)
    rec.clk = [0.0] * len(rec.times)
    rec.pos_valid = [True] * len(rec.times)
    rec.clk_valid = [True] * len(rec.times)
    pos, reason = satellite_position(rec, 900.0 * 10.5)
    assert pos is None and "gap" in reason
    pos, reason = satellite_position(rec, 900.0 * 1.5)
    assert pos is not None  # away from the gap it still interpolates
