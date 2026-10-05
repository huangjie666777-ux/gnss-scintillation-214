import subprocess
import sys
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import app


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_synthetic.py")],
                   check=True, cwd=ROOT)


def test_position_endpoint():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/position",
                    files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)})
    f1.close()
    f2.close()
    assert r.status_code == 200
    body = r.json()
    assert body["n_epochs"] == 5 and body["n_ok"] == 5
    ep = body["epochs"][0]
    assert ep["status"] == "ok"
    assert abs(ep["geodetic"]["lat_deg"] - 30.0) < 1e-4
    assert abs(ep["receiver_clock_s"] - 1.5e-4) < 1e-6
    assert len(ep["residuals_m"]) == 8


def test_reject_garbage():
    client = TestClient(app)
    r = client.post("/position",
                    files={"rinex": ("a.rnx", b"garbage"),
                           "sp3": ("b.sp3", b"garbage")})
    assert r.status_code == 422
    assert "rinex" in r.json()["detail"]


def _post_tec(client, biases=None):
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    if biases is None:
        biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    r = client.post("/tec", files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data={"station_lat_deg": "30.0", "station_lon_deg": "114.0",
                          "station_height_m": "50.0", "biases": json.dumps(biases)})
    f1.close()
    f2.close()
    return r


def test_tec_endpoint():
    client = TestClient(app)
    r = _post_tec(client)
    assert r.status_code == 200
    body = r.json()
    sats = body["satellites"]
    assert len(sats) == 8
    g02 = sats["G02"]["records"]
    ok = [x for x in g02 if x["status"] == "ok"]
    assert len(ok) == 5
    for x in ok:
        assert x["arc_id"] == "G02-A1"
        assert 5.0 < x["stec_tecu"] < 30.0
        assert abs(x["vtec_tecu"]) <= abs(x["stec_tecu"]) + 1e-9
        assert -90.0 <= x["pierce_point"]["lat_deg"] <= 90.0
    # G01 is below the 10 deg elevation cutoff at this station
    g01 = sats["G01"]["records"]
    assert all(x["status"] == "excluded" and "elevation" in x["reason"] for x in g01)
    # G03: LLI at epoch 3 splits arc; trailing 2-sample arc is invalid
    g03 = sats["G03"]["records"]
    assert sats["G03"]["n_arcs"] == 2
    assert g03[3]["status"] == "invalid" and "loss-of-lock" in g03[3]["reason"]
    assert g03[4]["status"] == "invalid" and "fewer than 3" in g03[4]["reason"]
    assert g03[0]["arc_id"] == "G03-A1"


def test_tec_missing_bias_invalid():
    client = TestClient(app)
    biases = json.loads((ROOT / "examples" / "biases.json").read_text())
    del biases["G05"]
    r = _post_tec(client, biases)
    assert r.status_code == 200
    g05 = r.json()["satellites"]["G05"]["records"]
    assert all(x["status"] == "invalid" and "bias" in x["reason"] for x in g05)


def test_tec_rejects_bad_inputs():
    client = TestClient(app)
    f1 = open(ROOT / "examples" / "obs.rnx", "rb")
    f2 = open(ROOT / "examples" / "eph.sp3", "rb")
    r = client.post("/tec", files={"rinex": ("obs.rnx", f1), "sp3": ("eph.sp3", f2)},
                    data={"station_lat_deg": "95.0", "station_lon_deg": "114.0",
                          "station_height_m": "50.0", "biases": "{}"})
    f1.close()
    f2.close()
    assert r.status_code == 422
    assert "station_lat_deg" in r.json()["detail"]
