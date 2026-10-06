import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import app
from gnss_scint214.reflect import build_grid


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_scint.py")],
                   check=True, cwd=ROOT)


def _post(client, path=ROOT / "examples" / "obs_scint.rnx"):
    with open(path, "rb") as f:
        return client.post("/scintillation", files={"rinex": ("obs.rnx", f)})


def test_scintillation_endpoint():
    client = TestClient(app)
    r = _post(client)
    assert r.status_code == 200
    body = r.json()
    assert body["n_epochs"] == 15000
    assert body["interval_s"] == 0.02
    sats = body["satellites"]
    assert set(sats) == {"G01", "G02"}

    g01 = sats["G01"]
    assert g01["n_arcs"] == 1
    assert g01["arcs"][0]["n_samples"] == 15000
    wins = g01["windows"]
    assert len(wins) == 5
    # first minute: filter warm-up (window start < 60 s after arc start)
    assert wins[0]["status"] == "failed" and "warm-up" in wins[0]["reason"]
    for w in wins[1:]:
        assert w["status"] == "ok"
        assert w["arc_id"] == "G01-A1"
        assert w["n_samples"] == 3000
        assert 0.0 < w["s4"] < 1.0
        assert w["sigma_phi_rad"] > 100.0  # 0.5 Hz phase oscillation passes

    g02 = sats["G02"]
    # outage at 90 s and LLI at 150 s -> 3 arcs
    assert g02["n_arcs"] == 3
    recs = g02["records"]
    excl = [x for x in recs if x["status"] == "excluded"]
    assert len(excl) == 3
    assert any("missing observation" in x["reason"] for x in excl)
    assert any("loss-of-lock" in x["reason"] for x in excl)
    wins2 = {w["window_start"]: w for w in g02["windows"]}
    by_min = [wins2[k] for k in sorted(wins2)]
    # minute 1 (60-120 s): outage leaves it incomplete
    assert by_min[1]["status"] == "failed" and "gap" in by_min[1]["reason"]
    # minute 2 (120-180 s): LLI epoch excluded -> 2999 of 3000 samples
    assert by_min[2]["status"] == "failed" and "2999" in by_min[2]["reason"]
    # minute 3 (180-240 s): single arc but warm-up after the 150 s break
    assert by_min[3]["status"] == "failed" and "warm-up" in by_min[3]["reason"]
    # minute 4 (240-300 s): complete, single arc, warm-up satisfied
    assert by_min[4]["status"] == "ok" and by_min[4]["arc_id"] == "G02-A3"


def test_scintillation_rejects_bad_interval():
    client = TestClient(app)
    text = (ROOT / "examples" / "obs_scint.rnx").read_text()
    bad = text.replace("    0.0200", "    1.0000", 1)
    r = client.post("/scintillation", files={"rinex": ("obs.rnx", bad)})
    assert r.status_code == 422
    assert "INTERVAL" in r.json()["detail"]


def test_scintillation_rejects_missing_obs_types():
    client = TestClient(app)
    with open(ROOT / "examples" / "obs_reflect.rnx", "rb") as f:
        r = client.post("/scintillation", files={"rinex": ("obs.rnx", f)})
    assert r.status_code == 422
    assert "L1C" in r.json()["detail"]


def test_grid_includes_legal_upper_bound():
    # (2.5 - 0.5) / 0.1 == 19.999... in floating point; 2.5 must survive
    grid = build_grid(0.5, 2.5, 0.1)
    assert grid[-1] == pytest.approx(2.5)
    assert len(grid) == 21
    grid2 = build_grid(1.0, 8.0, 0.005)
    assert grid2[-1] == pytest.approx(8.0)
