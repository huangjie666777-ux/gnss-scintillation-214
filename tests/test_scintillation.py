import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from main import app
from gnss_scint214.reflect import build_grid
from gnss_scint214.rinex import parse_rinex
from gnss_scint214.scintillation import compute_scintillation


@pytest.fixture(scope="module", autouse=True)
def synthetic():
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_scintillation.py")],
                   check=True, cwd=ROOT)


def test_scintillation_windows_and_arc_reset():
    path = ROOT / "examples" / "obs_scint.rnx"
    obs = parse_rinex(path.read_text(), strict_scintillation=True)
    assert len(obs.epochs) == 18000
    result = compute_scintillation(obs)
    g01 = result["satellites"]["G01"]["windows"]
    assert [w["status"] for w in g01] == [
        "failed", "ok", "ok", "ok", "ok", "ok"]
    assert "warm-up" in g01[0]["reason"]
    for window in g01[1:]:
        assert window["n_samples"] == 3000
        assert window["arc_id"] == "G01-A1"
        assert window["s4"] > 0.0
        assert window["sigma_phi_rad"] > 0.0

    g04 = result["satellites"]["G04"]
    assert g04["n_arcs"] == 2
    assert g04["windows"][2]["status"] == "failed"
    assert "incomplete" in g04["windows"][2]["reason"]
    assert g04["windows"][3]["status"] == "failed"
    assert "warm-up" in g04["windows"][3]["reason"]
    assert g04["windows"][4]["status"] == "ok"
    assert g04["windows"][5]["status"] == "ok"
    bad = [q for q in g04["quality"] if q["status"] == "excluded"]
    assert len(bad) == 1 and "L1C" in bad[0]["reason"]


def test_scintillation_endpoint():
    client = TestClient(app)
    with open(ROOT / "examples" / "obs_scint.rnx", "rb") as f:
        response = client.post("/scintillation",
                               files={"rinex": ("obs_scint.rnx", f)})
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["window_samples"] == 3000
    assert len(body["satellites"]) == 4
    assert any("noise correction" in note for note in body["notes"])


def test_legacy_epoch_limit_remains_100():
    text = (ROOT / "examples" / "obs.rnx").read_text()
    header, _, body = text.partition("END OF HEADER")
    block = "\n> " + body.strip().split("> ")[1]
    with pytest.raises(Exception):
        parse_rinex(header + "END OF HEADER\n" + block * 101)


def test_strict_interval_must_be_0_02_seconds():
    text = (ROOT / "examples" / "obs_scint.rnx").read_text()
    bad = text.replace("     0.020", "     0.025", 1)
    with pytest.raises(Exception, match="INTERVAL"):
        parse_rinex(bad, strict_scintillation=True)


def test_strict_mode_rejects_more_than_four_satellites_per_epoch():
    header = (
        "3.04                O                   G                   RINEX VERSION / TYPE\n"
        "        0.0000        0.0000        0.0000                  APPROX POSITION XYZ\n"
        "     0                                                      RCV CLOCK OFFS APPL\n"
        "G    2 S1C L1C                                              SYS / # / OBS TYPES\n"
        "     0.020                                                  INTERVAL\n"
        "  2024     6     1     0     0   0.0000000      GPS         TIME OF FIRST OBS\n"
        "                                                            END OF HEADER\n"
    )
    epoch = (
        "> 2024 06 01 00 00  0.0000000  0  5      \n"
        + "".join(f"G{i:02d}{42.0:14.3f}  {1.0:14.6f}  \n" for i in range(1, 6))
    )
    with pytest.raises(Exception, match="satellite count"):
        parse_rinex(header + epoch, strict_scintillation=True)


def test_reflect_grid_includes_exact_upper_bound():
    grid = build_grid(1.0, 8.0, 0.005)
    assert len(grid) == 1401
    assert abs(grid[-1] - 8.0) < 1e-12
