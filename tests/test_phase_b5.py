import importlib.util
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phase_b5", ROOT / "scripts" / "phase_b5_periodicity.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_spectrum_peak_morphology_detects_internal_peak():
    q = pd.DataFrame({
        "period_days": [300, 330, 350, 375, 400, 430, 500],
        "power": [0.1, 0.2, 0.3, 0.9, 0.25, 0.2, 0.1],
    })
    r = mod.spectrum_peak_morphology(q, 330, 400)
    assert r["period_days"] == 375.0
    assert r["is_local_peak"] is True
    assert r["at_band_edge"] is False
    assert r["prominence"] > 0


def test_spectrum_peak_morphology_flags_boundary_maximum():
    q = pd.DataFrame({
        "period_days": [300, 330, 340, 350, 360, 370, 380, 390, 400, 430],
        "power": [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05],
    })
    r = mod.spectrum_peak_morphology(q, 330, 400)
    assert r["period_days"] == 330.0
    assert r["at_band_edge"] is True
    assert r["is_local_peak"] is False
