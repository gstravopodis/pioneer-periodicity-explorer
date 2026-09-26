import importlib.util
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phase_b4", ROOT / "scripts" / "phase_b4_localization.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_peak_morphology_flags_band_edge_and_local_peak():
    q = pd.DataFrame({
        "lag_bin_start": [320, 330, 340, 350, 360, 370, 380, 390, 400],
        "lag_bin_end":   [329, 339, 349, 359, 369, 379, 389, 399, 409],
        "lag_bin_center":[324.5,334.5,344.5,354.5,364.5,374.5,384.5,394.5,404.5],
        "pairs": [1000]*9,
        "dcf": [0.1,0.5,0.4,0.3,0.2,0.15,0.12,0.1,0.08],
    })
    r = mod.peak_morphology(q, 330, 400)
    assert r["lag_calendar_days"] == 334.5
    assert r["at_band_edge"] is True
    assert r["is_local_peak"] is True
    assert r["prominence"] > 0
