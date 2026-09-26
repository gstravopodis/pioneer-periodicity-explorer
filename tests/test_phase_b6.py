import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phase_b6", ROOT / "scripts" / "phase_b6_periodicity_atlas.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_alias_warning_when_signal_and_window_are_same_frequency():
    r = mod.nearest_window_peak_alias(365.0, [{"period_days": 365.0, "power": 1.0, "prominence": 1.0}], 3650)
    assert r["within_one_rayleigh"] is True
    assert r["frequency_separation_rayleigh"] == 0.0


def test_alias_not_flagged_when_well_separated_in_frequency():
    r = mod.nearest_window_peak_alias(500.0, [{"period_days": 365.0, "power": 1.0, "prominence": 1.0}], 3650)
    assert r["within_one_rayleigh"] is False
    assert r["frequency_separation_rayleigh"] > 1.0
