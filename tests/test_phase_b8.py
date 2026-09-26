import importlib.util
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phase_b8", ROOT / "scripts" / "phase_b8_resolution_stress.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_rayleigh_period_resolution_shrinks_with_span():
    r5 = mod.rayleigh_period_context(470.0, 480.0, 5 * 365)
    r10 = mod.rayleigh_period_context(470.0, 480.0, 10 * 365)
    assert r10["approx_period_resolution_days"] < r5["approx_period_resolution_days"]
    assert r10["frequency_separation_rayleigh"] > r5["frequency_separation_rayleigh"]


def test_annual_windows_respects_requested_length():
    w = mod.annual_windows(pd.Timestamp("1973-01-01"), pd.Timestamp("1992-12-31"), years=10, step_years=1)
    assert len(w) == 11
    assert (w[0][1] - w[0][0]).days + 1 in (3652, 3653)


def test_scale_summary_reports_alias_fraction():
    rows = [
        {"signal": {"is_local_peak": True, "period_days": 470.0}, "sampling_window": {"period_days": 480.0},
         "signal_vs_sampling_window": {"within_one_rayleigh": True, "frequency_separation_rayleigh": 0.3,
                                        "absolute_period_difference_days": 10.0, "approx_period_resolution_days": 120.0}},
        {"signal": {"is_local_peak": True, "period_days": 470.0}, "sampling_window": {"period_days": 600.0},
         "signal_vs_sampling_window": {"within_one_rayleigh": False, "frequency_separation_rayleigh": 1.5,
                                        "absolute_period_difference_days": 130.0, "approx_period_resolution_days": 60.0}},
    ]
    s = mod._scale_summary(rows)
    assert s["local_peak_windows"] == 2
    assert s["alias_fraction"] == 0.5
