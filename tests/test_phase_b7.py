import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("phase_b7", ROOT / "scripts" / "phase_b7_sampling_tracking.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_rayleigh_separation_flags_same_frequency():
    r = mod.rayleigh_separation(365.0, 365.0, 1826)
    assert r["within_one_rayleigh"] is True
    assert r["frequency_separation_rayleigh"] == 0.0


def test_tracking_summary_counts_alias_windows_and_correlation():
    windows = []
    for sp, wp, alias in [(350.0, 352.0, True), (360.0, 362.0, True), (370.0, 373.0, True), (380.0, 384.0, True)]:
        windows.append({
            "signal": {"period_days": sp, "is_local_peak": True},
            "sampling_window": {"period_days": wp},
            "signal_vs_sampling_window": {"frequency_separation_rayleigh": 0.2, "within_one_rayleigh": alias},
        })
    r = mod.tracking_summary(windows)
    assert r["signal_local_peak_windows"] == 4
    assert r["alias_warning_windows"] == 4
    assert r["alias_fraction_among_local_peak_windows"] == 1.0
    assert r["signal_window_frequency_correlation"] > 0.99
