import numpy as np
import pandas as pd

from src.modern import (
    _classic_lomb_basis,
    _classic_lomb_power_batch,
    _center_or_linear_detrend_batch,
    ar1_lomb_band_max_test,
)


def test_classic_lomb_basis_finds_synthetic_period():
    t = np.arange(0, 800, dtype=float)
    true_period = 80.0
    y = np.sin(2 * np.pi * t / true_period)
    periods = np.linspace(60, 100, 161)
    c, s, c2, s2 = _classic_lomb_basis(t, periods)
    yy = _center_or_linear_detrend_batch(y, t, "none")
    p = _classic_lomb_power_batch(yy, c, s, c2, s2)
    found = periods[int(np.nanargmax(p))]
    assert abs(found - true_period) <= 0.5


def test_linear_detrend_removes_linear_component():
    t = np.arange(100, dtype=float)
    y = 2.0 + 0.5 * t + np.sin(t / 4)
    z = _center_or_linear_detrend_batch(y, t, "linear")
    slope = np.polyfit(t, z, 1)[0]
    assert abs(slope) < 1e-12
    assert abs(np.mean(z)) < 1e-12


def test_ar1_lomb_band_max_test_is_deterministic_and_bounded():
    rng = np.random.default_rng(7)
    dates = pd.date_range("2000-01-01", periods=500, freq="D")
    y = np.exp(0.5 * np.sin(2 * np.pi * np.arange(500) / 90.0) + rng.normal(0, 0.2, 500))
    y[::13] = np.nan
    df = pd.DataFrame({"date": dates, "x": y})
    a = ar1_lomb_band_max_test(df, "x", 70, 110, simulations=19, nfreq=24, batch_size=5, seed=123)
    b = ar1_lomb_band_max_test(df, "x", 70, 110, simulations=19, nfreq=24, batch_size=5, seed=123)
    assert a["status"] == "ok"
    assert a["observed"] == b["observed"]
    assert a["p_value_band_max_ar1_lomb"] == b["p_value_band_max_ar1_lomb"]
    assert 0 < a["p_value_band_max_ar1_lomb"] <= 1
    assert a["null_simulations_used"] == 19
