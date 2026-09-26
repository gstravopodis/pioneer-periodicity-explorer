import numpy as np
import pandas as pd

from src.modern import _calendar_windows, ar1_lomb_intermit_joint_test


def test_calendar_windows_use_calendar_years():
    w = _calendar_windows(pd.Timestamp("2000-03-01"), pd.Timestamp("2010-03-01"), window_years=4, step_years=2)
    assert len(w) == 4
    assert str(w[0][0].date()) == "2000-03-01"
    assert str(w[0][1].date()) == "2004-02-29"
    assert str(w[-1][0].date()) == "2006-03-01"


def test_intermittent_joint_test_recovers_same_window_period_and_is_deterministic():
    rng = np.random.default_rng(4)
    dates = pd.date_range("2000-01-01", "2015-12-31", freq="D")
    t = np.arange(len(dates), dtype=float)
    # A shared 180 d oscillation is confined mostly to an early 8-year block.
    env = ((dates >= pd.Timestamp("2001-01-01")) & (dates <= pd.Timestamp("2008-12-31"))).astype(float)
    sig = 1.2 * env * np.sin(2 * np.pi * t / 180.0)
    y10 = np.exp(sig + rng.normal(0, 0.25, len(t)))
    y11 = np.exp(sig + rng.normal(0, 0.25, len(t)))
    y10[::17] = np.nan; y11[::19] = np.nan
    d10 = pd.DataFrame({"date": dates, "x": y10})
    d11 = pd.DataFrame({"date": dates, "x": y11})
    a = ar1_lomb_intermit_joint_test(d10, d11, "x", 160, 200, simulations=19,
                                    nfreq=32, window_years=8, step_years=2, batch_size=5, seed=77)
    b = ar1_lomb_intermit_joint_test(d10, d11, "x", 160, 200, simulations=19,
                                    nfreq=32, window_years=8, step_years=2, batch_size=5, seed=77)
    assert a["status"] == "ok"
    assert a["observed"]["joint_same_window_same_period_max"] == b["observed"]["joint_same_window_same_period_max"]
    p = a["observed"]["joint_same_window_same_period_max"]["period_days"]
    assert abs(p - 180.0) < 3.0
    assert a["null"]["joint"]["p_value"] == b["null"]["joint"]["p_value"]
    assert 0 < a["null"]["joint"]["p_value"] <= 1
    assert a["null"]["joint"]["null_simulations_used"] == 19
