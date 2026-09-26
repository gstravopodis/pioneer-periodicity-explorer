import numpy as np
import pandas as pd

from src.modern import lomb_scargle_periodogram, band_peak, estimate_ar1


def test_lomb_scargle_preserves_gaps_and_recovers_470d():
    rng = np.random.default_rng(1)
    dates = pd.date_range("1972-01-01", periods=5000, freq="D")
    keep = rng.random(len(dates)) > 0.2
    t = np.arange(len(dates))[keep]
    df = pd.DataFrame({"date": dates[keep], "p_11_20_mev": np.sin(2*np.pi*t/470.0)})
    pg = lomb_scargle_periodogram(df, "p_11_20_mev", min_period_days=300, max_period_days=700, nfreq=3000)
    peak = band_peak(pg, 430, 520)
    assert abs(peak["period_days"] - 470) < 20


def test_ar1_estimate_reasonable():
    rng = np.random.default_rng(2)
    phi = 0.7
    x = np.zeros(5000)
    e = rng.normal(size=len(x))
    for i in range(1, len(x)):
        x[i] = phi*x[i-1] + e[i]
    assert abs(estimate_ar1(x) - phi) < 0.05


def test_calendar_crosscorrelation_uses_real_calendar_days():
    from src.modern import calendar_crosscorrelation
    # P11 impulse on Jan 1; P10 matching impulse on Jan 11. Positive thesis-oriented lag = +10 d.
    dates11 = pd.to_datetime(["2000-01-01", "2000-01-03", "2000-01-08"])
    dates10 = pd.to_datetime(["2000-01-11", "2000-01-13", "2000-01-18"])
    p11 = pd.DataFrame({"date": dates11, "p_11_20_mev": [1.0, 2.0, 4.0]})
    p10 = pd.DataFrame({"date": dates10, "p_11_20_mev": [1.0, 2.0, 4.0]})
    out = calendar_crosscorrelation(p10, p11, "p_11_20_mev", 8, 12)
    row = out.loc[out.pearson_raw.idxmax()]
    assert int(row.lag_calendar_days) == 10
    assert np.isclose(row.pearson_raw, 1.0)
    assert int(row.pairs) == 3


def test_calendar_autocorrelation_ignores_missing_dates_as_elapsed_time():
    from src.modern import calendar_autocorrelation
    # Same shape repeats exactly 10 calendar days later despite irregular sampling inside each block.
    dates = pd.to_datetime([
        "2000-01-01", "2000-01-03", "2000-01-08",
        "2000-01-11", "2000-01-13", "2000-01-18",
    ])
    df = pd.DataFrame({"date": dates, "p_11_20_mev": [1, 2, 4, 1, 2, 4]})
    out = calendar_autocorrelation(df, "p_11_20_mev", 8, 12)
    row = out.loc[out.pearson_raw.idxmax()]
    assert int(row.lag_calendar_days) == 10
    assert np.isclose(row.pearson_raw, 1.0)


def test_sampling_window_periodogram_returns_requested_range():
    from src.modern import sampling_window_periodogram
    dates = pd.date_range("2000-01-01", periods=1000, freq="D")
    keep = np.ones(len(dates), dtype=bool)
    keep[::7] = False
    df = pd.DataFrame({"date": dates[keep], "p_11_20_mev": np.ones(keep.sum())})
    pg = sampling_window_periodogram(df, "p_11_20_mev", 20, 500)
    assert not pg.empty
    assert pg.period_days.min() >= 20
    assert pg.period_days.max() <= 500


def test_dcf_recovers_calendar_lag_with_gaps():
    from src.modern import calendar_discrete_correlation, bin_discrete_correlation, strongest_dcf_bin
    rng = np.random.default_rng(12)
    dates = pd.date_range("2001-01-01", periods=1400, freq="D")
    base = np.exp(1.2 + 0.35*np.sin(2*np.pi*np.arange(len(dates))/173.0) + rng.normal(0, 0.08, len(dates)))
    lag = 40
    # Construct P10 as a delayed copy: P10(t+lag) = P11(t).
    b = base[:-lag]
    a = base[:-lag]
    da = dates[lag:]
    db = dates[:-lag]
    keep_a = rng.random(len(da)) > 0.15
    keep_b = rng.random(len(db)) > 0.18
    p10 = pd.DataFrame({"date": da[keep_a], "p_11_20_mev": a[keep_a]})
    p11 = pd.DataFrame({"date": db[keep_b], "p_11_20_mev": b[keep_b]})
    exact = calendar_discrete_correlation(p10, "p_11_20_mev", p11, 20, 60, transform="log10")
    binned = bin_discrete_correlation(exact, bin_days=5, origin=20)
    peak = strongest_dcf_bin(binned, 25, 55)
    assert abs(peak["lag_calendar_days"] - lag) <= 5
    assert peak["dcf"] > 0.8


def test_log_lomb_scargle_recovers_irregular_period():
    from src.modern import lomb_scargle_log_periodogram, band_peak
    rng = np.random.default_rng(13)
    dates = pd.date_range("2000-01-01", periods=3000, freq="D")
    keep = rng.random(len(dates)) > 0.25
    t = np.arange(len(dates))[keep]
    flux = 10 ** (1.0 + 0.2*np.sin(2*np.pi*t/480.0))
    df = pd.DataFrame({"date": dates[keep], "p_11_20_mev": flux})
    pg = lomb_scargle_log_periodogram(df, "p_11_20_mev", 400, 560, nfreq=1600)
    peak = band_peak(pg, 430, 520)
    assert abs(peak["period_days"] - 480) < 15


def test_dcf_auto_zero_lag_near_one():
    from src.modern import calendar_discrete_correlation
    dates = pd.date_range("2000-01-01", periods=100, freq="D")
    flux = np.exp(np.linspace(0.1, 2.0, len(dates)))
    df = pd.DataFrame({"date": dates, "p_11_20_mev": flux})
    dcf = calendar_discrete_correlation(df, "p_11_20_mev", None, 0, 0, transform="log10")
    assert np.isclose(float(dcf.iloc[0].dcf), 1.0, atol=1e-12)


def test_fast_dcf_matches_existing_dcf_on_log_grid():
    from src.modern import (
        dense_log_grid, fast_calendar_dcf_from_grids, calendar_discrete_correlation
    )
    rng = np.random.default_rng(21)
    dates = pd.date_range("2002-01-01", periods=500, freq="D")
    keep = rng.random(len(dates)) > 0.18
    flux = np.exp(1 + 0.3*np.sin(2*np.pi*np.arange(len(dates))/57.0))
    df = pd.DataFrame({"date": dates[keep], "p_11_20_mev": flux[keep]})
    _, grid = dense_log_grid(df, "p_11_20_mev")
    fast = fast_calendar_dcf_from_grids(grid, None, 20, 80)
    slow = calendar_discrete_correlation(df, "p_11_20_mev", None, 20, 80, transform="log10")
    merged = fast.merge(slow[["lag_calendar_days", "pairs", "dcf"]], on="lag_calendar_days", suffixes=("_fast", "_slow"))
    assert np.array_equal(merged.pairs_fast.to_numpy(), merged.pairs_slow.to_numpy())
    assert np.allclose(merged.dcf_fast, merged.dcf_slow, atol=1e-10, equal_nan=True)


def test_daily_ar1_estimator_uses_only_adjacent_calendar_pairs():
    from src.modern import dense_log_grid, estimate_daily_ar1_from_grid
    dates = pd.to_datetime(["2000-01-01", "2000-01-02", "2000-01-10", "2000-01-11", "2000-01-12"])
    df = pd.DataFrame({"date": dates, "p_11_20_mev": np.exp([1.0, 1.1, 2.0, 2.1, 2.2])})
    _, grid = dense_log_grid(df, "p_11_20_mev")
    out = estimate_daily_ar1_from_grid(grid)
    assert out["adjacent_pairs"] == 3


def test_ar1_dcf_surrogate_test_finds_strong_nonred_periodic_acf():
    from src.modern import ar1_dcf_max_surrogate_test
    rng = np.random.default_rng(22)
    n = 1600
    dates = pd.date_range("2000-01-01", periods=n, freq="D")
    t = np.arange(n)
    flux = 10 ** (1.0 + 0.8*np.sin(2*np.pi*t/80.0) + rng.normal(0, 0.05, n))
    keep = rng.random(n) > 0.12
    df = pd.DataFrame({"date": dates[keep], "p_11_20_mev": flux[keep]})
    out = ar1_dcf_max_surrogate_test(
        df, "p_11_20_mev", {"target": (70, 90)}, simulations=79,
        min_lag_days=50, max_lag_days=110, bin_days=5, seed=3
    )
    assert out["status"] == "ok"
    assert out["bands"]["target"]["observed"]["dcf"] > 0.5
    assert out["bands"]["target"]["p_value_band_max_ar1"] <= 0.05


def test_circular_shift_ccf_test_detects_delayed_shared_signal():
    from src.modern import circular_shift_ccf_max_test
    rng = np.random.default_rng(23)
    n = 1800
    lag = 60
    dates = pd.date_range("2000-01-01", periods=n, freq="D")
    x = rng.normal(size=n)
    # smooth random process so a shifted copy has a clear DCF peak but no strict periodicity
    x = np.convolve(x, np.ones(9)/9, mode="same")
    p11_vals = 10 ** (1.0 + x[:-lag])
    p10_vals = 10 ** (1.0 + x[:-lag])
    d11 = dates[:-lag]
    d10 = dates[lag:]
    k11 = rng.random(len(d11)) > 0.1
    k10 = rng.random(len(d10)) > 0.1
    p11 = pd.DataFrame({"date": d11[k11], "p_11_20_mev": p11_vals[k11]})
    p10 = pd.DataFrame({"date": d10[k10], "p_11_20_mev": p10_vals[k10]})
    out = circular_shift_ccf_max_test(
        p10, p11, "p_11_20_mev", {"target": (50, 70)}, simulations=79,
        min_lag_days=40, max_lag_days=80, bin_days=5, min_shift_days=150, seed=4
    )
    assert out["status"] == "ok"
    assert abs(out["bands"]["target"]["observed"]["lag_calendar_days"] - lag) <= 5
    assert out["bands"]["target"]["observed"]["dcf"] > 0.5
    assert out["bands"]["target"]["p_value_band_max_circular_shift"] <= 0.05
