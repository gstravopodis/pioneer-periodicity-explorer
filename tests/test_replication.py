import numpy as np
import pandas as pd

from src.replication import (
    next_power_of_two,
    thesis_fft,
    autocorrelation,
    xcorr_coeff,
    propagation_delay_days_from_common,
    shift_p11_backward_like_thesis,
    strongest_in_band,
)
from src.io import decimal_year_1998_to_datetime, apply_p11_cpi_quality_mask


def test_next_power_of_two_matches_thesis_example():
    assert next_power_of_two(7042) == 8192


def test_decimal_year_example_from_thesis_maps_to_day_63():
    d = decimal_year_1998_to_datetime(pd.Series([72.1721])).iloc[0]
    assert d.year == 1972
    assert d.dayofyear == 63


def test_periodic_signal_recovered_near_470_days():
    n = 7042
    t = np.arange(n)
    x = np.sin(2*np.pi*t/470.0)
    pg = thesis_fft(x)
    target = pg[(pg.period_days > 430) & (pg.period_days < 520)]
    assert not target.empty
    peak = target.loc[target.power.idxmax(), "period_days"]
    assert abs(peak - 470) < 30


def test_acf_zero_lag_is_one():
    x = np.arange(1, 100, dtype=float)
    acf = autocorrelation(x)
    zero = acf.loc[acf.lag_days.eq(0), "coefficient"].iloc[0]
    assert np.isclose(zero, 1.0)


def test_crosscorr_detects_known_positive_shift():
    x = np.zeros(200)
    x[50] = 1
    y = np.zeros(200)
    y[70] = 1
    lags, c = xcorr_coeff(x, y)
    lag = lags[np.argmax(c)]
    assert abs(lag) == 20


def test_p11_quality_mask_starts_day_239_1980():
    df = pd.DataFrame({
        "date": pd.to_datetime(["1980-08-25", "1980-08-26"]),
        "p_24_29_mev": [1.0, 2.0],
        "p_11_20_mev": [3.0, 4.0],
    })
    out = apply_p11_cpi_quality_mask(df)
    assert out.loc[0, "p_24_29_mev"] == 1.0
    assert np.isnan(out.loc[1, "p_24_29_mev"])
    assert out.loc[1, "p_11_20_mev"] == 4.0


def test_thesis_propagation_factor_about_1731_days_per_au_per_kms():
    common = pd.DataFrame({
        "x_au_10": [2.0], "y_au_10": [0.0], "z_au_10": [0.0],
        "x_au_11": [1.0], "y_au_11": [0.0], "z_au_11": [0.0],
        "solar_wind_speed_km_s_10": [400.0],
        "solar_wind_speed_km_s_11": [400.0],
    })
    dt = propagation_delay_days_from_common(common, thesis_exact=True)
    assert dt[0] == 4  # round(1 AU / 400 km/s * 1731) = 4 d


def test_shift_p11_matches_index_rule():
    n = 20
    common = pd.DataFrame({
        "date": pd.date_range("1980-01-01", periods=n),
        "p10": np.arange(n, dtype=float) + 1,
        "p11": np.arange(n, dtype=float) + 100,
        "x_au_10": np.full(n, 2.0), "y_au_10": np.zeros(n), "z_au_10": np.zeros(n),
        "x_au_11": np.full(n, 1.0), "y_au_11": np.zeros(n), "z_au_11": np.zeros(n),
        "solar_wind_speed_km_s_10": np.full(n, 400.0),
        "solar_wind_speed_km_s_11": np.full(n, 400.0),
    })
    out = shift_p11_backward_like_thesis(common)
    # delay is 4 d, so index 10 uses original P11 index 6.
    assert out.loc[10, "delay_days"] == 4
    assert out.loc[10, "p11_shifted"] == common.loc[6, "p11"]


def test_strongest_in_declared_band():
    corr = pd.DataFrame({"lag_days": [330, 340, 350], "coefficient": [0.2, 0.8, 0.4]})
    r = strongest_in_band(corr, 315, 365)
    assert r["lag_days"] == 340
    assert np.isclose(r["coefficient"], 0.8)


def test_thesis_crosscorr_sign_translation():
    from src.replication import thesis_crosscorr_lag_from_scipy
    assert thesis_crosscorr_lag_from_scipy(-341) == 341
    assert thesis_crosscorr_lag_from_scipy(214) == -214


def test_calendar_separation_for_gapped_sample_lag():
    from src.replication import calendar_separation_for_sample_lag
    dates = pd.to_datetime(["2000-01-01", "2000-01-02", "2000-01-05", "2000-01-06"])
    a = pd.DataFrame({"date": dates, "p_11_20_mev": [1, 2, 3, 4]})
    b = pd.DataFrame({"date": dates, "p_11_20_mev": [4, 3, 2, 1]})
    d = calendar_separation_for_sample_lag(a, b, "p_11_20_mev", 2)
    assert d["pairs"] == 2
    assert d["median_calendar_days"] == 4.0


def test_benchmark_match_tolerance():
    from src.replication import benchmark_match
    r = benchmark_match({"lag_days": 341, "coefficient": 0.2387}, 341, 0.24)
    assert r["replicated"] is True


def test_single_spacecraft_sample_lag_calendar_diagnostic():
    from src.replication import calendar_separation_for_single_sample_lag
    df = pd.DataFrame({
        "date": pd.to_datetime(["2000-01-01", "2000-01-02", "2000-01-05", "2000-01-06"]),
        "p_11_20_mev": [1.0, 2.0, 3.0, 4.0],
    })
    d = calendar_separation_for_single_sample_lag(df, "p_11_20_mev", 2)
    assert d["pairs"] == 2
    assert d["median_calendar_days"] == 4.0
