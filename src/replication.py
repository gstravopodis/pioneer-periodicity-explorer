from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import signal

AU_KM_PER_DAY_PER_KMS = 149_597_870.7 / 86_400.0  # 1731.456... days*km/s per AU
THESIS_AU_FACTOR = 1731.0  # value explicitly used in the dissertation code


def next_power_of_two(n: int) -> int:
    return 1 if n <= 1 else 2 ** int(np.ceil(np.log2(n)))


def thesis_fft(x, sample_days: float = 1.0) -> pd.DataFrame:
    """Reconstruct the dissertation FFT recipe.

    - finite samples only, preserving order;
    - zero-pad to next power of 2;
    - remove DC bin;
    - keep positive frequencies;
    - thesis power definition: 2*|FFT|^2.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 4:
        return pd.DataFrame(columns=["frequency_cpd", "period_days", "power", "phase_deg"])
    n = next_power_of_two(len(x))
    ft = np.fft.fft(x, n=n)
    freq = np.fft.fftfreq(n, d=sample_days)
    positive = freq > 0
    f = freq[positive]
    z = ft[positive]
    power = 2.0 * np.abs(z) ** 2
    phase = np.angle(z, deg=True)
    return pd.DataFrame({
        "frequency_cpd": f,
        "period_days": 1.0 / f,
        "power": power,
        "phase_deg": phase,
    }).sort_values("period_days").reset_index(drop=True)


def xcorr_coeff(x, y=None) -> tuple[np.ndarray, np.ndarray]:
    """Approximate MATLAB xcorr(...,'coeff') used in the dissertation.

    The thesis fed equal-length vectors after date matching.  The normalisation is
    constant across lags: sqrt(sum(x^2) sum(y^2)).  We intentionally do not
    demean because the 1998 workflow did not do so.
    """
    x = np.asarray(x, dtype=float)
    if y is None:
        y = x.copy()
    else:
        y = np.asarray(y, dtype=float)
    if len(x) != len(y):
        raise ValueError("Replication xcorr expects equal-length sequences after alignment.")
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    c = signal.correlate(x, y, mode="full", method="fft")
    denom = np.sqrt(np.dot(x, x) * np.dot(y, y))
    if denom > 0:
        c = c / denom
    lags = signal.correlation_lags(len(x), len(y), mode="full")
    return lags.astype(int), c


def autocorrelation(x) -> pd.DataFrame:
    lags, coeff = xcorr_coeff(x)
    return pd.DataFrame({"lag_days": lags, "coefficient": coeff})


def align_on_common_dates(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                          include_geometry: bool = False) -> pd.DataFrame:
    cols = ["date", channel]
    geom = ["x_au", "y_au", "z_au", "solar_wind_speed_km_s"]
    if include_geometry:
        cols += [c for c in geom if c in df10.columns and c in df11.columns]
    a = df10[cols].rename(columns={channel: "p10", **{c: f"{c}_10" for c in geom}}).copy()
    b = df11[cols].rename(columns={channel: "p11", **{c: f"{c}_11" for c in geom}}).copy()
    out = a.merge(b, on="date", how="inner")
    out = out.dropna(subset=["p10", "p11"])
    return out.sort_values("date").reset_index(drop=True)


def crosscorrelation_common_dates(df10: pd.DataFrame, df11: pd.DataFrame, channel: str) -> pd.DataFrame:
    common = align_on_common_dates(df10, df11, channel)
    lags, coeff = xcorr_coeff(common["p10"].to_numpy(), common["p11"].to_numpy())
    return pd.DataFrame({"lag_days": lags, "coefficient": coeff})


def _mean_solar_wind_speed(v10: np.ndarray, v11: np.ndarray) -> np.ndarray:
    v10 = np.nan_to_num(v10.astype(float), nan=0.0)
    v11 = np.nan_to_num(v11.astype(float), nan=0.0)
    both = (v10 != 0) & (v11 != 0)
    only10 = (v10 != 0) & (v11 == 0)
    only11 = (v10 == 0) & (v11 != 0)
    v = np.full(len(v10), 400.0, dtype=float)
    v[both] = (v10[both] + v11[both]) / 2.0
    v[only10] = v10[only10]
    v[only11] = v11[only11]
    return v


def propagation_delay_days_from_common(common: pd.DataFrame, thesis_exact: bool = True) -> np.ndarray:
    required = [
        "x_au_10", "y_au_10", "z_au_10", "solar_wind_speed_km_s_10",
        "x_au_11", "y_au_11", "z_au_11", "solar_wind_speed_km_s_11",
    ]
    missing = [c for c in required if c not in common.columns]
    if missing:
        raise ValueError(f"Missing geometry/solar-wind fields: {', '.join(missing)}")
    r10 = np.sqrt(common.x_au_10**2 + common.y_au_10**2 + common.z_au_10**2).to_numpy(float)
    r11 = np.sqrt(common.x_au_11**2 + common.y_au_11**2 + common.z_au_11**2).to_numpy(float)
    v = _mean_solar_wind_speed(
        common.solar_wind_speed_km_s_10.to_numpy(float),
        common.solar_wind_speed_km_s_11.to_numpy(float),
    )
    factor = THESIS_AU_FACTOR if thesis_exact else AU_KM_PER_DAY_PER_KMS
    return np.rint(((r10 - r11) / v) * factor).astype(int)


def propagation_delay_days(df10: pd.DataFrame, df11: pd.DataFrame, channel: str = "p_11_20_mev",
                           thesis_exact: bool = True) -> pd.DataFrame:
    common = align_on_common_dates(df10, df11, channel, include_geometry=True)
    delay = propagation_delay_days_from_common(common, thesis_exact=thesis_exact)
    r10 = np.sqrt(common.x_au_10**2 + common.y_au_10**2 + common.z_au_10**2)
    r11 = np.sqrt(common.x_au_11**2 + common.y_au_11**2 + common.z_au_11**2)
    v = _mean_solar_wind_speed(
        common.solar_wind_speed_km_s_10.to_numpy(float),
        common.solar_wind_speed_km_s_11.to_numpy(float),
    )
    return pd.DataFrame({
        "date": common.date,
        "r10_au": r10,
        "r11_au": r11,
        "solar_wind_speed_km_s": v,
        "delay_days": delay,
    })


def shift_p11_backward_like_thesis(common: pd.DataFrame, thesis_exact: bool = True) -> pd.DataFrame:
    """Reconstruct the index-shift loop printed on dissertation pp. 53–54.

    For each common-date index i, the P11 sample becomes p11[i-dt[i]], with 0
    inserted if that index is non-positive/out of range or the P10 sample is 0.
    This is intentionally index based because that is what the 1998 MATLAB code did.
    """
    delay = propagation_delay_days_from_common(common, thesis_exact=thesis_exact)
    p10 = common["p10"].to_numpy(float)
    p11 = common["p11"].to_numpy(float)
    shifted = np.zeros(len(common), dtype=float)
    src_idx = np.full(len(common), -1, dtype=int)
    for i in range(len(common) - 1, -1, -1):
        j = i - int(delay[i])
        if j < 0 or j >= len(common) or p10[i] == 0 or not np.isfinite(p11[j]):
            shifted[i] = 0.0
        else:
            shifted[i] = p11[j]
            src_idx[i] = j
    out = common[["date", "p10", "p11"]].copy()
    out["delay_days"] = delay
    out["p11_shifted"] = shifted
    out["p11_source_index"] = src_idx
    return out


def propagation_shifted_crosscorrelation(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                                         thesis_exact: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = align_on_common_dates(df10, df11, channel, include_geometry=True)
    shifted = shift_p11_backward_like_thesis(common, thesis_exact=thesis_exact)
    lags, coeff = xcorr_coeff(shifted["p10"].to_numpy(), shifted["p11_shifted"].to_numpy())
    corr = pd.DataFrame({"lag_days": lags, "coefficient": coeff})
    return corr, shifted




def coefficient_at_lag(corr: pd.DataFrame, lag: int) -> dict:
    """Return the coefficient at an exact integer sample lag, if present."""
    q = corr[corr.lag_days == int(lag)]
    if q.empty:
        return {"lag_days": int(lag), "coefficient": np.nan}
    row = q.iloc[0]
    return {"lag_days": int(lag), "coefficient": float(row.coefficient)}


def top_correlation_peaks(corr: pd.DataFrame, n: int = 12, min_abs_lag: int = 1) -> list[dict]:
    """Return strongest local maxima by coefficient for replication diagnostics.

    This intentionally ranks positive coefficients, matching the dissertation tables,
    and retains the signed lag.
    """
    q = corr[corr.lag_days.abs() >= int(min_abs_lag)].copy()
    if q.empty:
        return []
    vals = q.coefficient.to_numpy(float)
    idx, _ = signal.find_peaks(vals)
    peaks = q.iloc[idx] if len(idx) else q
    peaks = peaks.nlargest(n, "coefficient")
    return [
        {"lag_days": int(r.lag_days), "coefficient": float(r.coefficient)}
        for r in peaks.itertuples(index=False)
    ]


def alignment_diagnostics(df10: pd.DataFrame, df11: pd.DataFrame, channel: str) -> dict:
    common = align_on_common_dates(df10, df11, channel)
    if common.empty:
        return {"common_rows": 0}
    dates = pd.to_datetime(common["date"]).sort_values().reset_index(drop=True)
    delta = dates.diff().dt.days.dropna()
    return {
        "common_rows": int(len(common)),
        "start_date": str(dates.iloc[0].date()),
        "stop_date": str(dates.iloc[-1].date()),
        "calendar_span_days": int((dates.iloc[-1] - dates.iloc[0]).days),
        "one_day_steps": int((delta == 1).sum()),
        "gaps_gt_1_day": int((delta > 1).sum()),
        "max_gap_days": int(delta.max()) if len(delta) else 0,
        "mean_step_days": float(delta.mean()) if len(delta) else np.nan,
        "median_step_days": float(delta.median()) if len(delta) else np.nan,
        "note": "As in the 1998 MATLAB code, xcorr lags are sample-index lags; the thesis labels them as days.",
    }


def thesis_crosscorr_lag_from_scipy(lag: int) -> int:
    """Translate SciPy's correlation lag sign to the convention printed in the thesis.

    The recovered 1998 benchmark amplitudes match exactly after a sign reversal:
    thesis_lag = -scipy_lag.  We keep both conventions explicit rather than
    silently mutating the raw SciPy output.
    """
    return -int(lag)


def calendar_separation_for_sample_lag(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                                       sample_lag: int) -> dict:
    """Describe the *calendar* separations represented by an index/sample lag.

    The 1998 workflow aligned the two spacecraft on common measurement dates and
    then used MATLAB xcorr.  Its lag axis was labelled in days, but after gaps are
    removed an xcorr lag is formally a number of common-date samples.  For a lag
    of |k| samples this function reports the distribution of actual calendar-day
    differences between paired common dates k rows apart.
    """
    common = align_on_common_dates(df10, df11, channel)
    k = abs(int(sample_lag))
    if common.empty or k == 0 or k >= len(common):
        return {
            "sample_lag": int(sample_lag),
            "pairs": 0,
            "mean_calendar_days": np.nan,
            "median_calendar_days": np.nan,
            "min_calendar_days": np.nan,
            "max_calendar_days": np.nan,
        }
    dates = pd.to_datetime(common["date"]).reset_index(drop=True)
    delta = (dates.iloc[k:].reset_index(drop=True) - dates.iloc[:-k].reset_index(drop=True)).dt.days
    return {
        "sample_lag": int(sample_lag),
        "pairs": int(len(delta)),
        "mean_calendar_days": float(delta.mean()),
        "median_calendar_days": float(delta.median()),
        "min_calendar_days": int(delta.min()),
        "max_calendar_days": int(delta.max()),
    }


def benchmark_match(observed: dict, expected_lag: int, expected_coeff: float,
                    lag_tolerance: int = 0, coeff_tolerance: float = 0.015) -> dict:
    """Machine-readable replication check for a pre-declared historical benchmark."""
    lag = observed.get("lag_days")
    coeff = observed.get("coefficient")
    ok = (
        lag is not None and coeff is not None and
        np.isfinite(lag) and np.isfinite(coeff) and
        abs(int(lag) - int(expected_lag)) <= int(lag_tolerance) and
        abs(float(coeff) - float(expected_coeff)) <= float(coeff_tolerance)
    )
    return {
        "expected_lag_days": int(expected_lag),
        "expected_coefficient": float(expected_coeff),
        "observed_lag_days": None if lag is None or not np.isfinite(lag) else int(lag),
        "observed_coefficient": None if coeff is None or not np.isfinite(coeff) else float(coeff),
        "lag_tolerance_days": int(lag_tolerance),
        "coefficient_tolerance": float(coeff_tolerance),
        "replicated": bool(ok),
    }

def strongest_positive_lags(acf: pd.DataFrame, threshold: float = 0.1, n: int = 20) -> pd.DataFrame:
    q = acf[(acf.lag_days > 0) & (acf.coefficient >= threshold)].copy()
    if len(q) < 3:
        return q.nlargest(n, "coefficient")
    idx, _ = signal.find_peaks(q.coefficient.to_numpy())
    peaks = q.iloc[idx]
    return peaks.nlargest(n, "coefficient").sort_values("lag_days")


def strongest_in_band(corr: pd.DataFrame, lo: int, hi: int, absolute: bool = False) -> dict:
    q = corr[(corr.lag_days >= lo) & (corr.lag_days <= hi)].copy()
    if q.empty:
        return {"lag_days": np.nan, "coefficient": np.nan}
    score = q.coefficient.abs() if absolute else q.coefficient
    row = q.loc[score.idxmax()]
    return {"lag_days": int(row.lag_days), "coefficient": float(row.coefficient)}


def calendar_separation_for_single_sample_lag(df: pd.DataFrame, channel: str,
                                              sample_lag: int) -> dict:
    """Calendar-day separations represented by a sample/index lag in one spacecraft."""
    q = df[["date", channel]].dropna().copy()
    q["date"] = pd.to_datetime(q["date"])
    q = q.sort_values("date").drop_duplicates("date").reset_index(drop=True)
    k = abs(int(sample_lag))
    if q.empty or k == 0 or k >= len(q):
        return {
            "sample_lag": int(sample_lag), "pairs": 0,
            "mean_calendar_days": np.nan, "median_calendar_days": np.nan,
            "min_calendar_days": np.nan, "max_calendar_days": np.nan,
        }
    dates = q["date"]
    delta = (dates.iloc[k:].reset_index(drop=True) - dates.iloc[:-k].reset_index(drop=True)).dt.days
    return {
        "sample_lag": int(sample_lag),
        "pairs": int(len(delta)),
        "mean_calendar_days": float(delta.mean()),
        "median_calendar_days": float(delta.median()),
        "min_calendar_days": int(delta.min()),
        "max_calendar_days": int(delta.max()),
    }
