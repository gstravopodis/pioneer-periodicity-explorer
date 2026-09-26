from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import signal


def prepare_irregular_series(df: pd.DataFrame, channel: str, detrend: str = "linear") -> pd.DataFrame:
    """Prepare a gap-preserving series for modern irregular-time analysis.

    Unlike the 1998 replication path, missing dates are not interpolated here.
    """
    q = df[["date", channel]].dropna().sort_values("date").drop_duplicates("date").copy()
    if q.empty:
        return pd.DataFrame(columns=["date", channel, "t_days", "value_processed"])
    t = (q.date - q.date.min()).dt.total_seconds().to_numpy(float) / 86400.0
    y = q[channel].to_numpy(float)
    if detrend == "linear" and len(y) > 2:
        p = np.polyfit(t, y, 1)
        y = y - np.polyval(p, t)
    elif detrend == "constant":
        y = y - np.mean(y)
    q["t_days"] = t
    q["value_processed"] = y
    return q


def lomb_scargle_periodogram(df: pd.DataFrame, channel: str,
                             min_period_days: float = 20,
                             max_period_days: float = 4000,
                             nfreq: int = 10000,
                             detrend: str = "linear") -> pd.DataFrame:
    q = prepare_irregular_series(df, channel, detrend=detrend)
    if len(q) < 10:
        return pd.DataFrame(columns=["period_days", "power"])
    t = q.t_days.to_numpy(float)
    y = q.value_processed.to_numpy(float)
    y = y - y.mean()
    freqs = np.linspace(1 / max_period_days, 1 / min_period_days, nfreq)
    angular = 2 * np.pi * freqs
    p = signal.lombscargle(t, y, angular, normalize=True)
    return pd.DataFrame({"period_days": 1 / freqs, "power": p}).sort_values("period_days")


def band_peak(periodogram: pd.DataFrame, lo: float, hi: float) -> dict:
    q = periodogram[(periodogram.period_days >= lo) & (periodogram.period_days <= hi)]
    if q.empty:
        return {"period_days": np.nan, "power": np.nan}
    row = q.loc[q.power.idxmax()]
    return {"period_days": float(row.period_days), "power": float(row.power)}


def estimate_ar1(y: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 3:
        return 0.0
    y = y - y.mean()
    den = np.dot(y[:-1], y[:-1])
    if den == 0:
        return 0.0
    phi = float(np.dot(y[:-1], y[1:]) / den)
    return float(np.clip(phi, -0.98, 0.98))


def ar1_surrogate(phi: float, n: int, sigma: float, rng: np.random.Generator) -> np.ndarray:
    eps = rng.normal(0.0, sigma, n)
    x = np.zeros(n, dtype=float)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + eps[i]
    return x


def red_noise_band_test(df: pd.DataFrame, channel: str, lo: float, hi: float,
                        simulations: int = 500, detrend: str = "linear",
                        seed: int = 1998) -> dict:
    """Exploratory AR(1) Monte Carlo test for a pre-declared period band.

    This is a falsification aid, not yet a publication-grade global-significance
    implementation. It preserves the actual observation times and compares the
    maximum Lomb–Scargle power inside the declared band with AR(1) surrogates.
    """
    q = prepare_irregular_series(df, channel, detrend=detrend)
    if len(q) < 30:
        return {"observed_power": np.nan, "p_value": np.nan, "phi": np.nan, "simulations": 0}
    t = q.t_days.to_numpy(float)
    y = q.value_processed.to_numpy(float)
    phi = estimate_ar1(y)
    sigma = float(np.std(y) * np.sqrt(max(1e-9, 1 - phi**2)))
    freqs = np.linspace(1 / hi, 1 / lo, 800)
    ang = 2 * np.pi * freqs
    obs = signal.lombscargle(t, y - y.mean(), ang, normalize=True).max()
    rng = np.random.default_rng(seed)
    exceed = 0
    # Generate on observation index; this is intentionally a simple null model.
    for _ in range(int(simulations)):
        s = ar1_surrogate(phi, len(y), sigma, rng)
        p = signal.lombscargle(t, s - s.mean(), ang, normalize=True).max()
        exceed += int(p >= obs)
    pval = (exceed + 1) / (simulations + 1)
    return {"observed_power": float(obs), "p_value": float(pval), "phi": phi, "simulations": int(simulations)}


def _dated_values(df: pd.DataFrame, channel: str) -> pd.DataFrame:
    """Return one finite value per calendar date for calendar-aware analysis."""
    q = df[["date", channel]].copy()
    q["date"] = pd.to_datetime(q["date"]).dt.normalize()
    q[channel] = pd.to_numeric(q[channel], errors="coerce")
    q = q.dropna(subset=["date", channel]).sort_values("date")
    # The recovered daily archive should already be unique; fail-safe averaging is
    # intentionally avoided because it would silently change a historical product.
    if q["date"].duplicated().any():
        raise ValueError("Calendar-aware analysis requires at most one value per date.")
    return q.reset_index(drop=True)


def _pair_statistics(x: np.ndarray, y: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n == 0:
        return {
            "pairs": 0,
            "cosine_raw": np.nan,
            "pearson_raw": np.nan,
            "pearson_log10": np.nan,
            "log_pairs": 0,
        }
    den = np.sqrt(np.dot(x, x) * np.dot(y, y))
    cosine = float(np.dot(x, y) / den) if den > 0 else np.nan
    if n >= 3 and np.std(x) > 0 and np.std(y) > 0:
        pearson = float(np.corrcoef(x, y)[0, 1])
    else:
        pearson = np.nan
    pos = (x > 0) & (y > 0)
    lx = np.log10(x[pos])
    ly = np.log10(y[pos])
    if len(lx) >= 3 and np.std(lx) > 0 and np.std(ly) > 0:
        pearson_log = float(np.corrcoef(lx, ly)[0, 1])
    else:
        pearson_log = np.nan
    return {
        "pairs": int(n),
        "cosine_raw": cosine,
        "pearson_raw": pearson,
        "pearson_log10": pearson_log,
        "log_pairs": int(len(lx)),
    }


def calendar_autocorrelation(df: pd.DataFrame, channel: str,
                             min_lag_days: int = 1, max_lag_days: int = 2500,
                             step_days: int = 1) -> pd.DataFrame:
    """Calendar-aware ACF using exact date separations, not row-index lags.

    For positive lag L, pairs are x(t+L) with x(t), but only where both calendar
    dates are actually observed.  Results include the raw cosine coefficient
    (closest in spirit to the non-demeaned 1998 xcorr normalisation), ordinary
    Pearson r, and Pearson r after log10-transforming strictly positive fluxes.
    """
    q = _dated_values(df, channel).rename(columns={channel: "value"})
    base = q.rename(columns={"date": "date_earlier", "value": "y"})
    rows = []
    for lag in range(int(min_lag_days), int(max_lag_days) + 1, int(step_days)):
        later = q.copy()
        later["date_earlier"] = later["date"] - pd.to_timedelta(lag, unit="D")
        later = later.rename(columns={"value": "x"})[["date_earlier", "x"]]
        m = later.merge(base, on="date_earlier", how="inner")
        stats = _pair_statistics(m["x"].to_numpy(), m["y"].to_numpy()) if len(m) else _pair_statistics([], [])
        rows.append({"lag_calendar_days": lag, **stats})
    return pd.DataFrame(rows)


def calendar_crosscorrelation(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                              min_lag_days: int = -2500, max_lag_days: int = 2500,
                              step_days: int = 1) -> pd.DataFrame:
    """Calendar-aware P10/P11 CCF in an explicit thesis-oriented convention.

    A positive lag L pairs P10 at calendar date (t + L) with P11 at date t.
    This convention mirrors the sign used in the dissertation tables after the
    empirically established translation ``thesis_lag = -scipy_lag``.
    """
    p10 = _dated_values(df10, channel).rename(columns={channel: "p10"})
    p11 = _dated_values(df11, channel).rename(columns={channel: "p11", "date": "date_p11"})
    rows = []
    for lag in range(int(min_lag_days), int(max_lag_days) + 1, int(step_days)):
        a = p10.copy()
        a["date_p11"] = a["date"] - pd.to_timedelta(lag, unit="D")
        a = a[["date_p11", "p10"]]
        m = a.merge(p11, on="date_p11", how="inner")
        stats = _pair_statistics(m["p10"].to_numpy(), m["p11"].to_numpy()) if len(m) else _pair_statistics([], [])
        rows.append({"lag_calendar_days": lag, **stats})
    return pd.DataFrame(rows)


def strongest_calendar_lag(table: pd.DataFrame, lo: int, hi: int,
                           metric: str = "pearson_log10", min_pairs: int = 100) -> dict:
    q = table[(table.lag_calendar_days >= lo) & (table.lag_calendar_days <= hi)].copy()
    q = q[q["pairs"] >= int(min_pairs)]
    q = q[np.isfinite(q[metric])]
    if q.empty:
        return {"lag_calendar_days": np.nan, metric: np.nan, "pairs": 0}
    row = q.loc[q[metric].idxmax()]
    return {
        "lag_calendar_days": int(row.lag_calendar_days),
        metric: float(row[metric]),
        "pairs": int(row.pairs),
        "cosine_raw": float(row.cosine_raw) if np.isfinite(row.cosine_raw) else np.nan,
        "pearson_raw": float(row.pearson_raw) if np.isfinite(row.pearson_raw) else np.nan,
        "pearson_log10": float(row.pearson_log10) if np.isfinite(row.pearson_log10) else np.nan,
    }


def calendar_lag_exact(table: pd.DataFrame, lag_days: int) -> dict:
    q = table[table.lag_calendar_days == int(lag_days)]
    if q.empty:
        return {"lag_calendar_days": int(lag_days), "pairs": 0}
    row = q.iloc[0]
    return {
        "lag_calendar_days": int(lag_days),
        "pairs": int(row.pairs),
        "cosine_raw": float(row.cosine_raw) if np.isfinite(row.cosine_raw) else np.nan,
        "pearson_raw": float(row.pearson_raw) if np.isfinite(row.pearson_raw) else np.nan,
        "pearson_log10": float(row.pearson_log10) if np.isfinite(row.pearson_log10) else np.nan,
        "log_pairs": int(row.log_pairs),
    }


def sampling_window_periodogram(df: pd.DataFrame, channel: str,
                                min_period_days: float = 20,
                                max_period_days: float = 2500) -> pd.DataFrame:
    """Spectrum of the observation-availability window on a complete daily grid.

    The series is 1 on dates with a finite observation and 0 on missing dates.
    A peak near one year would warn that sampling itself contains annual structure.
    """
    q = _dated_values(df, channel)
    if len(q) < 10:
        return pd.DataFrame(columns=["period_days", "power"])
    idx = pd.date_range(q.date.min(), q.date.max(), freq="D")
    obs = pd.Series(0.0, index=idx)
    obs.loc[pd.DatetimeIndex(q.date)] = 1.0
    y = obs.to_numpy(float)
    y = y - y.mean()
    freq, power = signal.periodogram(y, fs=1.0, detrend=False, scaling="spectrum")
    mask = (freq > 0)
    freq = freq[mask]
    power = power[mask]
    period = 1.0 / freq
    keep = (period >= float(min_period_days)) & (period <= float(max_period_days))
    return pd.DataFrame({"period_days": period[keep], "power": power[keep]}).sort_values("period_days")


def strongest_period_in_band(periodogram: pd.DataFrame, lo: float, hi: float) -> dict:
    q = periodogram[(periodogram.period_days >= lo) & (periodogram.period_days <= hi)]
    if q.empty:
        return {"period_days": np.nan, "power": np.nan}
    row = q.loc[q.power.idxmax()]
    return {"period_days": float(row.period_days), "power": float(row.power)}


def _transformed_daily_grid(df: pd.DataFrame, channel: str, transform: str = "log10",
                            start: pd.Timestamp | None = None,
                            stop: pd.Timestamp | None = None) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """Put one transformed observation per date on a dense daily grid with NaNs for gaps.

    This does not interpolate.  It is only a convenient representation for exact
    calendar-lag calculations.  ``transform='log10'`` drops non-positive fluxes.
    """
    q = _dated_values(df, channel)
    if q.empty:
        return pd.DatetimeIndex([]), np.array([], dtype=float)
    values = q[channel].to_numpy(float)
    if transform == "log10":
        values = np.where(values > 0, np.log10(values), np.nan)
    elif transform == "raw":
        pass
    else:
        raise ValueError("transform must be 'raw' or 'log10'")
    q = q.assign(_value=values).dropna(subset=["_value"])
    if q.empty:
        return pd.DatetimeIndex([]), np.array([], dtype=float)
    lo = pd.Timestamp(start).normalize() if start is not None else q.date.min().normalize()
    hi = pd.Timestamp(stop).normalize() if stop is not None else q.date.max().normalize()
    idx = pd.date_range(lo, hi, freq="D")
    grid = np.full(len(idx), np.nan, dtype=float)
    offsets = (pd.DatetimeIndex(q.date) - lo).days.to_numpy(int)
    ok = (offsets >= 0) & (offsets < len(grid))
    grid[offsets[ok]] = q.loc[ok, "_value"].to_numpy(float)
    return idx, grid


def _standardize_observed(x: np.ndarray) -> np.ndarray:
    out = np.asarray(x, dtype=float).copy()
    finite = np.isfinite(out)
    if finite.sum() < 2:
        return out
    mu = float(np.mean(out[finite]))
    sigma = float(np.std(out[finite], ddof=0))
    if sigma <= 0:
        out[finite] = 0.0
    else:
        out[finite] = (out[finite] - mu) / sigma
    return out


def calendar_discrete_correlation(df_a: pd.DataFrame, channel: str,
                                  df_b: pd.DataFrame | None = None,
                                  min_lag_days: int = 1,
                                  max_lag_days: int = 600,
                                  transform: str = "log10") -> pd.DataFrame:
    """Discrete correlation function on the actual daily calendar grid.

    The implementation follows the spirit of Edelson & Krolik: each observed
    value is standardized using the global mean/std of its own series, pairwise
    products are associated with their true time separation, and no interpolation
    is used.  Here the archive timestamps are integer days, so one-day lag bins
    can be computed exactly and later aggregated to wider DCF bins.

    For cross-correlation a positive lag L means A(t+L) is paired with B(t),
    matching the dissertation-oriented sign convention used elsewhere.
    """
    if df_b is None:
        df_b = df_a
    qa = _dated_values(df_a, channel)
    qb = _dated_values(df_b, channel)
    if qa.empty or qb.empty:
        return pd.DataFrame(columns=["lag_calendar_days", "pairs", "dcf", "stderr", "sum_udcf", "sum_sq_udcf"])
    start = min(qa.date.min(), qb.date.min())
    stop = max(qa.date.max(), qb.date.max())
    _, a = _transformed_daily_grid(df_a, channel, transform, start=start, stop=stop)
    _, b = _transformed_daily_grid(df_b, channel, transform, start=start, stop=stop)
    a = _standardize_observed(a)
    b = _standardize_observed(b)
    ngrid = len(a)
    rows: list[dict] = []
    for lag in range(int(min_lag_days), int(max_lag_days) + 1):
        k = int(lag)
        if abs(k) >= ngrid:
            rows.append({"lag_calendar_days": k, "pairs": 0, "dcf": np.nan, "stderr": np.nan,
                         "sum_udcf": 0.0, "sum_sq_udcf": 0.0})
            continue
        if k > 0:
            x, y = a[k:], b[:-k]
        elif k < 0:
            kk = -k
            x, y = a[:-kk], b[kk:]
        else:
            x, y = a, b
        mask = np.isfinite(x) & np.isfinite(y)
        prod = x[mask] * y[mask]
        n = int(len(prod))
        if n == 0:
            mean = se = np.nan
            s = ss = 0.0
        else:
            s = float(np.sum(prod))
            ss = float(np.dot(prod, prod))
            mean = s / n
            if n > 1:
                var = max(0.0, (ss - (s * s) / n) / (n - 1))
                se = float(np.sqrt(var / n))
            else:
                se = np.nan
        rows.append({"lag_calendar_days": k, "pairs": n, "dcf": mean, "stderr": se,
                     "sum_udcf": s, "sum_sq_udcf": ss})
    return pd.DataFrame(rows)


def bin_discrete_correlation(dcf: pd.DataFrame, bin_days: int = 10,
                             origin: int | None = None) -> pd.DataFrame:
    """Aggregate exact-day DCF values into wider time-lag bins by pair count."""
    if dcf.empty:
        return pd.DataFrame(columns=["lag_bin_start", "lag_bin_end", "lag_bin_center", "pairs", "dcf", "stderr"])
    width = int(bin_days)
    if width < 1:
        raise ValueError("bin_days must be >= 1")
    q = dcf.copy()
    if origin is None:
        origin = int(np.floor(q.lag_calendar_days.min() / width) * width)
    q["_bin"] = np.floor((q.lag_calendar_days - origin) / width).astype(int)
    rows = []
    for b, g in q.groupby("_bin", sort=True):
        n = int(g.pairs.sum())
        lo = int(origin + b * width)
        hi = int(lo + width - 1)
        if n <= 0:
            rows.append({"lag_bin_start": lo, "lag_bin_end": hi, "lag_bin_center": (lo + hi) / 2,
                         "pairs": 0, "dcf": np.nan, "stderr": np.nan})
            continue
        s = float(g.sum_udcf.sum())
        ss = float(g.sum_sq_udcf.sum())
        mean = s / n
        if n > 1:
            var = max(0.0, (ss - (s * s) / n) / (n - 1))
            se = float(np.sqrt(var / n))
        else:
            se = np.nan
        rows.append({"lag_bin_start": lo, "lag_bin_end": hi, "lag_bin_center": (lo + hi) / 2,
                     "pairs": n, "dcf": mean, "stderr": se})
    return pd.DataFrame(rows)


def strongest_dcf_bin(binned: pd.DataFrame, lo: float, hi: float,
                      min_pairs: int = 100) -> dict:
    if binned.empty:
        return {"lag_calendar_days": np.nan, "dcf": np.nan, "stderr": np.nan, "pairs": 0}
    q = binned[(binned.lag_bin_center >= lo) & (binned.lag_bin_center <= hi) & (binned.pairs >= int(min_pairs))]
    q = q[np.isfinite(q.dcf)]
    if q.empty:
        return {"lag_calendar_days": np.nan, "dcf": np.nan, "stderr": np.nan, "pairs": 0}
    row = q.loc[q.dcf.idxmax()]
    return {
        "lag_calendar_days": float(row.lag_bin_center),
        "bin_start": int(row.lag_bin_start),
        "bin_end": int(row.lag_bin_end),
        "dcf": float(row.dcf),
        "stderr": float(row.stderr) if np.isfinite(row.stderr) else np.nan,
        "pairs": int(row.pairs),
    }


def lomb_scargle_log_periodogram(df: pd.DataFrame, channel: str,
                                 min_period_days: float = 20,
                                 max_period_days: float = 2500,
                                 nfreq: int = 6000,
                                 detrend: str = "linear") -> pd.DataFrame:
    """Generalized Lomb–Scargle-style spectrum of log10 positive flux values.

    Observation times remain irregular.  A floating mean is used when supported
    by the installed SciPy; older SciPy falls back to pre-centering.
    """
    q = _dated_values(df, channel)
    q = q[q[channel] > 0].copy()
    if len(q) < 10:
        return pd.DataFrame(columns=["period_days", "power"])
    t = (q.date - q.date.min()).dt.total_seconds().to_numpy(float) / 86400.0
    y = np.log10(q[channel].to_numpy(float))
    if detrend == "linear" and len(y) > 2:
        p = np.polyfit(t, y, 1)
        y = y - np.polyval(p, t)
    elif detrend == "constant":
        y = y - np.mean(y)
    freqs = np.linspace(1 / float(max_period_days), 1 / float(min_period_days), int(nfreq))
    ang = 2 * np.pi * freqs
    try:
        power = signal.lombscargle(t, y, ang, normalize=True, floating_mean=True)
    except TypeError:
        power = signal.lombscargle(t, y - np.mean(y), ang, normalize=True, precenter=False)
    return pd.DataFrame({"period_days": 1 / freqs, "power": power}).sort_values("period_days").reset_index(drop=True)


def dense_log_grid(df: pd.DataFrame, channel: str, *, start=None, stop=None,
                   detrend: str = "none") -> tuple[pd.DatetimeIndex, np.ndarray]:
    """Return log10-positive flux on a complete daily grid with NaNs for gaps.

    ``detrend='linear'`` removes a least-squares linear trend from the observed
    log10 values as a robustness variant.  No missing value is interpolated.
    """
    q = _dated_values(df, channel)
    q = q[q[channel] > 0].copy()
    if q.empty:
        return pd.DatetimeIndex([]), np.array([], dtype=float)
    lo = pd.Timestamp(start).normalize() if start is not None else q.date.min().normalize()
    hi = pd.Timestamp(stop).normalize() if stop is not None else q.date.max().normalize()
    idx = pd.date_range(lo, hi, freq="D")
    vals = np.log10(q[channel].to_numpy(float))
    t = (pd.DatetimeIndex(q.date) - lo).days.to_numpy(float)
    if detrend == "linear" and len(vals) > 2:
        p = np.polyfit(t, vals, 1)
        vals = vals - np.polyval(p, t)
    elif detrend in ("none", None):
        pass
    else:
        raise ValueError("detrend must be 'none' or 'linear'")
    grid = np.full(len(idx), np.nan, dtype=float)
    offsets = t.astype(int)
    ok = (offsets >= 0) & (offsets < len(grid))
    grid[offsets[ok]] = vals[ok]
    return idx, grid


def estimate_daily_ar1_from_grid(grid: np.ndarray) -> dict:
    """Estimate a one-calendar-day AR(1) coefficient from genuinely adjacent observations.

    Consecutive *rows* separated by gaps are deliberately not treated as one-day
    transitions.  The estimate is a zero-intercept regression after centering on
    the observed mean.  Returned diagnostics expose how many adjacent-day pairs
    informed the estimate.
    """
    x = np.asarray(grid, dtype=float)
    finite = np.isfinite(x)
    n_obs = int(finite.sum())
    if n_obs < 3 or len(x) < 2:
        return {"phi_1day": 0.0, "adjacent_pairs": 0, "observations": n_obs}
    mu = float(np.nanmean(x))
    prev = x[:-1] - mu
    nxt = x[1:] - mu
    pair = np.isfinite(prev) & np.isfinite(nxt)
    n_pair = int(pair.sum())
    if n_pair < 3:
        return {"phi_1day": 0.0, "adjacent_pairs": n_pair, "observations": n_obs}
    den = float(np.dot(prev[pair], prev[pair]))
    phi = float(np.dot(prev[pair], nxt[pair]) / den) if den > 0 else 0.0
    phi = float(np.clip(phi, -0.995, 0.995))
    return {"phi_1day": phi, "adjacent_pairs": n_pair, "observations": n_obs}


def simulate_daily_ar1(phi: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Generate a stationary unit-variance daily AR(1) process."""
    n = int(n)
    if n <= 0:
        return np.array([], dtype=float)
    phi = float(np.clip(phi, -0.995, 0.995))
    sigma_eps = float(np.sqrt(max(1e-12, 1.0 - phi * phi)))
    out = np.empty(n, dtype=float)
    out[0] = rng.normal()
    eps = rng.normal(0.0, sigma_eps, n - 1)
    for i in range(1, n):
        out[i] = phi * out[i - 1] + eps[i - 1]
    return out


def _standardized_zero_and_mask(grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(grid, dtype=float)
    mask = np.isfinite(x)
    z = np.zeros(len(x), dtype=float)
    if mask.sum() >= 2:
        vals = x[mask]
        sd = float(np.std(vals, ddof=0))
        if sd > 0:
            z[mask] = (vals - float(np.mean(vals))) / sd
        else:
            z[mask] = 0.0
    return z, mask.astype(float)


def fast_calendar_dcf_from_grids(grid_a: np.ndarray, grid_b: np.ndarray | None = None,
                                 min_lag_days: int = 1, max_lag_days: int = 600) -> pd.DataFrame:
    """Fast exact-calendar DCF using FFT correlation on dense NaN-preserving grids.

    Positive lag follows the project convention: A(t+L) is paired with B(t).
    Results are algebraically equivalent to the pair-product mean used by
    ``calendar_discrete_correlation`` for the same transformed input grids.
    """
    a = np.asarray(grid_a, dtype=float)
    b = a if grid_b is None else np.asarray(grid_b, dtype=float)
    if len(a) != len(b):
        raise ValueError("grid_a and grid_b must share the same calendar grid")
    if len(a) == 0:
        return pd.DataFrame(columns=["lag_calendar_days", "pairs", "dcf"])
    za, ma = _standardized_zero_and_mask(a)
    zb, mb = _standardized_zero_and_mask(b)
    sums = signal.correlate(za, zb, mode="full", method="fft")
    counts = signal.correlate(ma, mb, mode="full", method="fft")
    lags = signal.correlation_lags(len(za), len(zb), mode="full")
    lo, hi = int(min_lag_days), int(max_lag_days)
    keep = (lags >= lo) & (lags <= hi)
    ls = lags[keep].astype(int)
    ss = sums[keep]
    nn = np.rint(np.maximum(counts[keep], 0.0)).astype(int)
    dcf = np.full(len(ls), np.nan, dtype=float)
    valid = nn > 0
    dcf[valid] = ss[valid] / nn[valid]
    return pd.DataFrame({"lag_calendar_days": ls, "pairs": nn, "dcf": dcf})


def fast_bin_dcf(dcf: pd.DataFrame, bin_days: int = 10, origin: int | None = None) -> pd.DataFrame:
    """Pair-count-weighted bins for the compact fast DCF representation."""
    if dcf.empty:
        return pd.DataFrame(columns=["lag_bin_start", "lag_bin_end", "lag_bin_center", "pairs", "dcf"])
    width = int(bin_days)
    if width < 1:
        raise ValueError("bin_days must be >= 1")
    q = dcf.copy()
    if origin is None:
        origin = int(np.floor(q.lag_calendar_days.min() / width) * width)
    q["_bin"] = np.floor((q.lag_calendar_days - origin) / width).astype(int)
    rows = []
    for b, g in q.groupby("_bin", sort=True):
        lo = int(origin + b * width)
        hi = int(lo + width - 1)
        n = int(g.pairs.sum())
        good = (g.pairs > 0) & np.isfinite(g.dcf)
        if n > 0 and good.any():
            s = float(np.sum(g.loc[good, "dcf"].to_numpy(float) * g.loc[good, "pairs"].to_numpy(float)))
            mean = s / int(g.loc[good, "pairs"].sum())
        else:
            mean = np.nan
        rows.append({"lag_bin_start": lo, "lag_bin_end": hi, "lag_bin_center": (lo + hi) / 2,
                     "pairs": n, "dcf": mean})
    return pd.DataFrame(rows)


def dcf_band_max(binned: pd.DataFrame, lo: float, hi: float, min_pairs: int = 100) -> dict:
    """Maximum positive DCF bin in a predeclared calendar-lag band."""
    q = binned[(binned.lag_bin_center >= float(lo)) & (binned.lag_bin_center <= float(hi))]
    q = q[(q.pairs >= int(min_pairs)) & np.isfinite(q.dcf)]
    if q.empty:
        return {"lag_calendar_days": np.nan, "dcf": np.nan, "pairs": 0}
    row = q.loc[q.dcf.idxmax()]
    return {"lag_calendar_days": float(row.lag_bin_center), "dcf": float(row.dcf), "pairs": int(row.pairs)}


def ar1_dcf_max_surrogate_test(df: pd.DataFrame, channel: str, bands: dict[str, tuple[int, int]],
                               *, simulations: int = 499, min_lag_days: int = 250,
                               max_lag_days: int = 600, bin_days: int = 10,
                               detrend: str = "none", seed: int = 1998) -> dict:
    """Band-wise max-DCF significance against a daily AR(1) red-noise null.

    The null is simulated on the complete daily calendar and then sampled using
    the *actual* observation mask.  Therefore gaps are preserved and are never
    compressed into sample-index time.  The test is one-sided and internally
    corrects for searching all DCF bins inside each predeclared band by comparing
    the observed band maximum with the surrogate band maxima.
    """
    idx, grid = dense_log_grid(df, channel, detrend=detrend)
    if len(grid) < 10 or np.isfinite(grid).sum() < 30:
        return {"status": "insufficient_data", "simulations": 0}
    ar = estimate_daily_ar1_from_grid(grid)
    observed = fast_bin_dcf(
        fast_calendar_dcf_from_grids(grid, None, min_lag_days, max_lag_days),
        bin_days=bin_days, origin=min_lag_days,
    )
    obs = {name: dcf_band_max(observed, lo, hi) for name, (lo, hi) in bands.items()}
    null = {name: [] for name in bands}
    mask = np.isfinite(grid)
    rng = np.random.default_rng(seed)
    for _ in range(int(simulations)):
        s = simulate_daily_ar1(ar["phi_1day"], len(grid), rng)
        sg = np.where(mask, s, np.nan)
        b = fast_bin_dcf(
            fast_calendar_dcf_from_grids(sg, None, min_lag_days, max_lag_days),
            bin_days=bin_days, origin=min_lag_days,
        )
        for name, (lo, hi) in bands.items():
            m = dcf_band_max(b, lo, hi)["dcf"]
            if np.isfinite(m):
                null[name].append(float(m))
    result_bands = {}
    for name in bands:
        arr = np.asarray(null[name], dtype=float)
        o = float(obs[name]["dcf"]) if np.isfinite(obs[name]["dcf"]) else np.nan
        if len(arr) and np.isfinite(o):
            exceed = int(np.sum(arr >= o))
            p = (exceed + 1) / (len(arr) + 1)
            q95 = float(np.quantile(arr, 0.95))
            q99 = float(np.quantile(arr, 0.99))
        else:
            p = q95 = q99 = np.nan
        result_bands[name] = {
            "observed": obs[name],
            "p_value_band_max_ar1": float(p) if np.isfinite(p) else np.nan,
            "null_q95": q95,
            "null_q99": q99,
            "null_simulations_used": int(len(arr)),
        }
    return {
        "status": "ok",
        "detrend": detrend,
        "transform": "log10 positive flux",
        "calendar_days": {"start": str(idx.min().date()), "stop": str(idx.max().date()), "grid_days": int(len(idx))},
        "ar1": ar,
        "simulations": int(simulations),
        "bands": result_bands,
    }


def circular_shift_ccf_max_test(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                                bands: dict[str, tuple[int, int]], *, simulations: int = 499,
                                min_lag_days: int = 250, max_lag_days: int = 600,
                                bin_days: int = 10, min_shift_days: int = 700,
                                detrend: str = "none", seed: int = 1998) -> dict:
    """Cross-spacecraft max-DCF test using random circular calendar shifts.

    Both the P11 values and its missing-data mask are shifted together inside the
    common calendar span.  This preserves P11's internal temporal structure and
    sampling pattern while breaking its absolute alignment to P10.  It is a
    deliberately stricter complement to an independent-red-noise null for CCFs.
    """
    q10 = _dated_values(df10, channel)
    q11 = _dated_values(df11, channel)
    if q10.empty or q11.empty:
        return {"status": "insufficient_data", "simulations": 0}
    start = max(q10.date.min(), q11.date.min())
    stop = min(q10.date.max(), q11.date.max())
    if stop <= start:
        return {"status": "no_overlap", "simulations": 0}
    idx, g10 = dense_log_grid(df10, channel, start=start, stop=stop, detrend=detrend)
    _, g11 = dense_log_grid(df11, channel, start=start, stop=stop, detrend=detrend)
    n = len(idx)
    if n < 2 * int(min_shift_days) + 1:
        return {"status": "span_too_short_for_shift_null", "simulations": 0}
    observed = fast_bin_dcf(
        fast_calendar_dcf_from_grids(g10, g11, min_lag_days, max_lag_days),
        bin_days=bin_days, origin=min_lag_days,
    )
    obs = {name: dcf_band_max(observed, lo, hi) for name, (lo, hi) in bands.items()}
    allowed = np.arange(int(min_shift_days), n - int(min_shift_days) + 1, dtype=int)
    if len(allowed) == 0:
        return {"status": "no_allowed_shifts", "simulations": 0}
    rng = np.random.default_rng(seed)
    shifts = rng.choice(allowed, size=int(simulations), replace=True)
    null = {name: [] for name in bands}
    for s in shifts:
        rolled = np.roll(g11, int(s))
        b = fast_bin_dcf(
            fast_calendar_dcf_from_grids(g10, rolled, min_lag_days, max_lag_days),
            bin_days=bin_days, origin=min_lag_days,
        )
        for name, (lo, hi) in bands.items():
            m = dcf_band_max(b, lo, hi)["dcf"]
            if np.isfinite(m):
                null[name].append(float(m))
    result_bands = {}
    for name in bands:
        arr = np.asarray(null[name], dtype=float)
        o = float(obs[name]["dcf"]) if np.isfinite(obs[name]["dcf"]) else np.nan
        if len(arr) and np.isfinite(o):
            exceed = int(np.sum(arr >= o))
            p = (exceed + 1) / (len(arr) + 1)
            q95 = float(np.quantile(arr, 0.95))
            q99 = float(np.quantile(arr, 0.99))
        else:
            p = q95 = q99 = np.nan
        result_bands[name] = {
            "observed": obs[name],
            "p_value_band_max_circular_shift": float(p) if np.isfinite(p) else np.nan,
            "null_q95": q95,
            "null_q99": q99,
            "null_simulations_used": int(len(arr)),
        }
    return {
        "status": "ok",
        "detrend": detrend,
        "transform": "log10 positive flux",
        "common_calendar": {"start": str(idx.min().date()), "stop": str(idx.max().date()), "grid_days": int(n)},
        "min_abs_shift_days": int(min_shift_days),
        "simulations": int(simulations),
        "bands": result_bands,
    }


def _classic_lomb_basis(t_days: np.ndarray, periods_days: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Precompute the classical Lomb-Scargle sine/cosine basis.

    The implementation uses the standard phase shift tau so the sine and cosine
    terms are orthogonal at each tested frequency. It is intentionally separate
    from SciPy's floating-mean implementation; B9 uses the same statistic for the
    observed series and every surrogate, so the Monte-Carlo p-value is internally
    calibrated to this exact statistic.
    """
    t = np.asarray(t_days, dtype=float)
    p = np.asarray(periods_days, dtype=float)
    if t.ndim != 1 or p.ndim != 1 or len(t) < 3 or len(p) < 2:
        raise ValueError("t_days and periods_days must be one-dimensional with sufficient samples")
    w = 2.0 * np.pi / p
    wt2 = 2.0 * np.outer(t, w)
    tau = np.arctan2(np.sin(wt2).sum(axis=0), np.cos(wt2).sum(axis=0)) / (2.0 * w)
    phase = (t[:, None] - tau[None, :]) * w[None, :]
    c = np.cos(phase)
    s = np.sin(phase)
    c2 = np.sum(c * c, axis=0)
    s2 = np.sum(s * s, axis=0)
    return c, s, c2, s2


def _classic_lomb_power_batch(y: np.ndarray, c: np.ndarray, s: np.ndarray,
                              c2: np.ndarray, s2: np.ndarray) -> np.ndarray:
    """Return normalized classical Lomb-Scargle power for one or many rows."""
    a = np.asarray(y, dtype=float)
    one = a.ndim == 1
    if one:
        a = a[None, :]
    if a.ndim != 2 or a.shape[1] != c.shape[0]:
        raise ValueError("y shape does not match Lomb basis")
    denom = np.sum(a * a, axis=1)
    yc = a @ c
    ys = a @ s
    with np.errstate(divide="ignore", invalid="ignore"):
        power = ((yc * yc) / c2[None, :] + (ys * ys) / s2[None, :]) / denom[:, None]
    power[~np.isfinite(power)] = np.nan
    return power[0] if one else power


def _center_or_linear_detrend_batch(y: np.ndarray, t_days: np.ndarray, detrend: str) -> np.ndarray:
    """Apply the same mean/linear treatment to observed and surrogate samples."""
    a = np.asarray(y, dtype=float)
    one = a.ndim == 1
    if one:
        a = a[None, :]
    out = a - np.mean(a, axis=1, keepdims=True)
    if detrend == "linear":
        tt = np.asarray(t_days, dtype=float)
        tc = tt - np.mean(tt)
        den = float(np.dot(tc, tc))
        if den > 0:
            slopes = (out @ tc) / den
            out = out - slopes[:, None] * tc[None, :]
    elif detrend in ("none", None, "constant"):
        pass
    else:
        raise ValueError("detrend must be 'none', 'constant' or 'linear'")
    return out[0] if one else out


def ar1_lomb_band_max_test(df: pd.DataFrame, channel: str, lo: float, hi: float, *,
                           simulations: int = 1999, detrend: str = "none",
                           nfreq: int = 96, batch_size: int = 32,
                           seed: int = 1998) -> dict:
    """Band-max spectral significance against a mask-preserving daily AR(1) null.

    This is the B9 candidate-level spectral test. A daily AR(1) process is
    simulated on the complete calendar, then sampled through the exact positive-
    flux observation mask. The test statistic is the maximum *classical* Lomb-
    Scargle power anywhere inside the predeclared period band. The same basis and
    preprocessing are used for the observed series and every surrogate.

    The routine deliberately reports a Monte-Carlo band-max p-value, not a
    single-frequency false-alarm probability. Candidate selection occurred in
    earlier phases, so B9 remains post-selection evidence rather than an
    independent replication.
    """
    idx, grid = dense_log_grid(df, channel, detrend="none")
    mask = np.isfinite(grid)
    if len(grid) < 30 or int(mask.sum()) < 30:
        return {"status": "insufficient_data", "simulations": 0}
    t = np.arange(len(grid), dtype=float)[mask]
    y_obs = grid[mask]
    y_obs = _center_or_linear_detrend_batch(y_obs, t, detrend)
    periods = np.linspace(float(lo), float(hi), int(nfreq), dtype=float)
    c, s, c2, s2 = _classic_lomb_basis(t, periods)
    p_obs = _classic_lomb_power_batch(y_obs, c, s, c2, s2)
    if not np.isfinite(p_obs).any():
        return {"status": "invalid_observed_periodogram", "simulations": 0}
    imax = int(np.nanargmax(p_obs))
    obs_max = float(p_obs[imax])
    obs_period = float(periods[imax])

    ar = estimate_daily_ar1_from_grid(grid)
    phi = float(ar["phi_1day"])
    rng = np.random.default_rng(seed)
    maxima: list[float] = []
    n_total = int(simulations)
    bs = max(1, int(batch_size))
    for start in range(0, n_total, bs):
        b = min(bs, n_total - start)
        x = np.empty((b, len(grid)), dtype=float)
        x[:, 0] = rng.normal(size=b)
        eps_sd = float(np.sqrt(max(1e-12, 1.0 - phi * phi)))
        for k in range(1, len(grid)):
            x[:, k] = phi * x[:, k - 1] + rng.normal(0.0, eps_sd, size=b)
        yy = x[:, mask]
        yy = _center_or_linear_detrend_batch(yy, t, detrend)
        pp = _classic_lomb_power_batch(yy, c, s, c2, s2)
        maxima.extend(np.nanmax(pp, axis=1).astype(float).tolist())

    arr = np.asarray(maxima, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr):
        exceed = int(np.sum(arr >= obs_max))
        pval = float((exceed + 1) / (len(arr) + 1))
        q95 = float(np.quantile(arr, 0.95))
        q99 = float(np.quantile(arr, 0.99))
    else:
        exceed = 0
        pval = q95 = q99 = np.nan
    return {
        "status": "ok",
        "transform": "log10 positive flux",
        "detrend": detrend,
        "band_calendar_days": [float(lo), float(hi)],
        "frequency_grid_points": int(nfreq),
        "calendar": {"start": str(idx.min().date()), "stop": str(idx.max().date()), "grid_days": int(len(idx))},
        "observations": int(mask.sum()),
        "ar1": ar,
        "observed": {"period_days": obs_period, "power": obs_max},
        "p_value_band_max_ar1_lomb": pval,
        "null_q95": q95,
        "null_q99": q99,
        "exceedances": int(exceed),
        "null_simulations_used": int(len(arr)),
        "simulations": int(simulations),
        "statistic_note": "maximum classical Lomb-Scargle power inside the locked band; same statistic for observed and surrogates",
    }


def _calendar_windows(start: pd.Timestamp, stop: pd.Timestamp, *, window_years: int = 8,
                      step_years: int = 2) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Build inclusive calendar windows using year offsets, preserving leap years."""
    lo = pd.Timestamp(start).normalize()
    hi = pd.Timestamp(stop).normalize()
    out: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    cur = lo
    while True:
        end = cur + pd.DateOffset(years=int(window_years)) - pd.Timedelta(days=1)
        if end > hi:
            break
        out.append((cur, end))
        cur = cur + pd.DateOffset(years=int(step_years))
    return out


def ar1_lomb_intermit_joint_test(df10: pd.DataFrame, df11: pd.DataFrame, channel: str,
                                 lo: float, hi: float, *, simulations: int = 1999,
                                 detrend: str = "none", nfreq: int = 96,
                                 window_years: int = 8, step_years: int = 2,
                                 batch_size: int = 16, seed: int = 1998,
                                 min_observations: int = 180) -> dict:
    """Test for an intermittent, contemporaneous band-limited signal in P10/P11.

    The statistic is deliberately global over the predeclared time-frequency search:
    for each fixed calendar window, classical Lomb-Scargle power is evaluated on a
    shared period grid.  Three maxima are retained across *all* windows and all
    frequencies: P10 power, P11 power, and the joint geometric-mean power at the
    same period in the same window.  Monte-Carlo AR(1) surrogates are simulated on
    the complete common daily calendar and sampled through each spacecraft's exact
    positive-flux mask.  Thus the p-values include the look-elsewhere cost of the
    window/frequency scan itself.

    No missing values are interpolated.  The two spacecraft null processes are
    simulated independently, which is appropriate for a null of no shared
    heliophysical oscillation beyond each series' own red-noise persistence and
    observation mask.
    """
    q10 = _dated_values(df10, channel)
    q11 = _dated_values(df11, channel)
    q10 = q10[q10[channel] > 0]
    q11 = q11[q11[channel] > 0]
    if q10.empty or q11.empty:
        return {"status": "insufficient_data", "simulations": 0}

    common_start = max(pd.Timestamp(q10.date.min()), pd.Timestamp(q11.date.min())).normalize()
    common_stop = min(pd.Timestamp(q10.date.max()), pd.Timestamp(q11.date.max())).normalize()
    if common_stop <= common_start:
        return {"status": "no_common_calendar", "simulations": 0}

    idx10, g10 = dense_log_grid(df10, channel, start=common_start, stop=common_stop, detrend="none")
    idx11, g11 = dense_log_grid(df11, channel, start=common_start, stop=common_stop, detrend="none")
    if len(idx10) != len(idx11) or len(idx10) < 30:
        return {"status": "invalid_common_grid", "simulations": 0}

    windows = _calendar_windows(common_start, common_stop, window_years=window_years, step_years=step_years)
    periods = np.linspace(float(lo), float(hi), int(nfreq), dtype=float)
    prepared: list[dict] = []
    observed_windows: list[dict] = []

    best10 = {"power": -np.inf}
    best11 = {"power": -np.inf}
    bestj = {"power": -np.inf}

    for ws, we in windows:
        i0 = int((ws - common_start).days)
        i1 = int((we - common_start).days) + 1
        sub10 = g10[i0:i1]
        sub11 = g11[i0:i1]
        m10 = np.isfinite(sub10)
        m11 = np.isfinite(sub11)
        if int(m10.sum()) < int(min_observations) or int(m11.sum()) < int(min_observations):
            continue
        t10 = np.arange(i1 - i0, dtype=float)[m10]
        t11 = np.arange(i1 - i0, dtype=float)[m11]
        y10 = _center_or_linear_detrend_batch(sub10[m10], t10, detrend)
        y11 = _center_or_linear_detrend_batch(sub11[m11], t11, detrend)
        c10, s10, c210, s210 = _classic_lomb_basis(t10, periods)
        c11, s11, c211, s211 = _classic_lomb_basis(t11, periods)
        p10 = _classic_lomb_power_batch(y10, c10, s10, c210, s210)
        p11 = _classic_lomb_power_batch(y11, c11, s11, c211, s211)
        joint = np.sqrt(np.maximum(p10, 0.0) * np.maximum(p11, 0.0))
        k10 = int(np.nanargmax(p10)); k11 = int(np.nanargmax(p11)); kj = int(np.nanargmax(joint))
        winrec = {
            "window_start": str(ws.date()), "window_stop": str(we.date()),
            "p10_observations": int(m10.sum()), "p11_observations": int(m11.sum()),
            "p10_peak": {"period_days": float(periods[k10]), "power": float(p10[k10])},
            "p11_peak": {"period_days": float(periods[k11]), "power": float(p11[k11])},
            "joint_peak": {"period_days": float(periods[kj]), "geometric_mean_power": float(joint[kj]),
                           "p10_power": float(p10[kj]), "p11_power": float(p11[kj])},
        }
        observed_windows.append(winrec)
        if float(p10[k10]) > best10["power"]:
            best10 = {"power": float(p10[k10]), "period_days": float(periods[k10]),
                      "window_start": str(ws.date()), "window_stop": str(we.date())}
        if float(p11[k11]) > best11["power"]:
            best11 = {"power": float(p11[k11]), "period_days": float(periods[k11]),
                      "window_start": str(ws.date()), "window_stop": str(we.date())}
        if float(joint[kj]) > bestj["power"]:
            bestj = {"power": float(joint[kj]), "period_days": float(periods[kj]),
                     "p10_power": float(p10[kj]), "p11_power": float(p11[kj]),
                     "window_start": str(ws.date()), "window_stop": str(we.date())}
        prepared.append({
            "i0": i0, "i1": i1, "m10": m10, "m11": m11,
            "t10": t10, "t11": t11,
            "basis10": (c10, s10, c210, s210),
            "basis11": (c11, s11, c211, s211),
        })

    if not prepared:
        return {"status": "insufficient_window_data", "simulations": 0,
                "common_calendar": {"start": str(common_start.date()), "stop": str(common_stop.date())}}

    ar10 = estimate_daily_ar1_from_grid(g10)
    ar11 = estimate_daily_ar1_from_grid(g11)
    phi10 = float(ar10["phi_1day"]); phi11 = float(ar11["phi_1day"])
    rng = np.random.default_rng(seed)
    null10: list[float] = []
    null11: list[float] = []
    nullj: list[float] = []
    total = int(simulations)
    bs = max(1, int(batch_size))
    ngrid = len(g10)
    sd10 = float(np.sqrt(max(1e-12, 1.0 - phi10 * phi10)))
    sd11 = float(np.sqrt(max(1e-12, 1.0 - phi11 * phi11)))

    for st in range(0, total, bs):
        b = min(bs, total - st)
        x10 = np.empty((b, ngrid), dtype=float)
        x11 = np.empty((b, ngrid), dtype=float)
        x10[:, 0] = rng.normal(size=b)
        x11[:, 0] = rng.normal(size=b)
        for k in range(1, ngrid):
            x10[:, k] = phi10 * x10[:, k - 1] + rng.normal(0.0, sd10, size=b)
            x11[:, k] = phi11 * x11[:, k - 1] + rng.normal(0.0, sd11, size=b)
        mx10 = np.full(b, -np.inf, dtype=float)
        mx11 = np.full(b, -np.inf, dtype=float)
        mxj = np.full(b, -np.inf, dtype=float)
        for w in prepared:
            y10 = x10[:, w["i0"]:w["i1"]][:, w["m10"]]
            y11 = x11[:, w["i0"]:w["i1"]][:, w["m11"]]
            y10 = _center_or_linear_detrend_batch(y10, w["t10"], detrend)
            y11 = _center_or_linear_detrend_batch(y11, w["t11"], detrend)
            p10 = _classic_lomb_power_batch(y10, *w["basis10"])
            p11 = _classic_lomb_power_batch(y11, *w["basis11"])
            joint = np.sqrt(np.maximum(p10, 0.0) * np.maximum(p11, 0.0))
            mx10 = np.maximum(mx10, np.nanmax(p10, axis=1))
            mx11 = np.maximum(mx11, np.nanmax(p11, axis=1))
            mxj = np.maximum(mxj, np.nanmax(joint, axis=1))
        null10.extend(mx10.tolist()); null11.extend(mx11.tolist()); nullj.extend(mxj.tolist())

    def summarize(null: list[float], observed: float) -> dict:
        arr = np.asarray(null, dtype=float)
        arr = arr[np.isfinite(arr)]
        if len(arr) == 0:
            return {"p_value": np.nan, "null_simulations_used": 0, "exceedances": 0}
        exc = int(np.sum(arr >= float(observed)))
        return {
            "p_value": float((exc + 1) / (len(arr) + 1)),
            "exceedances": exc,
            "null_simulations_used": int(len(arr)),
            "null_q95": float(np.quantile(arr, 0.95)),
            "null_q99": float(np.quantile(arr, 0.99)),
        }

    return {
        "status": "ok",
        "transform": "log10 positive flux",
        "detrend": detrend,
        "band_calendar_days": [float(lo), float(hi)],
        "frequency_grid_points": int(nfreq),
        "window_years": int(window_years),
        "step_years": int(step_years),
        "common_calendar": {"start": str(common_start.date()), "stop": str(common_stop.date()),
                            "grid_days": int(len(g10))},
        "windows_used": int(len(prepared)),
        "ar1": {"p10": ar10, "p11": ar11},
        "observed": {"p10_max_over_time_frequency": best10,
                     "p11_max_over_time_frequency": best11,
                     "joint_same_window_same_period_max": bestj},
        "observed_windows": observed_windows,
        "null": {"p10": summarize(null10, best10["power"]),
                 "p11": summarize(null11, best11["power"]),
                 "joint": summarize(nullj, bestj["power"])},
        "statistic_note": (
            "Each p-value is calibrated to the maximum across every predeclared window and every frequency in the locked band; "
            "the joint statistic is the maximum geometric mean of P10/P11 Lomb power at the same period in the same window."
        ),
    }
