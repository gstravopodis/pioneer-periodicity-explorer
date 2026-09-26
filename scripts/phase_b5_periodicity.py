"""Phase B5: individual-spacecraft periodicity morphology and time localization.

B4 established that a significant cross-spacecraft lag is not automatically a
periodicity. B5 therefore tests the *individual* P10 and P11 time series for
periodic structure using two interpolation-free views:

1. calendar-lag ACF morphology on the NaN-preserving daily grid;
2. Lomb–Scargle morphology on the actual observation times.

The same B3/B4 candidate bands are carried forward. Sliding-window spectra are
descriptive localization only; no new p-values are created here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.io import load_auto, apply_p11_cpi_quality_mask
from src.modern import (
    dense_log_grid,
    fast_calendar_dcf_from_grids,
    fast_bin_dcf,
    dcf_band_max,
    lomb_scargle_log_periodogram,
)


def _finite_span(df: pd.DataFrame, channel: str):
    q = df[["date", channel]].copy()
    q["date"] = pd.to_datetime(q["date"], errors="coerce").dt.normalize()
    q[channel] = pd.to_numeric(q[channel], errors="coerce")
    q = q.dropna(subset=["date", channel])
    q = q[q[channel] > 0]
    if q.empty:
        return None, None
    return pd.Timestamp(q.date.min()), pd.Timestamp(q.date.max())


def _subset(df: pd.DataFrame, start: pd.Timestamp, stop: pd.Timestamp) -> pd.DataFrame:
    d = pd.to_datetime(df["date"], errors="coerce").dt.normalize()
    return df.loc[(d >= start) & (d <= stop)].copy()


def _window_starts(start: pd.Timestamp, stop: pd.Timestamp, years: int, step_years: int):
    cur = pd.Timestamp(start).normalize()
    while True:
        end = cur + pd.DateOffset(years=int(years)) - pd.Timedelta(days=1)
        if end > stop:
            break
        yield cur, end
        cur = cur + pd.DateOffset(years=int(step_years))


def dcf_peak_morphology(binned: pd.DataFrame, lo: int, hi: int, min_pairs: int = 100) -> dict:
    peak = dcf_band_max(binned, lo, hi, min_pairs=min_pairs)
    if not np.isfinite(peak.get("dcf", np.nan)):
        return {**peak, "at_band_edge": None, "is_local_peak": None, "prominence": np.nan}
    q = binned[(binned.pairs >= int(min_pairs)) & np.isfinite(binned.dcf)].copy().reset_index(drop=True)
    if q.empty:
        return {**peak, "at_band_edge": None, "is_local_peak": None, "prominence": np.nan}
    y = q.dcf.to_numpy(float)
    peaks, props = signal.find_peaks(y, prominence=0)
    centers = q.lag_bin_center.to_numpy(float)
    target = float(peak["lag_calendar_days"])
    idx = int(np.argmin(np.abs(centers - target)))
    pos = np.where(peaks == idx)[0]
    prom = float(props["prominences"][pos[0]]) if len(pos) else 0.0
    width = float(q.lag_bin_end.iloc[0] - q.lag_bin_start.iloc[0] + 1)
    at_edge = target <= lo + width / 2 or target >= hi - width / 2
    return {**peak, "at_band_edge": bool(at_edge), "is_local_peak": bool(len(pos)), "prominence": prom}


def spectrum_peak_morphology(periodogram: pd.DataFrame, lo: int, hi: int) -> dict:
    q = periodogram[np.isfinite(periodogram.power)].copy().sort_values("period_days").reset_index(drop=True)
    b = q[(q.period_days >= float(lo)) & (q.period_days <= float(hi))]
    if b.empty:
        return {"period_days": np.nan, "power": np.nan, "at_band_edge": None,
                "is_local_peak": None, "prominence": np.nan}
    row = b.loc[b.power.idxmax()]
    target = float(row.period_days)
    y = q.power.to_numpy(float)
    peaks, props = signal.find_peaks(y, prominence=0)
    idx = int(np.argmin(np.abs(q.period_days.to_numpy(float) - target)))
    pos = np.where(peaks == idx)[0]
    prom = float(props["prominences"][pos[0]]) if len(pos) else 0.0
    # Frequency-grid spectra do not have constant period spacing, so edge is
    # defined as the outer 2% of the predeclared period band.
    margin = 0.02 * (hi - lo)
    at_edge = target <= lo + margin or target >= hi - margin
    return {
        "period_days": target,
        "power": float(row.power),
        "at_band_edge": bool(at_edge),
        "is_local_peak": bool(len(pos)),
        "prominence": prom,
    }


def top_spectral_peaks(periodogram: pd.DataFrame, n: int = 6) -> list[dict]:
    q = periodogram[np.isfinite(periodogram.power)].copy().sort_values("period_days").reset_index(drop=True)
    if len(q) < 3:
        return []
    y = q.power.to_numpy(float)
    peaks, props = signal.find_peaks(y, prominence=0)
    rows = []
    for j, idx in enumerate(peaks):
        rows.append({
            "period_days": float(q.loc[idx, "period_days"]),
            "power": float(q.loc[idx, "power"]),
            "prominence": float(props["prominences"][j]),
        })
    rows.sort(key=lambda r: (r["prominence"], r["power"]), reverse=True)
    return rows[: int(n)]


def _acf_profile(df: pd.DataFrame, channel: str, start, stop, detrend: str,
                 min_lag: int, max_lag: int, bin_days: int) -> pd.DataFrame:
    _, grid = dense_log_grid(df, channel, start=start, stop=stop, detrend=detrend)
    if len(grid) == 0:
        return pd.DataFrame(columns=["lag_bin_start", "lag_bin_end", "lag_bin_center", "pairs", "dcf"])
    dcf = fast_calendar_dcf_from_grids(grid, None, min_lag, max_lag)
    return fast_bin_dcf(dcf, bin_days=bin_days, origin=min_lag)


def _ls_profile(df: pd.DataFrame, channel: str, start, stop, detrend: str,
                min_period: int, max_period: int, nfreq: int) -> pd.DataFrame:
    sub = _subset(df, start, stop)
    return lomb_scargle_log_periodogram(
        sub, channel,
        min_period_days=float(min_period),
        max_period_days=float(max_period),
        nfreq=int(nfreq),
        detrend=detrend,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--b4", type=Path, default=ROOT / "results" / "phase_b4_localization" / "localization_summary.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b5_periodicity")
    ap.add_argument("--min-period", type=int, default=200)
    ap.add_argument("--max-period", type=int, default=700)
    ap.add_argument("--bin-days", type=int, default=10)
    ap.add_argument("--nfreq", type=int, default=2400)
    ap.add_argument("--window-nfreq", type=int, default=1200)
    ap.add_argument("--window-years", type=int, default=5)
    ap.add_argument("--step-years", type=int, default=1)
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A first.")
    if not args.b4.exists():
        raise SystemExit("B4 localization_summary.json is missing. Complete B4 first.")

    p10 = load_auto(args.p10, args.p10.name)
    p11 = apply_p11_cpi_quality_mask(load_auto(args.p11, args.p11.name))
    b4 = json.loads(args.b4.read_text(encoding="utf-8"))
    candidates = b4.get("candidates", [])
    if not candidates:
        raise SystemExit("No B4 candidates found.")

    bands = {"annual_330_400": (330, 400), "one_point_three_year_430_520": (430, 520)}
    args.out.mkdir(parents=True, exist_ok=True)
    profiles_dir = args.out / "profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "phase": "B5 individual-spacecraft periodicity morphology and time localization",
        "purpose": "Distinguish a cross-spacecraft lag correlation from periodic structure in the individual P10 and P11 records.",
        "guardrails": [
            "A CCF lag is not a period. Periodicity requires structure within each spacecraft time series.",
            "B5 adds morphology and localization, not new multiple-testing-corrected p-values.",
            "Sliding 5-year windows overlap and are descriptive; repeated adjacent windows are not independent replications.",
            "Lomb-Scargle uses actual irregular observation times; no missing value is interpolated.",
        ],
        "period_scan_calendar_days": [int(args.min_period), int(args.max_period)],
        "acf_bin_days": int(args.bin_days),
        "window_years": int(args.window_years),
        "step_years": int(args.step_years),
        "candidates": [],
    }
    flat_rows = []

    # Preserve one record per B4 channel+band candidate.
    for ci, c in enumerate(candidates, start=1):
        channel = c["channel"]
        band_name = c["band"]
        lo, hi = bands[band_name]
        print(f"[{ci}/{len(candidates)}] B5 {channel} {band_name}", flush=True)
        rec = {
            "channel": channel,
            "label": c.get("label", channel),
            "band": band_name,
            "p11_quality_limited_after_1980_239": bool(c.get("p11_quality_limited_after_1980_239", False)),
            "spacecraft": {},
            "common_window_period_concordance": [],
        }

        full_results = {}
        for sc_name, df in (("p10", p10), ("p11", p11)):
            start, stop = _finite_span(df, channel)
            if start is None:
                rec["spacecraft"][sc_name] = {"status": "insufficient_data"}
                continue
            sc_rec = {
                "calendar": {"start": str(start.date()), "stop": str(stop.date())},
                "variants": {},
            }
            full_results[sc_name] = {}
            for vname, detrend in (("raw_log10", "none"), ("linear_detrended_log10", "linear")):
                acf = _acf_profile(df, channel, start, stop, detrend,
                                   args.min_period, args.max_period, args.bin_days)
                ls = _ls_profile(df, channel, start, stop, detrend,
                                 args.min_period, args.max_period, args.nfreq)
                acf_name = f"{channel}__{band_name}__{sc_name}__{vname}__acf.csv"
                ls_name = f"{channel}__{band_name}__{sc_name}__{vname}__lomb_scargle.csv"
                acf.to_csv(profiles_dir / acf_name, index=False)
                ls.to_csv(profiles_dir / ls_name, index=False)
                full_results[sc_name][vname] = {"acf": acf, "ls": ls}
                sc_rec["variants"][vname] = {
                    "acf_band_morphology": dcf_peak_morphology(acf, lo, hi),
                    "lomb_scargle_band_morphology": spectrum_peak_morphology(ls, lo, hi),
                    "top_spectral_peaks_scan": top_spectral_peaks(ls),
                    "acf_profile_csv": str(Path("profiles") / acf_name),
                    "lomb_scargle_profile_csv": str(Path("profiles") / ls_name),
                    "windows": [],
                }
                for ws, we in _window_starts(start, stop, args.window_years, args.step_years):
                    wls = _ls_profile(df, channel, ws, we, detrend,
                                      args.min_period, args.max_period, args.window_nfreq)
                    sm = spectrum_peak_morphology(wls, lo, hi)
                    sc_rec["variants"][vname]["windows"].append({
                        "window_start": str(ws.date()),
                        "window_stop": str(we.date()),
                        **sm,
                    })
                    flat_rows.append({
                        "channel": channel, "band": band_name, "spacecraft": sc_name,
                        "variant": vname, "window_start": str(ws.date()), "window_stop": str(we.date()),
                        **sm,
                    })
            rec["spacecraft"][sc_name] = sc_rec

        # Compare periods only in genuinely common five-year windows and only descriptively.
        p10_start, p10_stop = _finite_span(p10, channel)
        p11_start, p11_stop = _finite_span(p11, channel)
        if p10_start is not None and p11_start is not None:
            cs = max(p10_start, p11_start)
            ce = min(p10_stop, p11_stop)
            if ce > cs:
                for ws, we in _window_starts(cs, ce, args.window_years, args.step_years):
                    row = {"window_start": str(ws.date()), "window_stop": str(we.date()), "variants": {}}
                    for vname, detrend in (("raw_log10", "none"), ("linear_detrended_log10", "linear")):
                        s10 = spectrum_peak_morphology(
                            _ls_profile(p10, channel, ws, we, detrend, args.min_period, args.max_period, args.window_nfreq), lo, hi)
                        s11 = spectrum_peak_morphology(
                            _ls_profile(p11, channel, ws, we, detrend, args.min_period, args.max_period, args.window_nfreq), lo, hi)
                        p10_period = s10.get("period_days", np.nan)
                        p11_period = s11.get("period_days", np.nan)
                        diff = abs(p10_period - p11_period) if np.isfinite(p10_period) and np.isfinite(p11_period) else np.nan
                        meanp = (p10_period + p11_period) / 2 if np.isfinite(diff) else np.nan
                        row["variants"][vname] = {
                            "p10_period_days": p10_period,
                            "p11_period_days": p11_period,
                            "absolute_difference_days": float(diff) if np.isfinite(diff) else np.nan,
                            "relative_difference": float(diff / meanp) if np.isfinite(meanp) and meanp else np.nan,
                            "within_10_percent": bool(diff / meanp <= 0.10) if np.isfinite(meanp) and meanp else None,
                        }
                    rec["common_window_period_concordance"].append(row)
        summary["candidates"].append(rec)

    (args.out / "periodicity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame(flat_rows).to_csv(args.out / "sliding_lomb_scargle_summary.csv", index=False)
    print(json.dumps(summary, indent=2))
    print(f"\nB5 complete. See {args.out / 'periodicity_summary.json'}")


if __name__ == "__main__":
    main()
