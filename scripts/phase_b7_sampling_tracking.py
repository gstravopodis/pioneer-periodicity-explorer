"""Phase B7: sliding-window signal-vs-sampling co-tracking audit.

B6 showed that every publication-facing candidate still had at least one
one-Rayleigh sampling-window warning.  A single full-record proximity warning is
not enough to decide whether a spectral feature is an alias.  B7 asks a more
specific descriptive question: as the observation window changes through time,
does the apparent signal period move with the sampling-window peak?

For every CPI channel, spacecraft, robustness variant and predeclared band, B7
uses the exact common valid calendar span and overlapping 5-year windows.  In
each window it measures:

* the strongest Lomb-Scargle feature in the band;
* whether that feature is a genuine local maximum rather than a boundary max;
* the strongest sampling-window spectral feature in the same band;
* frequency separation in Rayleigh-resolution units;
* P10/P11 signal-period concordance on the same window.

The analysis is descriptive.  Overlapping windows are not independent tests,
and co-tracking with the sampling window is evidence of sampling sensitivity,
not proof that a feature is instrumental or non-physical.
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
from src.schema import PARTICLE_COLUMNS, CHANNEL_LABELS
from src.modern import lomb_scargle_log_periodogram, dense_log_grid

BANDS = {
    "annual_330_400": (330, 400),
    "one_point_three_year_430_520": (430, 520),
}


def _dated_positive_span(df: pd.DataFrame, channel: str):
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
    margin = 0.02 * (hi - lo)
    at_edge = target <= lo + margin or target >= hi - margin
    return {
        "period_days": target,
        "power": float(row.power),
        "at_band_edge": bool(at_edge),
        "is_local_peak": bool(len(pos)),
        "prominence": prom,
    }


def sampling_window_periodogram_positive(df: pd.DataFrame, channel: str, start: pd.Timestamp, stop: pd.Timestamp,
                                         min_period: int = 200, max_period: int = 700) -> pd.DataFrame:
    """FFT spectrum of the exact positive-flux observation mask in one window."""
    _, grid = dense_log_grid(df, channel, start=start, stop=stop, detrend="none")
    if len(grid) < 10:
        return pd.DataFrame(columns=["period_days", "power"])
    mask = np.isfinite(grid).astype(float)
    if mask.sum() < 10 or np.all(mask == mask[0]):
        return pd.DataFrame(columns=["period_days", "power"])
    y = mask - float(mask.mean())
    freq, power = signal.periodogram(y, fs=1.0, detrend=False, scaling="spectrum")
    keep = freq > 0
    freq, power = freq[keep], power[keep]
    period = 1.0 / freq
    keep = (period >= float(min_period)) & (period <= float(max_period))
    return pd.DataFrame({"period_days": period[keep], "power": power[keep]}).sort_values("period_days").reset_index(drop=True)


def rayleigh_separation(signal_period: float, window_period: float, span_days: int) -> dict:
    if not (np.isfinite(signal_period) and signal_period > 0 and np.isfinite(window_period) and window_period > 0 and span_days > 0):
        return {"absolute_period_difference_days": np.nan, "frequency_separation_rayleigh": np.nan,
                "within_one_rayleigh": None}
    ray = abs((1.0 / float(signal_period)) - (1.0 / float(window_period))) * float(span_days)
    return {
        "absolute_period_difference_days": abs(float(signal_period) - float(window_period)),
        "frequency_separation_rayleigh": float(ray),
        "within_one_rayleigh": bool(ray <= 1.0),
    }


def annual_windows(start: pd.Timestamp, stop: pd.Timestamp, years: int = 5, step_years: int = 1):
    cur = pd.Timestamp(start).normalize()
    stop = pd.Timestamp(stop).normalize()
    out = []
    while True:
        end = cur + pd.DateOffset(years=int(years)) - pd.Timedelta(days=1)
        if end > stop:
            break
        out.append((cur, end))
        cur = cur + pd.DateOffset(years=int(step_years))
    return out


def _window_spectrum(df: pd.DataFrame, channel: str, start: pd.Timestamp, stop: pd.Timestamp,
                     detrend: str, lo: int, hi: int, min_period: int, max_period: int, nfreq: int) -> dict:
    sub = _subset(df, start, stop)
    ls = lomb_scargle_log_periodogram(
        sub, channel, min_period_days=float(min_period), max_period_days=float(max_period),
        nfreq=int(nfreq), detrend=detrend,
    )
    win = sampling_window_periodogram_positive(df, channel, start, stop, min_period, max_period)
    sm = spectrum_peak_morphology(ls, lo, hi)
    wm = spectrum_peak_morphology(win, lo, hi)
    span = int((stop - start).days + 1)
    sep = rayleigh_separation(sm.get("period_days", np.nan), wm.get("period_days", np.nan), span)
    return {
        "signal": sm,
        "sampling_window": wm,
        "signal_vs_sampling_window": sep,
    }


def tracking_summary(windows: list[dict]) -> dict:
    usable = [w for w in windows if w.get("signal", {}).get("is_local_peak") is True and
              np.isfinite(w.get("signal", {}).get("period_days", np.nan)) and
              np.isfinite(w.get("sampling_window", {}).get("period_days", np.nan))]
    rays = [float(w["signal_vs_sampling_window"]["frequency_separation_rayleigh"])
            for w in usable if np.isfinite(w["signal_vs_sampling_window"].get("frequency_separation_rayleigh", np.nan))]
    aliases = [bool(w["signal_vs_sampling_window"]["within_one_rayleigh"])
               for w in usable if w["signal_vs_sampling_window"].get("within_one_rayleigh") is not None]
    sigf = np.asarray([1.0 / float(w["signal"]["period_days"]) for w in usable], dtype=float)
    winf = np.asarray([1.0 / float(w["sampling_window"]["period_days"]) for w in usable], dtype=float)
    if len(sigf) >= 4 and np.std(sigf) > 0 and np.std(winf) > 0:
        corr = float(np.corrcoef(sigf, winf)[0, 1])
    else:
        corr = np.nan
    return {
        "windows_total": int(len(windows)),
        "signal_local_peak_windows": int(len(usable)),
        "alias_warning_windows": int(sum(aliases)) if aliases else 0,
        "alias_fraction_among_local_peak_windows": float(np.mean(aliases)) if aliases else np.nan,
        "median_frequency_separation_rayleigh": float(np.median(rays)) if rays else np.nan,
        "signal_window_frequency_correlation": corr,
        "note": "High co-tracking or frequent <=1-Rayleigh proximity supports sampling sensitivity; neither is proof of aliasing.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b7_sampling_tracking")
    ap.add_argument("--channels", nargs="*", default=None)
    ap.add_argument("--window-years", type=int, default=5)
    ap.add_argument("--step-years", type=int, default=1)
    ap.add_argument("--min-period", type=int, default=200)
    ap.add_argument("--max-period", type=int, default=700)
    ap.add_argument("--nfreq", type=int, default=3500)
    ap.add_argument("--no-p11-mask", action="store_true")
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A acquisition first.")
    args.out.mkdir(parents=True, exist_ok=True)
    p10 = load_auto(args.p10, args.p10.name)
    p11raw = load_auto(args.p11, args.p11.name)
    p11 = p11raw if args.no_p11_mask else apply_p11_cpi_quality_mask(p11raw)

    requested = list(args.channels) if args.channels else list(PARTICLE_COLUMNS)
    channels = [c for c in requested if c in p10.columns and c in p11.columns]
    variants = [("raw_log10", "none"), ("linear_detrended_log10", "linear")]

    out = {
        "phase": "B7 sliding-window signal-vs-sampling co-tracking audit",
        "purpose": (
            "Test whether apparent periodicity peaks move with the observation-window spectrum through time, "
            "after B6 found at least one sampling-window proximity warning for every candidate."
        ),
        "guardrails": [
            "This phase is descriptive and creates no new significance p-values.",
            "Five-year windows overlap and are not independent replications.",
            "A <=1-Rayleigh signal/window separation is a warning, not proof of aliasing.",
            "High temporal co-tracking with the sampling-window frequency supports sampling sensitivity; low co-tracking does not prove a physical origin.",
            "P10 and P11 are compared on the exact same valid calendar windows for each channel.",
        ],
        "bands": BANDS,
        "window_years": int(args.window_years),
        "step_years": int(args.step_years),
        "channels": [],
    }
    flat = []

    for ci, channel in enumerate(channels, start=1):
        print(f"[{ci}/{len(channels)}] B7 {channel}", flush=True)
        s10, e10 = _dated_positive_span(p10, channel)
        s11, e11 = _dated_positive_span(p11, channel)
        if s10 is None or s11 is None:
            out["channels"].append({"channel": channel, "status": "insufficient_data"})
            continue
        start, stop = max(s10, s11), min(e10, e11)
        wins = annual_windows(start, stop, args.window_years, args.step_years)
        rec = {
            "channel": channel,
            "label": CHANNEL_LABELS.get(channel, channel),
            "p11_quality_limited_after_1980_239": bool(e11 < pd.Timestamp("1981-01-01")),
            "common_calendar": {"start": str(start.date()), "stop": str(stop.date()),
                                "span_days": int((stop - start).days + 1)},
            "bands": {},
        }
        for band_name, (lo, hi) in BANDS.items():
            bnode = {"variants": {}}
            for variant, detrend in variants:
                vwindows = []
                for ws, we in wins:
                    r10 = _window_spectrum(p10, channel, ws, we, detrend, lo, hi,
                                           args.min_period, args.max_period, args.nfreq)
                    r11 = _window_spectrum(p11, channel, ws, we, detrend, lo, hi,
                                           args.min_period, args.max_period, args.nfreq)
                    p10p = float(r10["signal"].get("period_days", np.nan))
                    p11p = float(r11["signal"].get("period_days", np.nan))
                    if np.isfinite(p10p) and np.isfinite(p11p) and max(p10p, p11p) > 0:
                        absdiff = abs(p10p - p11p)
                        reldiff = absdiff / ((p10p + p11p) / 2.0)
                    else:
                        absdiff = reldiff = np.nan
                    row = {
                        "window_start": str(ws.date()), "window_stop": str(we.date()),
                        "p10": r10, "p11": r11,
                        "p10_p11_signal_period_difference_days": float(absdiff) if np.isfinite(absdiff) else np.nan,
                        "p10_p11_signal_period_relative_difference": float(reldiff) if np.isfinite(reldiff) else np.nan,
                        "p10_p11_signal_periods_within_10_percent": bool(reldiff <= 0.10) if np.isfinite(reldiff) else None,
                    }
                    vwindows.append(row)
                    for sc, rr in (("p10", r10), ("p11", r11)):
                        flat.append({
                            "channel": channel, "label": CHANNEL_LABELS.get(channel, channel),
                            "band": band_name, "variant": variant, "spacecraft": sc,
                            "window_start": str(ws.date()), "window_stop": str(we.date()),
                            "signal_period_days": rr["signal"].get("period_days"),
                            "signal_power": rr["signal"].get("power"),
                            "signal_local_peak": rr["signal"].get("is_local_peak"),
                            "sampling_window_period_days": rr["sampling_window"].get("period_days"),
                            "window_period_power": rr["sampling_window"].get("power"),
                            "frequency_separation_rayleigh": rr["signal_vs_sampling_window"].get("frequency_separation_rayleigh"),
                            "sampling_alias_warning": rr["signal_vs_sampling_window"].get("within_one_rayleigh"),
                        })
                p10_track = tracking_summary([{**w["p10"]} for w in vwindows])
                p11_track = tracking_summary([{**w["p11"]} for w in vwindows])
                concordant = [w for w in vwindows if w.get("p10_p11_signal_periods_within_10_percent") is True and
                              w["p10"]["signal"].get("is_local_peak") is True and
                              w["p11"]["signal"].get("is_local_peak") is True]
                bnode["variants"][variant] = {
                    "p10_tracking_summary": p10_track,
                    "p11_tracking_summary": p11_track,
                    "cross_spacecraft_concordant_local_peak_windows": int(len(concordant)),
                    "windows_total": int(len(vwindows)),
                    "windows": vwindows,
                }
            rec["bands"][band_name] = bnode
        out["channels"].append(rec)

    out_path = args.out / "sampling_tracking_summary.json"
    csv_path = args.out / "sampling_tracking_windows.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    pd.DataFrame(flat).to_csv(csv_path, index=False)

    # Compact audit flags for rapid review; deliberately no score or ranking.
    flags = []
    for rec in out["channels"]:
        if rec.get("status"):
            continue
        for band_name, bnode in rec["bands"].items():
            row = {"channel": rec["channel"], "label": rec["label"], "band": band_name,
                   "p11_quality_limited_after_1980_239": rec["p11_quality_limited_after_1980_239"]}
            for variant in ("raw_log10", "linear_detrended_log10"):
                v = bnode["variants"][variant]
                for sc in ("p10", "p11"):
                    ts = v[f"{sc}_tracking_summary"]
                    row[f"{variant}__{sc}_local_peak_windows"] = ts["signal_local_peak_windows"]
                    row[f"{variant}__{sc}_alias_fraction"] = ts["alias_fraction_among_local_peak_windows"]
                    row[f"{variant}__{sc}_median_rayleigh"] = ts["median_frequency_separation_rayleigh"]
                    row[f"{variant}__{sc}_signal_window_frequency_corr"] = ts["signal_window_frequency_correlation"]
                row[f"{variant}__cross_spacecraft_concordant_local_peak_windows"] = v["cross_spacecraft_concordant_local_peak_windows"]
                row[f"{variant}__windows_total"] = v["windows_total"]
            flags.append(row)
    flags_path = args.out / "sampling_tracking_flags.json"
    flags_path.write_text(json.dumps({
        "phase": "B7 sampling-tracking flags (descriptive, no ranking)",
        "interpretation": [
            "No field in this review is a significance p-value.",
            "A high alias fraction or positive signal/window frequency co-tracking supports sampling sensitivity.",
            "Overlapping five-year windows are descriptive and must not be counted as independent replications.",
        ],
        "candidates": flags,
    }, indent=2), encoding="utf-8")

    print(f"\nPhase B7 complete: {out_path}")
    print(f"Window table: {csv_path}")
    print(f"Tracking flags: {flags_path}")


if __name__ == "__main__":
    main()
