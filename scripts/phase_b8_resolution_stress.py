"""Phase B8: multi-scale Rayleigh-resolution stress test.

B7 used overlapping 5-year windows. At periods of a few hundred days, one
Rayleigh frequency resolution for a five-year segment corresponds to a very
wide interval in period space. Consequently, a <=1-Rayleigh proximity warning
can cover much of the predeclared 330-400 d or 430-520 d bands and is not by
itself a discriminating alias diagnostic.

B8 repeats the signal-vs-sampling-window comparison at progressively longer
common-calendar windows. Longer windows narrow the Rayleigh frequency width;
features that remain unresolved from the sampling window at 10-15 yr are more
concerning than features flagged only in five-year windows.

This phase is descriptive and creates no new p-values.
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


def sampling_window_periodogram_positive(df: pd.DataFrame, channel: str,
                                         start: pd.Timestamp, stop: pd.Timestamp,
                                         min_period: int = 200, max_period: int = 700) -> pd.DataFrame:
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


def rayleigh_period_context(signal_period: float, window_period: float, span_days: int) -> dict:
    """Express signal/window separation and local Rayleigh width in frequency units.

    ``approx_period_resolution_days`` is the first-order conversion P^2/T and is
    included to make clear how broad a one-Rayleigh interval is in period space.
    """
    if not (np.isfinite(signal_period) and signal_period > 0 and
            np.isfinite(window_period) and window_period > 0 and span_days > 0):
        return {
            "absolute_period_difference_days": np.nan,
            "frequency_separation_rayleigh": np.nan,
            "within_one_rayleigh": None,
            "approx_period_resolution_days": np.nan,
        }
    p = float(signal_period)
    t = float(span_days)
    sep = abs((1.0 / p) - (1.0 / float(window_period))) * t
    return {
        "absolute_period_difference_days": abs(p - float(window_period)),
        "frequency_separation_rayleigh": float(sep),
        "within_one_rayleigh": bool(sep <= 1.0),
        "approx_period_resolution_days": float((p * p) / t),
    }


def annual_windows(start: pd.Timestamp, stop: pd.Timestamp, years: int, step_years: int = 1):
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


def _window_result(df: pd.DataFrame, channel: str, start: pd.Timestamp, stop: pd.Timestamp,
                   detrend: str, lo: int, hi: int, min_period: int, max_period: int,
                   nfreq: int) -> dict:
    sub = _subset(df, start, stop)
    ls = lomb_scargle_log_periodogram(
        sub, channel, min_period_days=float(min_period), max_period_days=float(max_period),
        nfreq=int(nfreq), detrend=detrend,
    )
    win = sampling_window_periodogram_positive(df, channel, start, stop, min_period, max_period)
    sm = spectrum_peak_morphology(ls, lo, hi)
    wm = spectrum_peak_morphology(win, lo, hi)
    span = int((stop - start).days + 1)
    return {
        "window_start": str(start.date()),
        "window_stop": str(stop.date()),
        "span_days": span,
        "signal": sm,
        "sampling_window": wm,
        "signal_vs_sampling_window": rayleigh_period_context(
            sm.get("period_days", np.nan), wm.get("period_days", np.nan), span
        ),
    }


def _scale_summary(results: list[dict]) -> dict:
    usable = [r for r in results if r["signal"].get("is_local_peak") is True and
              np.isfinite(r["signal"].get("period_days", np.nan)) and
              np.isfinite(r["sampling_window"].get("period_days", np.nan))]
    if not usable:
        return {
            "windows_total": len(results), "local_peak_windows": 0,
            "alias_fraction": np.nan, "median_rayleigh": np.nan,
            "median_period_difference_days": np.nan,
            "median_approx_period_resolution_days": np.nan,
        }
    alias = [bool(r["signal_vs_sampling_window"]["within_one_rayleigh"]) for r in usable]
    rays = [float(r["signal_vs_sampling_window"]["frequency_separation_rayleigh"]) for r in usable]
    pdiff = [float(r["signal_vs_sampling_window"]["absolute_period_difference_days"]) for r in usable]
    pres = [float(r["signal_vs_sampling_window"]["approx_period_resolution_days"]) for r in usable]
    return {
        "windows_total": int(len(results)),
        "local_peak_windows": int(len(usable)),
        "alias_warning_windows": int(sum(alias)),
        "alias_fraction": float(np.mean(alias)),
        "median_rayleigh": float(np.median(rays)),
        "median_period_difference_days": float(np.median(pdiff)),
        "median_approx_period_resolution_days": float(np.median(pres)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b8_resolution_stress")
    ap.add_argument("--channels", nargs="*", default=None)
    ap.add_argument("--window-years", nargs="+", type=int, default=[5, 8, 10, 12, 15])
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
    scales = sorted(set(int(x) for x in args.window_years if int(x) >= 2))

    out = {
        "phase": "B8 multi-scale Rayleigh-resolution stress test",
        "purpose": (
            "Determine whether B7's universal <=1-Rayleigh warnings were genuinely informative or largely a consequence "
            "of the coarse frequency resolution of five-year windows."
        ),
        "guardrails": [
            "This phase is descriptive and creates no new p-values.",
            "A five-year Rayleigh interval is broad in period space; a <=1-Rayleigh warning can cover much of a 330-520 d band.",
            "Longer windows improve frequency resolution but reduce time localization and the number of available windows.",
            "A persistent <=1-Rayleigh proximity at longer scales strengthens sampling concern; separation >1 Rayleigh weakens that specific concern but does not prove a physical origin.",
            "Channels limited by the Pioneer 11 post-1980 reliability mask may not support the longer window lengths.",
        ],
        "bands": BANDS,
        "window_years": scales,
        "step_years": int(args.step_years),
        "channels": [],
    }
    flat = []

    for ci, channel in enumerate(channels, start=1):
        print(f"[{ci}/{len(channels)}] B8 {channel}", flush=True)
        s10, e10 = _dated_positive_span(p10, channel)
        s11, e11 = _dated_positive_span(p11, channel)
        if s10 is None or s11 is None:
            out["channels"].append({"channel": channel, "status": "insufficient_data"})
            continue
        start, stop = max(s10, s11), min(e10, e11)
        common_span_days = int((stop - start).days + 1)
        rec = {
            "channel": channel,
            "label": CHANNEL_LABELS.get(channel, channel),
            "p11_quality_limited_after_1980_239": bool(e11 < pd.Timestamp("1981-01-01")),
            "common_calendar": {"start": str(start.date()), "stop": str(stop.date()), "span_days": common_span_days},
            "bands": {},
        }
        for band_name, (lo, hi) in BANDS.items():
            bnode = {"variants": {}}
            for variant, detrend in variants:
                vnode = {"scales": {}}
                for years in scales:
                    wins = annual_windows(start, stop, years=years, step_years=args.step_years)
                    scale_node = {"window_years": years, "windows_total": len(wins), "p10": [], "p11": []}
                    for ws, we in wins:
                        r10 = _window_result(p10, channel, ws, we, detrend, lo, hi,
                                             args.min_period, args.max_period, args.nfreq)
                        r11 = _window_result(p11, channel, ws, we, detrend, lo, hi,
                                             args.min_period, args.max_period, args.nfreq)
                        scale_node["p10"].append(r10)
                        scale_node["p11"].append(r11)
                    scale_node["p10_summary"] = _scale_summary(scale_node["p10"])
                    scale_node["p11_summary"] = _scale_summary(scale_node["p11"])
                    vnode["scales"][str(years)] = scale_node
                    for sc in ("p10", "p11"):
                        ss = scale_node[f"{sc}_summary"]
                        flat.append({
                            "channel": channel, "label": CHANNEL_LABELS.get(channel, channel),
                            "band": band_name, "variant": variant, "spacecraft": sc,
                            "window_years": years, "common_span_days": common_span_days,
                            **ss,
                        })
                bnode["variants"][variant] = vnode
            rec["bands"][band_name] = bnode
        out["channels"].append(rec)

    out_path = args.out / "resolution_stress_summary.json"
    csv_path = args.out / "resolution_stress_flags.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    pd.DataFrame(flat).to_csv(csv_path, index=False)

    compact = []
    for row in flat:
        if row["window_years"] in scales:
            compact.append(row)
    flags_path = args.out / "resolution_stress_flags.json"
    flags_path.write_text(json.dumps({
        "phase": "B8 resolution-stress flags (descriptive, no ranking)",
        "interpretation": [
            "Compare the same channel/band across window lengths; do not compare raw alias fractions without considering Rayleigh width.",
            "Five-year alias_fraction=1.0 is weak evidence when the median approximate period resolution is comparable to the band width.",
            "Long-window separation >1 Rayleigh is evidence against that specific sampling-window coincidence, not proof of periodicity.",
        ],
        "rows": compact,
    }, indent=2), encoding="utf-8")

    print(f"\nPhase B8 complete: {out_path}")
    print(f"Flags JSON: {flags_path}")
    print(f"Flat table: {csv_path}")


if __name__ == "__main__":
    main()
