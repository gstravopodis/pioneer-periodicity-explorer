"""Phase B4: full-lag morphology and time localization of B3 candidates.

B3 establishes whether correlation inside predeclared bands is unusual under the
chosen nulls. B4 asks a different question: does a surviving band maximum look
like a distinct lag feature, and is it persistent or confined to particular
calendar intervals?

No new global significance claim is made here. The phase is descriptive and is
specifically designed to prevent a band-edge maximum from being mislabeled as a
periodicity when it could simply be part of a broad/monotonic correlation tail.
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
from src.modern import dense_log_grid, fast_calendar_dcf_from_grids, fast_bin_dcf, dcf_band_max


def _positive_dates(df: pd.DataFrame, channel: str) -> pd.Series:
    q = df[["date", channel]].copy()
    q["date"] = pd.to_datetime(q["date"], errors="coerce").dt.normalize()
    q[channel] = pd.to_numeric(q[channel], errors="coerce")
    q = q.dropna(subset=["date", channel])
    q = q[q[channel] > 0]
    return q["date"]


def _common_span(df10: pd.DataFrame, df11: pd.DataFrame, channel: str):
    d10 = _positive_dates(df10, channel)
    d11 = _positive_dates(df11, channel)
    if d10.empty or d11.empty:
        return None, None
    start = max(d10.min(), d11.min())
    stop = min(d10.max(), d11.max())
    if stop <= start:
        return None, None
    return pd.Timestamp(start), pd.Timestamp(stop)


def _binned_ccf(df10: pd.DataFrame, df11: pd.DataFrame, channel: str, *,
                start: pd.Timestamp, stop: pd.Timestamp, detrend: str,
                min_lag: int, max_lag: int, bin_days: int) -> pd.DataFrame:
    _, g10 = dense_log_grid(df10, channel, start=start, stop=stop, detrend=detrend)
    _, g11 = dense_log_grid(df11, channel, start=start, stop=stop, detrend=detrend)
    dcf = fast_calendar_dcf_from_grids(g10, g11, min_lag, max_lag)
    return fast_bin_dcf(dcf, bin_days=bin_days, origin=min_lag)


def peak_morphology(binned: pd.DataFrame, lo: int, hi: int, min_pairs: int = 100) -> dict:
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
    width = int(round(float(q.lag_bin_end.iloc[0] - q.lag_bin_start.iloc[0] + 1)))
    at_edge = target <= lo + width / 2 or target >= hi - width / 2
    return {
        **peak,
        "at_band_edge": bool(at_edge),
        "is_local_peak": bool(len(pos)),
        "prominence": prom,
    }


def top_prominent_peaks(binned: pd.DataFrame, n: int = 6, min_pairs: int = 100) -> list[dict]:
    q = binned[(binned.pairs >= int(min_pairs)) & np.isfinite(binned.dcf)].copy().reset_index(drop=True)
    if len(q) < 3:
        return []
    y = q.dcf.to_numpy(float)
    peaks, props = signal.find_peaks(y, prominence=0)
    rows = []
    for j, idx in enumerate(peaks):
        rows.append({
            "lag_calendar_days": float(q.loc[idx, "lag_bin_center"]),
            "dcf": float(q.loc[idx, "dcf"]),
            "pairs": int(q.loc[idx, "pairs"]),
            "prominence": float(props["prominences"][j]),
        })
    rows.sort(key=lambda r: (r["prominence"], r["dcf"]), reverse=True)
    return rows[: int(n)]


def _window_starts(start: pd.Timestamp, stop: pd.Timestamp, years: int, step_years: int):
    cur = pd.Timestamp(start).normalize()
    while True:
        end = cur + pd.DateOffset(years=int(years)) - pd.Timedelta(days=1)
        if end > stop:
            break
        yield cur, end
        cur = cur + pd.DateOffset(years=int(step_years))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--review", type=Path, default=ROOT / "results" / "phase_b3_publication" / "candidate_review.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b4_localization")
    ap.add_argument("--min-lag", type=int, default=200)
    ap.add_argument("--max-lag", type=int, default=700)
    ap.add_argument("--bin-days", type=int, default=10)
    ap.add_argument("--window-years", type=int, default=5)
    ap.add_argument("--step-years", type=int, default=1)
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A first.")
    if not args.review.exists():
        raise SystemExit("Publication-resolution candidate_review.json is missing. Complete B3 publication run first.")

    p10 = load_auto(args.p10, args.p10.name)
    p11 = apply_p11_cpi_quality_mask(load_auto(args.p11, args.p11.name))
    review = json.loads(args.review.read_text(encoding="utf-8"))
    candidates = review.get("robust_cross_spacecraft_fdr05_both_variants", [])
    if not candidates:
        raise SystemExit("No robust cross-spacecraft candidates were found in the B3 review.")

    args.out.mkdir(parents=True, exist_ok=True)
    profiles_dir = args.out / "profiles"
    profiles_dir.mkdir(parents=True, exist_ok=True)

    bands = {"annual_330_400": (330, 400), "one_point_three_year_430_520": (430, 520)}
    unique = []
    seen = set()
    for c in candidates:
        key = (c["channel"], c["band"])
        if key not in seen:
            unique.append(c)
            seen.add(key)

    summary = {
        "phase": "B4 full-lag morphology and time localization",
        "purpose": "Characterize B3 survivors without turning a band maximum into a periodicity claim.",
        "guardrails": [
            "B3 significance means correlation inside a predeclared band is unusual under the chosen null; it does not prove an oscillatory period.",
            "A maximum at a band boundary is a warning that the band may be sampling a broad correlation tail rather than a distinct peak.",
            "Sliding-window results are descriptive localization only and are not separately multiple-testing corrected.",
        ],
        "lag_scan_calendar_days": [int(args.min_lag), int(args.max_lag)],
        "bin_days": int(args.bin_days),
        "window_years": int(args.window_years),
        "step_years": int(args.step_years),
        "candidates": [],
    }
    window_rows = []

    for ci, c in enumerate(unique, start=1):
        channel = c["channel"]
        band_name = c["band"]
        lo, hi = bands[band_name]
        start, stop = _common_span(p10, p11, channel)
        print(f"[{ci}/{len(unique)}] B4 {channel} {band_name}", flush=True)
        if start is None:
            continue
        rec = {
            "channel": channel,
            "label": c.get("label", channel),
            "band": band_name,
            "p11_quality_limited_after_1980_239": bool(c.get("p11_quality_limited_after_1980_239", False)),
            "common_calendar": {"start": str(start.date()), "stop": str(stop.date())},
            "variants": {},
        }
        for vname, detrend in (("raw_log10", "none"), ("linear_detrended_log10", "linear")):
            b = _binned_ccf(p10, p11, channel, start=start, stop=stop, detrend=detrend,
                            min_lag=args.min_lag, max_lag=args.max_lag, bin_days=args.bin_days)
            profile_name = f"{channel}__{band_name}__{vname}.csv"
            b.to_csv(profiles_dir / profile_name, index=False)
            rec["variants"][vname] = {
                "band_peak_morphology": peak_morphology(b, lo, hi),
                "top_prominent_peaks_scan": top_prominent_peaks(b),
                "profile_csv": str(Path("profiles") / profile_name),
                "windows": [],
            }
            for ws, we in _window_starts(start, stop, args.window_years, args.step_years):
                wb = _binned_ccf(p10, p11, channel, start=ws, stop=we, detrend=detrend,
                                 min_lag=args.min_lag, max_lag=args.max_lag, bin_days=args.bin_days)
                pk = peak_morphology(wb, lo, hi, min_pairs=30)
                wr = {
                    "window_start": str(ws.date()), "window_stop": str(we.date()),
                    "lag_calendar_days": pk.get("lag_calendar_days"),
                    "dcf": pk.get("dcf"), "pairs": pk.get("pairs"),
                    "at_band_edge": pk.get("at_band_edge"),
                    "is_local_peak": pk.get("is_local_peak"),
                    "prominence": pk.get("prominence"),
                }
                rec["variants"][vname]["windows"].append(wr)
                window_rows.append({"channel": channel, "band": band_name, "variant": vname, **wr})
        summary["candidates"].append(rec)

    (args.out / "localization_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame(window_rows).to_csv(args.out / "sliding_window_summary.csv", index=False)
    print(json.dumps(summary, indent=2))
    print(f"\nB4 complete. See {args.out / 'localization_summary.json'}")


if __name__ == "__main__":
    main()
