"""Phase B2: multichannel uneven-sampling screening.

Phase A remains frozen.  B2 adds a discrete correlation function (DCF) on true
calendar lags and a generalized Lomb–Scargle screen on log10 flux.  It is a
screening/falsification stage, not a final significance claim.  Expensive red-
noise Monte Carlo tests are intentionally reserved for candidates that survive
this stage.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.io import load_auto, apply_p11_cpi_quality_mask
from src.schema import PARTICLE_COLUMNS, CHANNEL_LABELS
from src.modern import (
    calendar_discrete_correlation,
    bin_discrete_correlation,
    strongest_dcf_bin,
    lomb_scargle_log_periodogram,
    band_peak,
    sampling_window_periodogram,
    strongest_period_in_band,
)

BANDS = {
    "annual_330_400": (330, 400),
    "one_point_three_year_430_520": (430, 520),
}


def _finite_count(df: pd.DataFrame, channel: str) -> int:
    return int(pd.to_numeric(df[channel], errors="coerce").notna().sum()) if channel in df.columns else 0


def _screen_channel(p10: pd.DataFrame, p11: pd.DataFrame, channel: str, out: Path,
                    bin_days: int, save_tables: bool) -> dict:
    dcf10_exact = calendar_discrete_correlation(p10, channel, None, 250, 600, transform="log10")
    dcf11_exact = calendar_discrete_correlation(p11, channel, None, 250, 600, transform="log10")
    dcf_cc_exact = calendar_discrete_correlation(p10, channel, p11, 250, 600, transform="log10")
    dcf10 = bin_discrete_correlation(dcf10_exact, bin_days=bin_days, origin=250)
    dcf11 = bin_discrete_correlation(dcf11_exact, bin_days=bin_days, origin=250)
    dcfcc = bin_discrete_correlation(dcf_cc_exact, bin_days=bin_days, origin=250)

    gls10 = lomb_scargle_log_periodogram(p10, channel, 250, 600, nfreq=3500)
    gls11 = lomb_scargle_log_periodogram(p11, channel, 250, 600, nfreq=3500)
    win10 = sampling_window_periodogram(p10, channel, 250, 600)
    win11 = sampling_window_periodogram(p11, channel, 250, 600)

    if save_tables:
        slug = channel
        dcf10.to_csv(out / f"p10_{slug}_dcf_log10.csv", index=False)
        dcf11.to_csv(out / f"p11_{slug}_dcf_log10.csv", index=False)
        dcfcc.to_csv(out / f"p10_p11_{slug}_dcf_log10.csv", index=False)
        gls10.to_csv(out / f"p10_{slug}_gls_log10.csv", index=False)
        gls11.to_csv(out / f"p11_{slug}_gls_log10.csv", index=False)

    bands = {}
    for name, (lo, hi) in BANDS.items():
        g10 = band_peak(gls10, lo, hi)
        g11 = band_peak(gls11, lo, hi)
        w10 = strongest_period_in_band(win10, lo, hi)
        w11 = strongest_period_in_band(win11, lo, hi)
        bands[name] = {
            "p10_dcf_acf": strongest_dcf_bin(dcf10, lo, hi),
            "p11_dcf_acf": strongest_dcf_bin(dcf11, lo, hi),
            "p10_p11_dcf_ccf": strongest_dcf_bin(dcfcc, lo, hi),
            "p10_gls_log10": g10,
            "p11_gls_log10": g11,
            "p10_sampling_window": w10,
            "p11_sampling_window": w11,
            "p10_gls_minus_window_days": (
                float(g10["period_days"] - w10["period_days"])
                if np.isfinite(g10["period_days"]) and np.isfinite(w10["period_days"]) else np.nan
            ),
            "p11_gls_minus_window_days": (
                float(g11["period_days"] - w11["period_days"])
                if np.isfinite(g11["period_days"]) and np.isfinite(w11["period_days"]) else np.nan
            ),
        }

    return {
        "channel": channel,
        "label": CHANNEL_LABELS.get(channel, channel),
        "finite_values": {"p10": _finite_count(p10, channel), "p11": _finite_count(p11, channel)},
        "dcf_definition": {
            "transform": "log10 positive flux",
            "bin_days": int(bin_days),
            "positive_cross_lag": "P10(t+L) paired with P11(t)",
            "interpolation": False,
        },
        "bands": bands,
    }


def _flatten(summary: dict) -> pd.DataFrame:
    rows = []
    for ch in summary["channels"]:
        for band_name, b in ch["bands"].items():
            rows.append({
                "channel": ch["channel"],
                "label": ch["label"],
                "band": band_name,
                "p10_n": ch["finite_values"]["p10"],
                "p11_n": ch["finite_values"]["p11"],
                "p10_dcf_lag": b["p10_dcf_acf"].get("lag_calendar_days"),
                "p10_dcf": b["p10_dcf_acf"].get("dcf"),
                "p11_dcf_lag": b["p11_dcf_acf"].get("lag_calendar_days"),
                "p11_dcf": b["p11_dcf_acf"].get("dcf"),
                "ccf_dcf_lag": b["p10_p11_dcf_ccf"].get("lag_calendar_days"),
                "ccf_dcf": b["p10_p11_dcf_ccf"].get("dcf"),
                "p10_gls_period": b["p10_gls_log10"].get("period_days"),
                "p10_gls_power": b["p10_gls_log10"].get("power"),
                "p10_window_period": b["p10_sampling_window"].get("period_days"),
                "p10_gls_window_delta": b["p10_gls_minus_window_days"],
                "p11_gls_period": b["p11_gls_log10"].get("period_days"),
                "p11_gls_power": b["p11_gls_log10"].get("power"),
                "p11_window_period": b["p11_sampling_window"].get("period_days"),
                "p11_gls_window_delta": b["p11_gls_minus_window_days"],
            })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--channel", default="p_11_20_mev")
    ap.add_argument("--all-channels", action="store_true")
    ap.add_argument("--bin-days", type=int, default=10)
    ap.add_argument("--save-tables", action="store_true")
    ap.add_argument("--no-p11-mask", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b2_screening")
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 daily CPI files are missing. Complete Phase A acquisition first.")
    args.out.mkdir(parents=True, exist_ok=True)
    p10 = load_auto(args.p10, args.p10.name)
    p11_raw = load_auto(args.p11, args.p11.name)
    p11 = p11_raw if args.no_p11_mask else apply_p11_cpi_quality_mask(p11_raw)

    channels = PARTICLE_COLUMNS if args.all_channels else [args.channel]
    channels = [c for c in channels if c in p10.columns and c in p11.columns]
    if not channels:
        raise SystemExit("No requested particle channels are present in both datasets.")

    out = {
        "phase": "B2 multichannel uneven-sampling screening",
        "purpose": (
            "Screen annual and ~1.3-year bands with an interpolation-free DCF and log10 generalized "
            "Lomb-Scargle before spending computation on red-noise significance tests."
        ),
        "guardrail": (
            "B2 is not a significance claim. Peaks close to the sampling-window spectrum remain suspect; "
            "only B3 red-noise/surrogate tests may promote a candidate."
        ),
        "bands_calendar_days": BANDS,
        "channels": [],
    }
    for i, channel in enumerate(channels, start=1):
        print(f"[{i}/{len(channels)}] screening {channel} ...", flush=True)
        out["channels"].append(_screen_channel(p10, p11, channel, args.out, args.bin_days, args.save_tables))

    json_path = args.out / "multichannel_screening_summary.json"
    csv_path = args.out / "multichannel_screening_summary.csv"
    json_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    flat = _flatten(out)
    flat.to_csv(csv_path, index=False)

    # Compact console ranking for the two predeclared bands.
    for band in BANDS:
        q = flat[flat.band == band].copy()
        q["screen_score"] = q[["p10_dcf", "p11_dcf", "ccf_dcf"]].max(axis=1, skipna=True)
        q = q.sort_values("screen_score", ascending=False)
        print(f"\n{band} — top DCF screen values")
        cols = ["channel", "p10_dcf_lag", "p10_dcf", "p11_dcf_lag", "p11_dcf", "ccf_dcf_lag", "ccf_dcf"]
        print(q[cols].head(8).to_string(index=False))

    print(f"\nPhase B2 complete. Full JSON: {json_path}")
    print(f"Flat CSV: {csv_path}")


if __name__ == "__main__":
    main()
