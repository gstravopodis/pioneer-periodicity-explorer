"""Phase B1: calendar-aware falsification of the historical lag interpretation.

This script intentionally leaves the successful Phase-A replication untouched.
It asks a different question: what happens when lags are measured in actual
calendar days rather than row/sample indices after missing dates have been
removed?
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.io import load_auto, apply_p11_cpi_quality_mask
from src.replication import calendar_separation_for_single_sample_lag
from src.modern import (
    calendar_autocorrelation,
    calendar_crosscorrelation,
    strongest_calendar_lag,
    calendar_lag_exact,
    lomb_scargle_periodogram,
    band_peak,
    sampling_window_periodogram,
    strongest_period_in_band,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--channel", default="p_11_20_mev")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b_calendar")
    ap.add_argument("--max-lag", type=int, default=2500)
    ap.add_argument("--no-p11-mask", action="store_true")
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 daily CPI files are missing. Complete Phase A acquisition first.")
    args.out.mkdir(parents=True, exist_ok=True)

    p10 = load_auto(args.p10, args.p10.name)
    p11_raw = load_auto(args.p11, args.p11.name)
    p11 = p11_raw if args.no_p11_mask else apply_p11_cpi_quality_mask(p11_raw)

    for name, df in (("p10", p10), ("p11", p11)):
        if args.channel not in df.columns:
            raise SystemExit(f"{args.channel} not found in {name}")

    ac10 = calendar_autocorrelation(p10, args.channel, 1, args.max_lag)
    ac11 = calendar_autocorrelation(p11, args.channel, 1, args.max_lag)
    cc = calendar_crosscorrelation(p10, p11, args.channel, -args.max_lag, args.max_lag)
    ac10.to_csv(args.out / f"p10_{args.channel}_calendar_acf.csv", index=False)
    ac11.to_csv(args.out / f"p11_{args.channel}_calendar_acf.csv", index=False)
    cc.to_csv(args.out / f"p10_p11_{args.channel}_calendar_ccf.csv", index=False)

    ls10 = lomb_scargle_periodogram(p10, args.channel, min_period_days=20, max_period_days=args.max_lag)
    ls11 = lomb_scargle_periodogram(p11, args.channel, min_period_days=20, max_period_days=args.max_lag)
    ls10.to_csv(args.out / f"p10_{args.channel}_lomb_scargle.csv", index=False)
    ls11.to_csv(args.out / f"p11_{args.channel}_lomb_scargle.csv", index=False)

    sw10 = sampling_window_periodogram(p10, args.channel, 20, args.max_lag)
    sw11 = sampling_window_periodogram(p11, args.channel, 20, args.max_lag)
    sw10.to_csv(args.out / f"p10_{args.channel}_sampling_window.csv", index=False)
    sw11.to_csv(args.out / f"p11_{args.channel}_sampling_window.csv", index=False)

    summary = {
        "phase": "B1 calendar-aware gap/annual-artifact falsification",
        "channel": args.channel,
        "definitions": {
            "calendar_acf": "pairs observations separated by an exact number of calendar days; missing dates remain missing",
            "calendar_ccf_positive_lag": "positive L pairs P10(t+L) with P11(t), matching the dissertation-oriented lag sign",
            "cosine_raw": "non-demeaned pairwise normalized dot product; closest in spirit to the historical xcorr coefficient",
            "pearson_raw": "Pearson correlation on raw flux pairs",
            "pearson_log10": "Pearson correlation after log10 transform of strictly positive flux pairs",
        },
        "historical_sample_lags_as_calendar_time": {
            "p10_acf_469_samples": calendar_separation_for_single_sample_lag(p10, args.channel, 469),
            "p11_acf_467_samples": calendar_separation_for_single_sample_lag(p11, args.channel, 467),
            "p11_acf_1697_samples": calendar_separation_for_single_sample_lag(p11, args.channel, 1697),
        },
        "calendar_acf": {
            "p10": {
                "annual_band_330_400": strongest_calendar_lag(ac10, 330, 400),
                "one_point_three_year_band_430_520": strongest_calendar_lag(ac10, 430, 520),
                "exact_340": calendar_lag_exact(ac10, 340),
                "exact_365": calendar_lag_exact(ac10, 365),
                "exact_469": calendar_lag_exact(ac10, 469),
                "exact_470": calendar_lag_exact(ac10, 470),
            },
            "p11": {
                "annual_band_330_400": strongest_calendar_lag(ac11, 330, 400),
                "one_point_three_year_band_430_520": strongest_calendar_lag(ac11, 430, 520),
                "exact_340": calendar_lag_exact(ac11, 340),
                "exact_365": calendar_lag_exact(ac11, 365),
                "exact_467": calendar_lag_exact(ac11, 467),
                "exact_470": calendar_lag_exact(ac11, 470),
            },
        },
        "calendar_crosscorrelation": {
            "annual_band_330_400": strongest_calendar_lag(cc, 330, 400),
            "one_point_three_year_band_430_520": strongest_calendar_lag(cc, 430, 520),
            "exact_340": calendar_lag_exact(cc, 340),
            "exact_341": calendar_lag_exact(cc, 341),
            "exact_365": calendar_lag_exact(cc, 365),
            "exact_470": calendar_lag_exact(cc, 470),
        },
        "lomb_scargle_actual_observation_times": {
            "p10": {
                "annual_band_330_400": band_peak(ls10, 330, 400),
                "one_point_three_year_band_430_520": band_peak(ls10, 430, 520),
            },
            "p11": {
                "annual_band_330_400": band_peak(ls11, 330, 400),
                "one_point_three_year_band_430_520": band_peak(ls11, 430, 520),
            },
        },
        "sampling_window_spectrum": {
            "p10": {
                "annual_band_330_400": strongest_period_in_band(sw10, 330, 400),
                "one_point_three_year_band_430_520": strongest_period_in_band(sw10, 430, 520),
            },
            "p11": {
                "annual_band_330_400": strongest_period_in_band(sw11, 330, 400),
                "one_point_three_year_band_430_520": strongest_period_in_band(sw11, 430, 520),
            },
        },
        "interpretation_guardrail": (
            "Phase A established exact computational reproducibility of the 1998 sample-index results. "
            "Phase B1 must not call a sample-index lag a physical period unless the feature survives exact-calendar and irregular-time analysis."
        ),
    }

    out_json = args.out / "calendar_reanalysis_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
