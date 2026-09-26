"""Headless Phase-A replication runner.

Example:
    python scripts/replicate.py --p10 data/p10cpi.dat --p11 data/p11cpi.dat \
        --channel p_11_20_mev --out results/phase_a

The command writes machine-readable CSV/JSON outputs and SHA-256 hashes of the
inputs so the exact replication run can be frozen for a paper.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.io import load_auto, apply_p11_cpi_quality_mask
from src.replication import (
    thesis_fft,
    autocorrelation,
    crosscorrelation_common_dates,
    propagation_shifted_crosscorrelation,
    strongest_in_band,
    coefficient_at_lag,
    top_correlation_peaks,
    alignment_diagnostics,
    thesis_crosscorr_lag_from_scipy,
    calendar_separation_for_sample_lag,
    benchmark_match,
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path)
    ap.add_argument("--p11", type=Path)
    ap.add_argument("--channel", default="p_11_20_mev")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_a")
    ap.add_argument("--no-p11-mask", action="store_true")
    args = ap.parse_args()

    if args.p10 is None and args.p11 is None:
        ap.error("Provide --p10 and/or --p11")
    args.out.mkdir(parents=True, exist_ok=True)

    frames = {}
    provenance = {}
    if args.p10:
        frames["p10"] = load_auto(args.p10, args.p10.name)
        provenance["p10"] = {"path": str(args.p10), "sha256": sha256(args.p10)}
    if args.p11:
        p11 = load_auto(args.p11, args.p11.name)
        frames["p11"] = p11 if args.no_p11_mask else apply_p11_cpi_quality_mask(p11)
        provenance["p11"] = {"path": str(args.p11), "sha256": sha256(args.p11), "quality_mask": not args.no_p11_mask}

    summary = {"channel": args.channel, "inputs": provenance, "spacecraft": {}}
    for key, df in frames.items():
        if args.channel not in df.columns:
            raise SystemExit(f"{args.channel} not found in {key}")
        s = df[args.channel].dropna().to_numpy(float)
        fft = thesis_fft(s)
        acf = autocorrelation(s)
        fft.to_csv(args.out / f"{key}_{args.channel}_fft.csv", index=False)
        acf.to_csv(args.out / f"{key}_{args.channel}_acf.csv", index=False)
        sc_summary = {
            "rows": int(len(df)),
            "finite_channel_values": int(len(s)),
            "acf_target_430_520_days": strongest_in_band(acf, 430, 520),
        }
        if args.channel == "p_11_20_mev":
            if key == "p10":
                obs = coefficient_at_lag(acf, 469)
                sc_summary["thesis_benchmark"] = {
                    "expected": [{"lag_days": 469, "coefficient": 0.19}],
                    "observed_exact_lags": [obs],
                    "replication_checks": [benchmark_match(obs, 469, 0.19)],
                }
            elif key == "p11":
                obs1697 = coefficient_at_lag(acf, 1697)
                obs467 = coefficient_at_lag(acf, 467)
                sc_summary["thesis_benchmark"] = {
                    "expected": [
                        {"lag_days": 1697, "coefficient": 0.21},
                        {"lag_days": 467, "coefficient": 0.11},
                    ],
                    "observed_exact_lags": [obs1697, obs467],
                    "replication_checks": [
                        benchmark_match(obs1697, 1697, 0.21),
                        benchmark_match(obs467, 467, 0.11),
                    ],
                }
        summary["spacecraft"][key] = sc_summary

    if "p10" in frames and "p11" in frames:
        cc = crosscorrelation_common_dates(frames["p10"], frames["p11"], args.channel)
        cc.to_csv(args.out / f"p10_p11_{args.channel}_ccf.csv", index=False)
        raw_peaks = top_correlation_peaks(cc, n=16)
        thesis_peaks = [
            {
                "thesis_lag_samples": thesis_crosscorr_lag_from_scipy(p["lag_days"]),
                "scipy_lag_samples": p["lag_days"],
                "coefficient": p["coefficient"],
            }
            for p in raw_peaks
        ]
        cc_summary = {
            "lag_convention": {
                "raw": "SciPy signal.correlate lag sign",
                "thesis": "thesis_lag = -scipy_lag",
                "reason": "The 1998 P10/P11 benchmark coefficients reproduce at the opposite SciPy lag signs; amplitudes and absolute lags match.",
            },
            "target_315_365_samples_scipy_positive": strongest_in_band(cc, 315, 365),
            "target_315_365_samples_scipy_negative": strongest_in_band(cc, -365, -315),
            "top_positive_correlation_peaks_raw_scipy": raw_peaks,
            "top_positive_correlation_peaks_thesis_convention": thesis_peaks,
            "alignment": alignment_diagnostics(frames["p10"], frames["p11"], args.channel),
        }
        if args.channel == "p_11_20_mev":
            thesis_cc = [(1905, 0.88), (341, 0.24), (1524, 0.13), (-214, 0.08)]
            observed_thesis = []
            checks = []
            calendar = []
            for lag, coeff in thesis_cc:
                raw = coefficient_at_lag(cc, -lag)
                translated = {
                    "lag_days": thesis_crosscorr_lag_from_scipy(raw["lag_days"]),
                    "coefficient": raw["coefficient"],
                    "raw_scipy_lag_samples": raw["lag_days"],
                }
                observed_thesis.append(translated)
                checks.append(benchmark_match(translated, lag, coeff))
                calendar.append(calendar_separation_for_sample_lag(frames["p10"], frames["p11"], args.channel, lag))
            cc_summary["thesis_benchmark"] = {
                "expected": [{"lag_days": lag, "coefficient": coeff} for lag, coeff in thesis_cc],
                "observed_thesis_convention": observed_thesis,
                "replication_checks": checks,
                "all_four_replicated": all(c["replicated"] for c in checks),
                "calendar_separation_diagnostics": calendar,
                "interpretation_note": (
                    "The historical xcorr lag is an index/sample lag after common-date filtering. "
                    "Because the aligned series contains gaps, the thesis label 'days' is not always equal to actual calendar days. "
                    "The calendar diagnostics quantify that distinction without altering the historical replication."
                ),
            }
        summary["crosscorrelation"] = cc_summary
        required = {"x_au", "y_au", "z_au", "solar_wind_speed_km_s"}
        have_required = required.issubset(frames["p10"].columns) and required.issubset(frames["p11"].columns)
        if have_required:
            # Do not treat all-NaN placeholder columns as real geometry/plasma data.
            finite_required = all(
                frames[k][c].notna().any()
                for k in ("p10", "p11")
                for c in required
            )
        else:
            finite_required = False
        if finite_required:
            shifted_cc, shifted = propagation_shifted_crosscorrelation(frames["p10"], frames["p11"], args.channel)
            shifted_cc.to_csv(args.out / f"p10_p11_{args.channel}_ccf_shifted.csv", index=False)
            shifted.to_csv(args.out / f"p10_p11_{args.channel}_shifted_series.csv", index=False)
            summary["propagation_shifted_crosscorrelation"] = {
                "target_315_365_days": strongest_in_band(shifted_cc, 315, 365),
                "aligned_rows": int(len(shifted)),
                "median_delay_days": float(shifted.delay_days.median()),
            }
        else:
            summary["propagation_shifted_crosscorrelation"] = {
                "status": "not_run",
                "reason": "XYZ coordinates and/or solar-wind speed are not yet available from the reconstructed particle backbone.",
            }

    (args.out / "replication_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
