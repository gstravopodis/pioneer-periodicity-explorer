"""Phase B3: calendar-aware red-noise / surrogate significance testing.

This phase tests the fixed annual (330-400 d) and ~1.3-year (430-520 d)
bands across all recovered CPI channels.  It deliberately keeps Phase A frozen
and does not reinterpret historical sample-index lags as calendar days.

Two null models are used:
  * ACF: daily AR(1) red noise simulated on the full calendar, then sampled with
    the exact historical missing-data mask.  The test statistic is the maximum
    10-day-binned DCF value anywhere inside each predeclared band.
  * P10/P11 CCF: random circular calendar shifts of the entire P11 value+mask
    series inside the common observing span.  This preserves P11's own red-noise
    structure and sampling pattern while breaking absolute P10/P11 alignment.

Both raw log10 flux and a linearly detrended log10 robustness variant are run.
P-values are band-max corrected internally and then adjusted across channels
with Holm (strict family-wise control) and Benjamini-Hochberg FDR.
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
from src.modern import ar1_dcf_max_surrogate_test, circular_shift_ccf_max_test

BANDS = {
    "annual_330_400": (330, 400),
    "one_point_three_year_430_520": (430, 520),
}


def _holm(p: list[float]) -> list[float]:
    arr = np.asarray(p, dtype=float)
    out = np.full(len(arr), np.nan)
    good = np.isfinite(arr)
    vals = arr[good]
    if len(vals) == 0:
        return out.tolist()
    order = np.argsort(vals)
    sv = vals[order]
    m = len(sv)
    adj_sorted = np.empty(m, dtype=float)
    running = 0.0
    for i, val in enumerate(sv):
        a = min(1.0, (m - i) * float(val))
        running = max(running, a)
        adj_sorted[i] = running
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    adj = adj_sorted[inv]
    out[np.where(good)[0]] = adj
    return out.tolist()


def _bh(p: list[float]) -> list[float]:
    arr = np.asarray(p, dtype=float)
    out = np.full(len(arr), np.nan)
    good = np.isfinite(arr)
    vals = arr[good]
    if len(vals) == 0:
        return out.tolist()
    order = np.argsort(vals)
    sv = vals[order]
    m = len(sv)
    q = np.empty(m, dtype=float)
    running = 1.0
    for i in range(m - 1, -1, -1):
        rank = i + 1
        running = min(running, m * float(sv[i]) / rank)
        q[i] = min(1.0, running)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    adj = q[inv]
    out[np.where(good)[0]] = adj
    return out.tolist()


def _attach_multiple_testing(channels: list[dict]) -> None:
    for variant in ("raw_log10", "linear_detrended_log10"):
        for band in BANDS:
            for branch, pkey in (
                ("p10_acf", "p_value_band_max_ar1"),
                ("p11_acf", "p_value_band_max_ar1"),
                ("p10_p11_ccf", "p_value_band_max_circular_shift"),
            ):
                ps = []
                refs = []
                for ch in channels:
                    node = ch["variants"][variant][branch]
                    if node.get("status") != "ok":
                        ps.append(np.nan)
                        refs.append(None)
                        continue
                    b = node["bands"][band]
                    ps.append(float(b.get(pkey, np.nan)))
                    refs.append(b)
                holm = _holm(ps)
                bh = _bh(ps)
                for ref, h, q in zip(refs, holm, bh):
                    if ref is not None:
                        ref["p_holm_across_12_channels"] = h
                        ref["q_bh_fdr_across_12_channels"] = q


def _flatten(out: dict) -> pd.DataFrame:
    rows = []
    for ch in out["channels"]:
        for variant, v in ch["variants"].items():
            for band in BANDS:
                a10 = v["p10_acf"].get("bands", {}).get(band, {})
                a11 = v["p11_acf"].get("bands", {}).get(band, {})
                cc = v["p10_p11_ccf"].get("bands", {}).get(band, {})
                rows.append({
                    "channel": ch["channel"],
                    "label": ch["label"],
                    "variant": variant,
                    "band": band,
                    "p10_observed_lag": a10.get("observed", {}).get("lag_calendar_days"),
                    "p10_observed_dcf": a10.get("observed", {}).get("dcf"),
                    "p10_p_ar1": a10.get("p_value_band_max_ar1"),
                    "p10_holm": a10.get("p_holm_across_12_channels"),
                    "p10_bh": a10.get("q_bh_fdr_across_12_channels"),
                    "p11_observed_lag": a11.get("observed", {}).get("lag_calendar_days"),
                    "p11_observed_dcf": a11.get("observed", {}).get("dcf"),
                    "p11_p_ar1": a11.get("p_value_band_max_ar1"),
                    "p11_holm": a11.get("p_holm_across_12_channels"),
                    "p11_bh": a11.get("q_bh_fdr_across_12_channels"),
                    "ccf_observed_lag": cc.get("observed", {}).get("lag_calendar_days"),
                    "ccf_observed_dcf": cc.get("observed", {}).get("dcf"),
                    "ccf_p_shift": cc.get("p_value_band_max_circular_shift"),
                    "ccf_holm": cc.get("p_holm_across_12_channels"),
                    "ccf_bh": cc.get("q_bh_fdr_across_12_channels"),
                })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--simulations", type=int, default=199,
                    help="Surrogates per channel/variant. 199 is a fast first pass; use 1999+ for final inference.")
    ap.add_argument("--bin-days", type=int, default=10)
    ap.add_argument("--min-shift-days", type=int, default=700)
    ap.add_argument("--seed", type=int, default=1998)
    ap.add_argument("--no-p11-mask", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b3_significance")
    ap.add_argument("--channels", nargs="*", default=None,
                    help="Optional subset. Default: all 12 CPI particle channels.")
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A acquisition first.")
    args.out.mkdir(parents=True, exist_ok=True)
    p10 = load_auto(args.p10, args.p10.name)
    p11_raw = load_auto(args.p11, args.p11.name)
    p11 = p11_raw if args.no_p11_mask else apply_p11_cpi_quality_mask(p11_raw)

    requested = list(args.channels) if args.channels else list(PARTICLE_COLUMNS)
    channels = [c for c in requested if c in p10.columns and c in p11.columns]
    if not channels:
        raise SystemExit("No requested channels exist in both P10 and P11 processed files.")

    out = {
        "phase": "B3 red-noise and alignment-surrogate significance",
        "purpose": (
            "Test the two fixed calendar bands after B2 without interpolation, while preserving actual gaps. "
            "ACF uses a daily AR(1) null sampled through the historical mask; cross-spacecraft CCF uses "
            "random circular calendar shifts of P11 to preserve its own red-noise and sampling structure."
        ),
        "guardrails": [
            "P-values are based on the maximum DCF bin anywhere inside each predeclared band, not the single observed peak bin.",
            "Holm and Benjamini-Hochberg adjustments are applied across channels separately for each band, statistic and robustness variant.",
            "A small CCF shift-null p-value means absolute P10/P11 alignment is unusual under this surrogate scheme; it does not by itself identify a physical propagation mechanism.",
            "The B3 default of 199 simulations is a computational first pass. Any candidate intended for publication should be rerun with >=1999 simulations and independently time-localized.",
        ],
        "bands_calendar_days": BANDS,
        "simulations_per_test": int(args.simulations),
        "bin_days": int(args.bin_days),
        "min_abs_circular_shift_days": int(args.min_shift_days),
        "channels": [],
    }

    variants = [("raw_log10", "none"), ("linear_detrended_log10", "linear")]
    for ci, channel in enumerate(channels, start=1):
        print(f"[{ci}/{len(channels)}] B3 {channel}", flush=True)
        rec = {"channel": channel, "label": CHANNEL_LABELS.get(channel, channel), "variants": {}}
        for vi, (vname, detrend) in enumerate(variants):
            base_seed = int(args.seed + ci * 10000 + vi * 1000)
            print(f"  - {vname}: P10 ACF", flush=True)
            a10 = ar1_dcf_max_surrogate_test(
                p10, channel, BANDS, simulations=args.simulations, min_lag_days=250,
                max_lag_days=600, bin_days=args.bin_days, detrend=detrend, seed=base_seed + 1,
            )
            print(f"  - {vname}: P11 ACF", flush=True)
            a11 = ar1_dcf_max_surrogate_test(
                p11, channel, BANDS, simulations=args.simulations, min_lag_days=250,
                max_lag_days=600, bin_days=args.bin_days, detrend=detrend, seed=base_seed + 2,
            )
            print(f"  - {vname}: P10/P11 circular-shift CCF", flush=True)
            cc = circular_shift_ccf_max_test(
                p10, p11, channel, BANDS, simulations=args.simulations, min_lag_days=250,
                max_lag_days=600, bin_days=args.bin_days, min_shift_days=args.min_shift_days,
                detrend=detrend, seed=base_seed + 3,
            )
            rec["variants"][vname] = {"p10_acf": a10, "p11_acf": a11, "p10_p11_ccf": cc}
        out["channels"].append(rec)

    _attach_multiple_testing(out["channels"])
    flat = _flatten(out)
    out_path = args.out / "significance_summary.json"
    csv_path = args.out / "significance_summary.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    flat.to_csv(csv_path, index=False)

    for band in BANDS:
        print(f"\n{band} — raw-log10 cross-spacecraft results ranked by circular-shift p")
        q = flat[(flat.variant == "raw_log10") & (flat.band == band)].copy()
        q = q.sort_values(["ccf_p_shift", "ccf_observed_dcf"], ascending=[True, False])
        cols = ["channel", "ccf_observed_lag", "ccf_observed_dcf", "ccf_p_shift", "ccf_holm", "ccf_bh"]
        print(q[cols].to_string(index=False))

    print(f"\nPhase B3 complete: {out_path}")
    print(f"Flat table: {csv_path}")
    if args.simulations < 1999:
        print("NOTE: This is a first-pass surrogate run. Rerun surviving candidates with --simulations 1999 or higher before publication claims.")


if __name__ == "__main__":
    main()
