"""Phase B9: candidate-focused mask-preserving red-noise spectral significance.

B8 showed that five-year sampling-window warnings were often resolution-limited
and that longer windows separate some candidate peaks from the observation
window. B9 now asks a different question: does maximum spectral power inside a
locked band exceed what is expected from daily AR(1) red noise observed through
the exact historical mask?

This is post-selection evidence, not an independent replication. The locked
candidate set is written into the output before the tests are run.
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
from src.schema import CHANNEL_LABELS
from src.modern import ar1_lomb_band_max_test

LOCKED_CANDIDATES = [
    {"channel": "he_20_24_mev_n", "band": "one_point_three_year_430_520", "range": (430, 520),
     "role": "primary_after_B8"},
    {"channel": "he_11_20_mev_n", "band": "one_point_three_year_430_520", "range": (430, 520),
     "role": "secondary_long_baseline_comparator"},
]


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
        running = max(running, min(1.0, (m - i) * float(val)))
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
        running = min(running, m * float(sv[i]) / (i + 1))
        q[i] = min(1.0, running)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    adj = q[inv]
    out[np.where(good)[0]] = adj
    return out.tolist()


def _attach_adjustments(records: list[dict]) -> None:
    # One family per robustness variant: 2 channels x 2 spacecraft = 4 tests.
    for variant in ("raw_log10", "linear_detrended_log10"):
        refs = []
        ps = []
        for rec in records:
            for sc in ("p10", "p11"):
                node = rec["variants"][variant][sc]
                refs.append(node)
                ps.append(float(node.get("p_value_band_max_ar1_lomb", np.nan)))
        hs = _holm(ps)
        qs = _bh(ps)
        for ref, h, q in zip(refs, hs, qs):
            ref["p_holm_across_locked_4_tests"] = h
            ref["q_bh_fdr_across_locked_4_tests"] = q


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b9_spectral_significance")
    ap.add_argument("--simulations", type=int, default=1999)
    ap.add_argument("--nfreq", type=int, default=96)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--seed", type=int, default=1998)
    ap.add_argument("--no-p11-mask", action="store_true")
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A acquisition first.")
    args.out.mkdir(parents=True, exist_ok=True)
    p10 = load_auto(args.p10, args.p10.name)
    p11raw = load_auto(args.p11, args.p11.name)
    p11 = p11raw if args.no_p11_mask else apply_p11_cpi_quality_mask(p11raw)

    out = {
        "phase": "B9 candidate-focused mask-preserving AR1 spectral significance",
        "purpose": (
            "Test whether maximum spectral power inside the B8-locked ~1.3-year candidate band exceeds "
            "daily AR(1) red-noise expectations when the exact historical observation mask is reapplied."
        ),
        "guardrails": [
            "B9 candidates were selected after examining earlier phases; these p-values are post-selection evidence, not an independent replication.",
            "The statistic is the maximum classical Lomb-Scargle power anywhere inside the locked 430-520 d band, not a single chosen frequency.",
            "Daily AR(1) surrogates are generated on the complete calendar and sampled through the exact positive-flux observation mask.",
            "Holm and Benjamini-Hochberg adjustments are applied across the locked family of four spacecraft/channel tests separately for each robustness variant.",
            "A significant B9 result still requires time-frequency localization before a physical quasi-periodicity claim."
        ],
        "locked_candidates": LOCKED_CANDIDATES,
        "simulations_per_test": int(args.simulations),
        "frequency_grid_points_per_band": int(args.nfreq),
        "records": [],
    }

    variants = [("raw_log10", "none"), ("linear_detrended_log10", "linear")]
    for ci, cand in enumerate(LOCKED_CANDIDATES, start=1):
        channel = cand["channel"]
        lo, hi = cand["range"]
        print(f"[{ci}/{len(LOCKED_CANDIDATES)}] B9 {channel} {lo}-{hi} d", flush=True)
        rec = {
            "channel": channel,
            "label": CHANNEL_LABELS.get(channel, channel),
            "band": cand["band"],
            "role": cand["role"],
            "variants": {},
        }
        for vi, (vname, detrend) in enumerate(variants):
            rec["variants"][vname] = {}
            for si, (sc, data) in enumerate((("p10", p10), ("p11", p11))):
                print(f"  - {vname} {sc}", flush=True)
                seed = int(args.seed + ci * 10000 + vi * 1000 + si * 100)
                rec["variants"][vname][sc] = ar1_lomb_band_max_test(
                    data, channel, lo, hi, simulations=args.simulations,
                    detrend=detrend, nfreq=args.nfreq, batch_size=args.batch_size,
                    seed=seed,
                )
        out["records"].append(rec)

    _attach_adjustments(out["records"])

    rows = []
    for rec in out["records"]:
        for variant, v in rec["variants"].items():
            for sc, node in v.items():
                rows.append({
                    "channel": rec["channel"], "label": rec["label"], "role": rec["role"],
                    "variant": variant, "spacecraft": sc,
                    "observed_period_days": node.get("observed", {}).get("period_days"),
                    "observed_power": node.get("observed", {}).get("power"),
                    "p_ar1_band_max_lomb": node.get("p_value_band_max_ar1_lomb"),
                    "holm_locked4": node.get("p_holm_across_locked_4_tests"),
                    "bh_locked4": node.get("q_bh_fdr_across_locked_4_tests"),
                    "phi_1day": node.get("ar1", {}).get("phi_1day"),
                    "observations": node.get("observations"),
                    "simulations": node.get("null_simulations_used"),
                })
    flat = pd.DataFrame(rows)
    out_path = args.out / "spectral_significance_summary.json"
    csv_path = args.out / "spectral_significance_summary.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    flat.to_csv(csv_path, index=False)
    print("\n" + flat.to_string(index=False))
    print(f"\nPhase B9 complete: {out_path}")


if __name__ == "__main__":
    main()
