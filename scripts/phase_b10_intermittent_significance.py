"""Phase B10: intermittent same-window/same-period red-noise significance.

B9 found no full-record global ~1.3-year spectral significance for the two locked
long-baseline helium candidates. B10 tests the distinct, predeclared alternative
that the signal is intermittent: the maximum over fixed 8-year windows and the
locked 430-520 d band is calibrated against mask-preserving daily AR(1) surrogates.
A joint P10/P11 statistic requires power at the same period in the same window.
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
from src.modern import ar1_lomb_intermit_joint_test

LOCKED_CANDIDATES = [
    {"channel": "he_20_24_mev_n", "band": "one_point_three_year_430_520", "range": (430, 520),
     "role": "primary_after_B8"},
    {"channel": "he_11_20_mev_n", "band": "one_point_three_year_430_520", "range": (430, 520),
     "role": "secondary_long_baseline_comparator"},
]


def _holm(p: list[float]) -> list[float]:
    arr = np.asarray(p, dtype=float); out = np.full(len(arr), np.nan)
    good = np.isfinite(arr); vals = arr[good]
    if len(vals) == 0: return out.tolist()
    order = np.argsort(vals); sv = vals[order]; m = len(sv)
    adj_sorted = np.empty(m); running = 0.0
    for i, val in enumerate(sv):
        running = max(running, min(1.0, (m - i) * float(val))); adj_sorted[i] = running
    inv = np.empty(m, dtype=int); inv[order] = np.arange(m)
    out[np.where(good)[0]] = adj_sorted[inv]
    return out.tolist()


def _bh(p: list[float]) -> list[float]:
    arr = np.asarray(p, dtype=float); out = np.full(len(arr), np.nan)
    good = np.isfinite(arr); vals = arr[good]
    if len(vals) == 0: return out.tolist()
    order = np.argsort(vals); sv = vals[order]; m = len(sv); q = np.empty(m); running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, m * float(sv[i]) / (i + 1)); q[i] = min(1.0, running)
    inv = np.empty(m, dtype=int); inv[order] = np.arange(m)
    out[np.where(good)[0]] = q[inv]
    return out.tolist()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b10_intermittent_significance")
    ap.add_argument("--simulations", type=int, default=1999)
    ap.add_argument("--nfreq", type=int, default=96)
    ap.add_argument("--window-years", type=int, default=8)
    ap.add_argument("--step-years", type=int, default=2)
    ap.add_argument("--batch-size", type=int, default=16)
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
        "phase": "B10 intermittent same-window/same-period AR1 spectral significance",
        "purpose": (
            "Test the distinct hypothesis left open by B9: an intermittent ~1.3-year signal may be weak globally but recur "
            "during limited epochs in both spacecraft. The look-elsewhere scan over fixed windows and frequencies is included in the null statistic."
        ),
        "guardrails": [
            "B9 found no full-record global spectral significance for either locked candidate; B10 therefore tests intermittency, not persistent periodicity.",
            "Candidates remain post-selection from earlier phases; B10 is not independent replication.",
            "Windows are fixed at 8 calendar years with a 2-year step before B10 results are examined.",
            "The joint statistic requires P10 and P11 Lomb power at the same period in the same window and is maximized over the entire locked search.",
            "Daily AR(1) surrogates are generated independently for P10/P11 on the common complete calendar and sampled through the exact positive-flux masks.",
            "No interpolation is used. A non-significant joint result is evidence against the remaining intermittent shared-signal hypothesis under this null model."
        ],
        "locked_candidates": LOCKED_CANDIDATES,
        "window_years": int(args.window_years), "step_years": int(args.step_years),
        "simulations_per_test": int(args.simulations), "frequency_grid_points_per_band": int(args.nfreq),
        "records": [],
    }

    variants = [("raw_log10", "none"), ("linear_detrended_log10", "linear")]
    for ci, cand in enumerate(LOCKED_CANDIDATES, start=1):
        channel = cand["channel"]; lo, hi = cand["range"]
        print(f"[{ci}/{len(LOCKED_CANDIDATES)}] B10 {channel} {lo}-{hi} d", flush=True)
        rec = {"channel": channel, "label": CHANNEL_LABELS.get(channel, channel),
               "band": cand["band"], "role": cand["role"], "variants": {}}
        for vi, (vname, detrend) in enumerate(variants):
            seed = int(args.seed + ci * 10000 + vi * 1000)
            rec["variants"][vname] = ar1_lomb_intermit_joint_test(
                p10, p11, channel, lo, hi, simulations=args.simulations, detrend=detrend,
                nfreq=args.nfreq, window_years=args.window_years, step_years=args.step_years,
                batch_size=args.batch_size, seed=seed,
            )
        out["records"].append(rec)

    # Primary inferential family: one joint P10/P11 test per locked candidate, separately by robustness variant.
    for variant in ("raw_log10", "linear_detrended_log10"):
        nodes = [rec["variants"][variant]["null"]["joint"] for rec in out["records"]]
        ps = [float(n.get("p_value", np.nan)) for n in nodes]
        hs = _holm(ps); qs = _bh(ps)
        for n, h, q in zip(nodes, hs, qs):
            n["p_holm_across_locked_2_joint_tests"] = h
            n["q_bh_fdr_across_locked_2_joint_tests"] = q

    rows = []
    for rec in out["records"]:
        for variant, node in rec["variants"].items():
            obs = node.get("observed", {})
            for stat in ("p10", "p11", "joint"):
                nn = node.get("null", {}).get(stat, {})
                if stat == "p10": oo = obs.get("p10_max_over_time_frequency", {})
                elif stat == "p11": oo = obs.get("p11_max_over_time_frequency", {})
                else: oo = obs.get("joint_same_window_same_period_max", {})
                rows.append({
                    "channel": rec["channel"], "label": rec["label"], "role": rec["role"],
                    "variant": variant, "statistic": stat,
                    "window_start": oo.get("window_start"), "window_stop": oo.get("window_stop"),
                    "period_days": oo.get("period_days"), "observed_power": oo.get("power"),
                    "p_ar1_time_frequency_max": nn.get("p_value"),
                    "holm_locked2_joint": nn.get("p_holm_across_locked_2_joint_tests") if stat == "joint" else np.nan,
                    "bh_locked2_joint": nn.get("q_bh_fdr_across_locked_2_joint_tests") if stat == "joint" else np.nan,
                    "null_q95": nn.get("null_q95"), "null_q99": nn.get("null_q99"),
                })
    flat = pd.DataFrame(rows)
    out_path = args.out / "intermittent_significance_summary.json"
    csv_path = args.out / "intermittent_significance_summary.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    flat.to_csv(csv_path, index=False)
    print("\n" + flat.to_string(index=False))
    print(f"\nPhase B10 complete: {out_path}")


if __name__ == "__main__":
    main()
