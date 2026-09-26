"""Post-process Phase B3 significance output into a robustness/candidate matrix.

This does not create new statistical tests. It only summarizes already-computed
B3 results and makes two interpretation constraints explicit:
  1) with N simulations, the smallest Monte-Carlo p-value is 1/(N+1);
  2) P11 CPI channels 6-13 are masked after day 239 of 1980, so those channels
     do not provide a full 1973-1992 cross-spacecraft confirmation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.schema import P11_AFFECTED_COLUMNS

BANDS = ("annual_330_400", "one_point_three_year_430_520")
VARIANTS = ("raw_log10", "linear_detrended_log10")
BRANCHES = ("p10_acf", "p11_acf", "p10_p11_ccf")


def _pkey(branch: str) -> str:
    return "p_value_band_max_circular_shift" if branch == "p10_p11_ccf" else "p_value_band_max_ar1"


def build_review(obj: dict) -> tuple[dict, pd.DataFrame]:
    sims = int(obj.get("simulations_per_test", 0) or 0)
    p_floor = 1.0 / (sims + 1) if sims >= 1 else None
    rows = []
    for ch in obj.get("channels", []):
        channel = ch.get("channel")
        label = ch.get("label", channel)
        affected = channel in set(P11_AFFECTED_COLUMNS)
        for band in BANDS:
            rec = {
                "channel": channel,
                "label": label,
                "band": band,
                "p11_quality_limited_after_1980_239": bool(affected),
            }
            for variant in VARIANTS:
                for branch in BRANCHES:
                    node = ch.get("variants", {}).get(variant, {}).get(branch, {})
                    b = node.get("bands", {}).get(band, {})
                    obs = b.get("observed", {})
                    prefix = f"{variant}__{branch}"
                    rec[f"{prefix}__lag_days"] = obs.get("lag_calendar_days")
                    rec[f"{prefix}__dcf"] = obs.get("dcf")
                    rec[f"{prefix}__p"] = b.get(_pkey(branch))
                    rec[f"{prefix}__holm"] = b.get("p_holm_across_12_channels")
                    rec[f"{prefix}__fdr"] = b.get("q_bh_fdr_across_12_channels")
            raw_fdr = rec.get("raw_log10__p10_p11_ccf__fdr")
            det_fdr = rec.get("linear_detrended_log10__p10_p11_ccf__fdr")
            raw_p = rec.get("raw_log10__p10_p11_ccf__p")
            det_p = rec.get("linear_detrended_log10__p10_p11_ccf__p")
            rec["ccf_fdr_sig_raw"] = raw_fdr is not None and raw_fdr < 0.05
            rec["ccf_fdr_sig_detrended"] = det_fdr is not None and det_fdr < 0.05
            rec["ccf_robust_both_variants_fdr05"] = rec["ccf_fdr_sig_raw"] and rec["ccf_fdr_sig_detrended"]
            rec["ccf_hits_mc_floor_raw"] = p_floor is not None and raw_p is not None and abs(raw_p - p_floor) < 1e-12
            rec["ccf_hits_mc_floor_detrended"] = p_floor is not None and det_p is not None and abs(det_p - p_floor) < 1e-12
            rows.append(rec)
    df = pd.DataFrame(rows)
    robust = df[df["ccf_robust_both_variants_fdr05"]].copy()
    result = {
        "phase": "B3 review (descriptive post-processing only)",
        "simulations_per_test": sims,
        "monte_carlo_p_floor": p_floor,
        "holm_floor_if_12_equal_min_p": (12 * p_floor) if p_floor is not None else None,
        "interpretation": [
            "This review does not add new tests or new p-values.",
            "With 199 simulations the p-value floor is 0.005, so a 12-channel Holm adjustment cannot fall below 0.06 even for a zero-exceedance test.",
            "FDR-significant first-pass findings are screening candidates until rerun at publication-resolution simulation counts.",
            "P11 columns 6-13 are quality-masked after day 239 of 1980; affected channels have a shorter P10/P11 overlap and should not be described as 20-year cross-spacecraft confirmations.",
        ],
        "robust_cross_spacecraft_fdr05_both_variants": robust[[
            "channel", "label", "band", "p11_quality_limited_after_1980_239",
            "raw_log10__p10_p11_ccf__lag_days", "raw_log10__p10_p11_ccf__dcf", "raw_log10__p10_p11_ccf__p", "raw_log10__p10_p11_ccf__fdr",
            "linear_detrended_log10__p10_p11_ccf__lag_days", "linear_detrended_log10__p10_p11_ccf__dcf", "linear_detrended_log10__p10_p11_ccf__p", "linear_detrended_log10__p10_p11_ccf__fdr",
        ]].to_dict(orient="records"),
    }
    return result, df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=ROOT / "results" / "phase_b3_significance" / "significance_summary.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b3_significance")
    args = ap.parse_args()
    obj = json.loads(args.input.read_text(encoding="utf-8"))
    review, df = build_review(obj)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "candidate_review.json").write_text(json.dumps(review, indent=2), encoding="utf-8")
    df.to_csv(args.out / "candidate_review.csv", index=False)
    print(json.dumps(review, indent=2))
    print(f"\nReview written to {args.out / 'candidate_review.json'}")


if __name__ == "__main__":
    main()
