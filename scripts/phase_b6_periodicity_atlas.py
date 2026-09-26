"""Phase B6: all-channel periodicity atlas and sampling-window audit.

B5 was intentionally conditioned on B4 cross-spacecraft survivors. That is useful
for following up unusual P10/P11 lags, but it is not a valid selection rule for a
survey of periodicity inside the individual spacecraft records: a real periodic
feature may exist in P10 or P11 even when the cross-spacecraft CCF is weak.

B6 therefore removes the CCF-survivor selection and audits all 12 CPI particle
channels in the two predeclared calendar bands. For each channel, spacecraft,
robustness variant and band it reports:

* calendar-time ACF morphology;
* irregular-time Lomb-Scargle morphology;
* the sampling-window spectrum on the exact positive-flux observation mask;
* the distance between the signal spectral peak and the nearest sampling-window
  peak in Rayleigh-frequency units;
* the publication-resolution B3 ACF p/Holm/FDR values when available;
* P10/P11 spectral concordance on the exact common valid calendar span.

No new p-values are created. B6 is a morphology/alias audit used to select a
small set for later publication-grade periodogram and time-frequency inference.
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
from src.modern import dense_log_grid, fast_calendar_dcf_from_grids, fast_bin_dcf, dcf_band_max, lomb_scargle_log_periodogram

BANDS = {
    "annual_330_400": (330, 400),
    "one_point_three_year_430_520": (430, 520),
}


def _finite_span(df: pd.DataFrame, channel: str):
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


def dcf_peak_morphology(binned: pd.DataFrame, lo: int, hi: int, min_pairs: int = 100) -> dict:
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
    width = float(q.lag_bin_end.iloc[0] - q.lag_bin_start.iloc[0] + 1)
    at_edge = target <= lo + width / 2 or target >= hi - width / 2
    return {**peak, "at_band_edge": bool(at_edge), "is_local_peak": bool(len(pos)), "prominence": prom}


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


def top_spectral_peaks(periodogram: pd.DataFrame, n: int = 8) -> list[dict]:
    q = periodogram[np.isfinite(periodogram.power)].copy().sort_values("period_days").reset_index(drop=True)
    if len(q) < 3:
        return []
    y = q.power.to_numpy(float)
    peaks, props = signal.find_peaks(y, prominence=0)
    rows = []
    for j, idx in enumerate(peaks):
        rows.append({
            "period_days": float(q.loc[idx, "period_days"]),
            "power": float(q.loc[idx, "power"]),
            "prominence": float(props["prominences"][j]),
        })
    rows.sort(key=lambda r: (r["prominence"], r["power"]), reverse=True)
    return rows[: int(n)]


def sampling_window_periodogram_positive(df: pd.DataFrame, channel: str, start: pd.Timestamp, stop: pd.Timestamp,
                                         min_period: int = 200, max_period: int = 700) -> pd.DataFrame:
    """FFT spectrum of the exact positive-flux observation mask used by log10 LS."""
    idx, grid = dense_log_grid(df, channel, start=start, stop=stop, detrend="none")
    if len(grid) < 10:
        return pd.DataFrame(columns=["period_days", "power"])
    mask = np.isfinite(grid).astype(float)
    if mask.sum() < 10 or np.all(mask == mask[0]):
        return pd.DataFrame(columns=["period_days", "power"])
    y = mask - float(mask.mean())
    freq, power = signal.periodogram(y, fs=1.0, detrend=False, scaling="spectrum")
    keep = freq > 0
    freq = freq[keep]
    power = power[keep]
    period = 1.0 / freq
    keep = (period >= float(min_period)) & (period <= float(max_period))
    return pd.DataFrame({"period_days": period[keep], "power": power[keep]}).sort_values("period_days").reset_index(drop=True)


def nearest_window_peak_alias(signal_period_days: float, window_peaks: list[dict], span_days: int) -> dict:
    """Compare a signal period with the nearest local spectral-window peak.

    Rayleigh frequency resolution is 1/T.  A separation <=1 Rayleigh unit is a
    warning that the two frequencies are not cleanly resolved over the available
    span; it is not proof that the signal is an alias.
    """
    if not np.isfinite(signal_period_days) or signal_period_days <= 0 or not window_peaks or span_days <= 0:
        return {
            "nearest_window_period_days": np.nan,
            "absolute_period_difference_days": np.nan,
            "frequency_separation_rayleigh": np.nan,
            "within_one_rayleigh": None,
        }
    f = 1.0 / float(signal_period_days)
    best = min(window_peaks, key=lambda r: abs((1.0 / float(r["period_days"])) - f))
    wp = float(best["period_days"])
    ray = abs((1.0 / wp) - f) * float(span_days)
    return {
        "nearest_window_period_days": wp,
        "nearest_window_power": float(best.get("power", np.nan)),
        "nearest_window_prominence": float(best.get("prominence", np.nan)),
        "absolute_period_difference_days": abs(float(signal_period_days) - wp),
        "frequency_separation_rayleigh": float(ray),
        "within_one_rayleigh": bool(ray <= 1.0),
    }


def _acf_profile(df: pd.DataFrame, channel: str, start, stop, detrend: str,
                 min_lag: int, max_lag: int, bin_days: int) -> pd.DataFrame:
    _, grid = dense_log_grid(df, channel, start=start, stop=stop, detrend=detrend)
    if len(grid) == 0:
        return pd.DataFrame(columns=["lag_bin_start", "lag_bin_end", "lag_bin_center", "pairs", "dcf"])
    dcf = fast_calendar_dcf_from_grids(grid, None, min_lag, max_lag)
    return fast_bin_dcf(dcf, bin_days=bin_days, origin=min_lag)


def _ls_profile(df: pd.DataFrame, channel: str, start, stop, detrend: str,
                min_period: int, max_period: int, nfreq: int) -> pd.DataFrame:
    return lomb_scargle_log_periodogram(
        _subset(df, start, stop), channel,
        min_period_days=float(min_period), max_period_days=float(max_period),
        nfreq=int(nfreq), detrend=detrend,
    )


def _load_b3(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _b3_index(b3: dict) -> dict:
    out = {}
    for ch in b3.get("channels", []):
        name = ch.get("channel")
        for variant, v in ch.get("variants", {}).items():
            for sc, branch in (("p10", "p10_acf"), ("p11", "p11_acf")):
                node = v.get(branch, {})
                for band, br in node.get("bands", {}).items():
                    out[(name, variant, sc, band)] = {
                        "p_value_band_max_ar1": br.get("p_value_band_max_ar1"),
                        "p_holm_across_12_channels": br.get("p_holm_across_12_channels"),
                        "q_bh_fdr_across_12_channels": br.get("q_bh_fdr_across_12_channels"),
                        "observed_dcf": br.get("observed", {}).get("dcf"),
                        "observed_lag_calendar_days": br.get("observed", {}).get("lag_calendar_days"),
                    }
    return out


def _analyze_spacecraft(df: pd.DataFrame, channel: str, start: pd.Timestamp, stop: pd.Timestamp,
                        detrend: str, band: tuple[int, int], args) -> dict:
    lo, hi = band
    acf = _acf_profile(df, channel, start, stop, detrend, args.min_period, args.max_period, args.bin_days)
    ls = _ls_profile(df, channel, start, stop, detrend, args.min_period, args.max_period, args.nfreq)
    win = sampling_window_periodogram_positive(df, channel, start, stop, args.min_period, args.max_period)
    acf_m = dcf_peak_morphology(acf, lo, hi)
    ls_m = spectrum_peak_morphology(ls, lo, hi)
    win_m = spectrum_peak_morphology(win, lo, hi)
    wpeaks = top_spectral_peaks(win, n=12)
    span_days = int((pd.Timestamp(stop) - pd.Timestamp(start)).days + 1)
    alias = nearest_window_peak_alias(float(ls_m.get("period_days", np.nan)), wpeaks, span_days)
    return {
        "calendar": {"start": str(start.date()), "stop": str(stop.date()), "span_days": span_days},
        "acf_band_morphology": acf_m,
        "lomb_scargle_band_morphology": ls_m,
        "sampling_window_band_morphology": win_m,
        "sampling_window_top_peaks_scan": wpeaks[:8],
        "signal_vs_sampling_window": alias,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p10", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv")
    ap.add_argument("--p11", type=Path, default=ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv")
    ap.add_argument("--b3", type=Path, default=ROOT / "results" / "phase_b3_publication" / "significance_summary.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "phase_b6_periodicity_atlas")
    ap.add_argument("--min-period", type=int, default=200)
    ap.add_argument("--max-period", type=int, default=700)
    ap.add_argument("--bin-days", type=int, default=10)
    ap.add_argument("--nfreq", type=int, default=3000)
    ap.add_argument("--channels", nargs="*", default=None)
    args = ap.parse_args()

    if not args.p10.exists() or not args.p11.exists():
        raise SystemExit("Processed P10/P11 CPI files are missing. Complete Phase A first.")

    p10 = load_auto(args.p10, args.p10.name)
    p11 = apply_p11_cpi_quality_mask(load_auto(args.p11, args.p11.name))
    b3 = _load_b3(args.b3)
    b3idx = _b3_index(b3)

    requested = list(args.channels) if args.channels else list(PARTICLE_COLUMNS)
    channels = [c for c in requested if c in p10.columns and c in p11.columns]
    if not channels:
        raise SystemExit("No requested channels exist in both P10 and P11 processed files.")

    args.out.mkdir(parents=True, exist_ok=True)
    out = {
        "phase": "B6 all-channel individual-spacecraft periodicity atlas and sampling-window audit",
        "purpose": (
            "Audit periodicity in all CPI channels without conditioning on cross-spacecraft CCF survival, "
            "and flag spectral peaks that are not cleanly separated from the observation-window spectrum."
        ),
        "guardrails": [
            "B6 creates no new p-values; B3 publication-resolution ACF significance is attached when available.",
            "A Lomb-Scargle local peak is morphology, not significance.",
            "A signal peak within one Rayleigh frequency resolution of a sampling-window peak is flagged as potentially unresolved from the observation window; this is a warning, not proof of aliasing.",
            "P10/P11 concordance is descriptive and uses the exact common valid calendar span for each channel.",
            "Channels affected by the Pioneer 11 post-1980 CPI reliability mask have a shorter common span.",
        ],
        "period_scan_calendar_days": [int(args.min_period), int(args.max_period)],
        "bands": BANDS,
        "b3_publication_summary_found": bool(b3),
        "channels": [],
    }
    flat = []

    for ci, channel in enumerate(channels, start=1):
        print(f"[{ci}/{len(channels)}] B6 {channel}", flush=True)
        s10, e10 = _finite_span(p10, channel)
        s11, e11 = _finite_span(p11, channel)
        if s10 is None or s11 is None:
            out["channels"].append({"channel": channel, "status": "insufficient_data"})
            continue
        common_start, common_stop = max(s10, s11), min(e10, e11)
        rec = {
            "channel": channel,
            "label": CHANNEL_LABELS.get(channel, channel),
            "p11_quality_limited_after_1980_239": bool(e11 < pd.Timestamp("1981-01-01")),
            "common_calendar": {"start": str(common_start.date()), "stop": str(common_stop.date()),
                                "span_days": int((common_stop - common_start).days + 1)},
            "bands": {},
        }
        for band_name, band in BANDS.items():
            bnode = {"spacecraft": {}, "common_overlap_concordance": {}}
            for variant, detrend in (("raw_log10", "none"), ("linear_detrended_log10", "linear")):
                vsc = {}
                for sc_name, df, s, e in (("p10", p10, s10, e10), ("p11", p11, s11, e11)):
                    node = _analyze_spacecraft(df, channel, s, e, detrend, band, args)
                    node["b3_acf_significance"] = b3idx.get((channel, variant, sc_name, band_name))
                    vsc[sc_name] = node

                # Exact common-span LS comparison: same start/stop for both spacecraft.
                c10 = _analyze_spacecraft(p10, channel, common_start, common_stop, detrend, band, args)
                c11 = _analyze_spacecraft(p11, channel, common_start, common_stop, detrend, band, args)
                p10p = float(c10["lomb_scargle_band_morphology"].get("period_days", np.nan))
                p11p = float(c11["lomb_scargle_band_morphology"].get("period_days", np.nan))
                if np.isfinite(p10p) and np.isfinite(p11p) and max(p10p, p11p) > 0:
                    absdiff = abs(p10p - p11p)
                    reldiff = absdiff / ((p10p + p11p) / 2.0)
                else:
                    absdiff = reldiff = np.nan
                concord = {
                    "p10_period_days": p10p,
                    "p11_period_days": p11p,
                    "absolute_difference_days": float(absdiff) if np.isfinite(absdiff) else np.nan,
                    "relative_difference": float(reldiff) if np.isfinite(reldiff) else np.nan,
                    "within_10_percent": bool(reldiff <= 0.10) if np.isfinite(reldiff) else None,
                    "p10_is_local_peak": c10["lomb_scargle_band_morphology"].get("is_local_peak"),
                    "p11_is_local_peak": c11["lomb_scargle_band_morphology"].get("is_local_peak"),
                    "p10_at_band_edge": c10["lomb_scargle_band_morphology"].get("at_band_edge"),
                    "p11_at_band_edge": c11["lomb_scargle_band_morphology"].get("at_band_edge"),
                    "p10_sampling_alias_warning": c10["signal_vs_sampling_window"].get("within_one_rayleigh"),
                    "p11_sampling_alias_warning": c11["signal_vs_sampling_window"].get("within_one_rayleigh"),
                    "p10_signal_vs_sampling_window": c10["signal_vs_sampling_window"],
                    "p11_signal_vs_sampling_window": c11["signal_vs_sampling_window"],
                }
                bnode["spacecraft"][variant] = vsc
                bnode["common_overlap_concordance"][variant] = concord

                for sc_name in ("p10", "p11"):
                    n = vsc[sc_name]
                    sig = n.get("b3_acf_significance") or {}
                    flat.append({
                        "channel": channel,
                        "label": CHANNEL_LABELS.get(channel, channel),
                        "band": band_name,
                        "variant": variant,
                        "spacecraft": sc_name,
                        "acf_lag_days": n["acf_band_morphology"].get("lag_calendar_days"),
                        "acf_dcf": n["acf_band_morphology"].get("dcf"),
                        "acf_local_peak": n["acf_band_morphology"].get("is_local_peak"),
                        "acf_band_edge": n["acf_band_morphology"].get("at_band_edge"),
                        "b3_acf_p": sig.get("p_value_band_max_ar1"),
                        "b3_acf_holm": sig.get("p_holm_across_12_channels"),
                        "b3_acf_fdr": sig.get("q_bh_fdr_across_12_channels"),
                        "ls_period_days": n["lomb_scargle_band_morphology"].get("period_days"),
                        "ls_power": n["lomb_scargle_band_morphology"].get("power"),
                        "ls_local_peak": n["lomb_scargle_band_morphology"].get("is_local_peak"),
                        "ls_band_edge": n["lomb_scargle_band_morphology"].get("at_band_edge"),
                        "nearest_window_period_days": n["signal_vs_sampling_window"].get("nearest_window_period_days"),
                        "window_separation_rayleigh": n["signal_vs_sampling_window"].get("frequency_separation_rayleigh"),
                        "sampling_alias_warning": n["signal_vs_sampling_window"].get("within_one_rayleigh"),
                    })
            rec["bands"][band_name] = bnode
        out["channels"].append(rec)

    out_path = args.out / "periodicity_atlas_summary.json"
    csv_path = args.out / "periodicity_atlas_summary.csv"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    pd.DataFrame(flat).to_csv(csv_path, index=False)

    # Compact review: no score/rank, only explicit evidence flags.
    review = []
    for rec in out["channels"]:
        if rec.get("status"):
            continue
        for band_name, bnode in rec["bands"].items():
            row = {"channel": rec["channel"], "label": rec["label"], "band": band_name,
                   "p11_quality_limited_after_1980_239": rec["p11_quality_limited_after_1980_239"]}
            for variant in ("raw_log10", "linear_detrended_log10"):
                co = bnode["common_overlap_concordance"][variant]
                row[f"{variant}__common_periods_within_10_percent"] = co["within_10_percent"]
                row[f"{variant}__both_common_ls_local_peaks"] = bool(co["p10_is_local_peak"] and co["p11_is_local_peak"])
                row[f"{variant}__either_common_sampling_alias_warning"] = bool(co["p10_sampling_alias_warning"] or co["p11_sampling_alias_warning"])
                for sc_name in ("p10", "p11"):
                    n = bnode["spacecraft"][variant][sc_name]
                    sig = n.get("b3_acf_significance") or {}
                    row[f"{variant}__{sc_name}_b3_acf_holm"] = sig.get("p_holm_across_12_channels")
                    row[f"{variant}__{sc_name}_ls_period_days"] = n["lomb_scargle_band_morphology"].get("period_days")
                    row[f"{variant}__{sc_name}_ls_local_peak"] = n["lomb_scargle_band_morphology"].get("is_local_peak")
                    row[f"{variant}__{sc_name}_sampling_alias_warning"] = n["signal_vs_sampling_window"].get("within_one_rayleigh")
            review.append(row)
    review_path = args.out / "candidate_flags.json"
    review_path.write_text(json.dumps({
        "phase": "B6 candidate flags (descriptive, no ranking)",
        "interpretation": [
            "These are evidence flags, not a score or a significance claim for Lomb-Scargle peaks.",
            "B3 Holm values refer to calendar-ACF band-max red-noise tests; B6 adds independent spectral morphology and sampling-window diagnostics.",
            "A publication-facing periodicity candidate should ideally have B3-supported ACF structure, non-edge local LS peaks in both spacecraft over the same valid span, concordant periods, and no one-Rayleigh sampling-window warning before dedicated LS/time-frequency significance testing.",
        ],
        "candidates": review,
    }, indent=2), encoding="utf-8")

    print(f"\nPhase B6 complete: {out_path}")
    print(f"Flat table: {csv_path}")
    print(f"Candidate flags: {review_path}")


if __name__ == "__main__":
    main()
