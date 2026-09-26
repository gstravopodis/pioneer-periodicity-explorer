from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .acquisition import HapiDatasetInfo
from .schema import THESIS_COLUMNS, NASA_FIELD_ALIASES, PARTICLE_COLUMNS


def _norm(text: object) -> str:
    return re.sub(r"[^A-Z0-9]+", "", str(text).upper())


def _parameter_text(p: dict) -> str:
    fields = [p.get("name", ""), p.get("description", ""), p.get("units", ""), p.get("label", "")]
    return " ".join(str(x) for x in fields if x is not None)


def _metadata_by_name(info: HapiDatasetInfo | dict) -> dict[str, dict]:
    pars = info.parameters if isinstance(info, HapiDatasetInfo) else info.get("parameters", [])
    return {str(p.get("name", "")): p for p in pars}


def find_time_column(df: pd.DataFrame, info: HapiDatasetInfo | dict | None = None) -> str:
    meta = _metadata_by_name(info) if info is not None else {}
    for c in df.columns:
        p = meta.get(str(c), {})
        if str(p.get("type", "")).lower() == "isotime":
            return str(c)
    for c in df.columns:
        if _norm(c) in {"TIME", "EPOCH", "DATETIME", "TIMESTAMP", "UT"}:
            return str(c)
    # HAPI requires the first parameter to be time; use it only as a last resort.
    if len(df.columns):
        return str(df.columns[0])
    raise ValueError("Could not identify HAPI time column")


def _channel_score(target: str, name: str, p: dict) -> int:
    n = _norm(name)
    text = _norm(_parameter_text(p))
    score = 0
    # First prefer the field aliases used by SPDF/OMNI exports.
    if n in {_norm(k) for k, v in NASA_FIELD_ALIASES.items() if v == target}:
        score += 100

    specs = {
        "p_11_20_mev": ("PROTON", "11", "20"),
        "he_11_20_mev_n": ("HE", "11", "20"),
        "p_20_24_mev": ("PROTON", "20", "24"),
        "he_20_24_mev_n": ("HE", "20", "24"),
        "p_24_29_mev": ("PROTON", "24", "29"),
        "he_24_29_mev_n": ("HE", "24", "29"),
        "p_29_67_mev": ("PROTON", "29", "67"),
        "he_29_67_mev_n": ("HE", "29", "67"),
        "e_7_17_mev": ("ELECTRON", "7", "17"),
        "minimum_ionizing_x2": ("MINIMUMION", "2", ""),
        "ions_gt67_mev_n": ("ION", "67", ""),
        "ions_zgt5_gt67_mev_n": ("Z", "5", "67"),
    }
    words = specs[target]
    if all((not w) or (w in text) for w in words):
        score += 30
    if "MEV" in text:
        score += 2
    if target.startswith("he_") and any(x in text for x in ("ALPHA", "HELIUM", "HE4")):
        score += 8
    if target.startswith("p_") and ("PROTON" in text or re.search(r"(^|[^A-Z])P([^A-Z]|$)", _parameter_text(p).upper())):
        score += 8
    if target == "ions_gt67_mev_n" and ("ZGT5" in text or "Z>5" in _parameter_text(p).upper()):
        score -= 20
    return score


def resolve_cpi_channels(df: pd.DataFrame, info: HapiDatasetInfo | dict) -> dict[str, str]:
    """Map PDS/HAPI CPI fields to thesis particle columns.

    Resolution is deliberately conservative: a field must score positively and
    ties are rejected.  We prefer a hard failure to silently assigning the wrong
    energy channel.
    """
    meta = _metadata_by_name(info)
    result: dict[str, str] = {}
    used: set[str] = set()
    for target in PARTICLE_COLUMNS:
        ranked = []
        for c in df.columns:
            if str(c) in used:
                continue
            s = _channel_score(target, str(c), meta.get(str(c), {}))
            if s > 0:
                ranked.append((s, str(c)))
        ranked.sort(reverse=True)
        if not ranked:
            continue
        if len(ranked) > 1 and ranked[0][0] == ranked[1][0]:
            raise ValueError(f"Ambiguous HAPI CPI mapping for {target}: {ranked[:3]}")
        result[target] = ranked[0][1]
        used.add(ranked[0][1])
    return result


def _coord_score(axis: str, name: str, p: dict) -> int:
    n = _norm(name)
    text = _norm(_parameter_text(p))
    units = _norm(p.get("units", ""))
    score = 0
    aliases = {"X": "x", "Y": "y", "Z": "z"}
    if n == axis.upper() or n in {f"{axis.upper()}AU", f"HCI{axis.upper()}", f"HC{axis.upper()}"}:
        score += 50
    if "AU" in units or "AU" in text or "ASTRONOMICALUNIT" in text:
        score += 10
    if any(k in text for k in ("HELIOCENTRIC", "ECLIPTIC", "CARTESIAN", "SPACECRAFTPOSITION")):
        score += 8
    # Make the requested axis explicit and reject obvious other-axis fields.
    if axis.upper() in n:
        score += 5
    for other in set(aliases) - {axis.upper()}:
        if n == other or n.startswith(other):
            score -= 20
    return score


def resolve_coordinates(df: pd.DataFrame, info: HapiDatasetInfo | dict) -> dict[str, str]:
    meta = _metadata_by_name(info)
    out = {}
    used: set[str] = set()
    for axis, target in [("X", "x_au"), ("Y", "y_au"), ("Z", "z_au")]:
        ranked = []
        for c in df.columns:
            if str(c) in used:
                continue
            s = _coord_score(axis, str(c), meta.get(str(c), {}))
            if s > 0:
                ranked.append((s, str(c)))
        ranked.sort(reverse=True)
        if ranked:
            if len(ranked) > 1 and ranked[0][0] == ranked[1][0]:
                raise ValueError(f"Ambiguous coordinate mapping for {axis}: {ranked[:3]}")
            out[target] = ranked[0][1]
            used.add(ranked[0][1])
    return out


def resolve_bulk_speed(df: pd.DataFrame, info: HapiDatasetInfo | dict) -> str:
    meta = _metadata_by_name(info)
    ranked = []
    for c in df.columns:
        p = meta.get(str(c), {})
        n = _norm(c)
        text = _norm(_parameter_text(p))
        units = _norm(p.get("units", ""))
        score = 0
        if n in {"SW_SPEED", "SWSPEED", "VSW", "SPEED", "BULKVELOCITY", "BULKSPEED"}:
            score += 100
        if "KM/S" in str(p.get("units", "")).upper() or "KMS" in units or "KMPERSEC" in units:
            score += 15
        if "BULKVELOCITY" in text or "SOLARWINDSPEED" in text or "BULKSPEED" in text:
            score += 35
        if any(x in text for x in ("DENSITY", "TEMPERATURE", "ANGLE")):
            score -= 50
        size = p.get("size")
        if isinstance(size, list) and np.prod(size) > 1:
            score -= 30
        if score > 0:
            ranked.append((score, str(c)))
    ranked.sort(reverse=True)
    if not ranked:
        raise ValueError("Could not identify a scalar solar-wind bulk-speed field")
    if len(ranked) > 1 and ranked[0][0] == ranked[1][0]:
        raise ValueError(f"Ambiguous solar-wind speed mapping: {ranked[:3]}")
    return ranked[0][1]


def thesis_decimal_year(dates: Iterable[pd.Timestamp]) -> np.ndarray:
    d = pd.to_datetime(pd.Series(dates), errors="coerce")
    # Historical encoding uses YY + DOY/365 even in leap years.
    return (d.dt.year % 100).to_numpy(float) + d.dt.dayofyear.to_numpy(float) / 365.0


def canonicalise_cpi_daily(df: pd.DataFrame, info: HapiDatasetInfo | dict) -> tuple[pd.DataFrame, dict]:
    time_col = find_time_column(df, info)
    ch = resolve_cpi_channels(df, info)
    coords = resolve_coordinates(df, info)
    mapping = {**ch, **coords}
    out = pd.DataFrame({"date": pd.to_datetime(df[time_col], errors="coerce", utc=True).dt.tz_localize(None).dt.normalize()})
    for target, source in mapping.items():
        out[target] = pd.to_numeric(df[source], errors="coerce")
    out = out.dropna(subset=["date"]).sort_values("date")
    # Some HAPI collections can expose multiple records per calendar date; the
    # historical input had one row/day.  Particle count/rate and coordinates are
    # collapsed by arithmetic mean, matching a daily-mean reconstruction.
    value_cols = [c for c in out.columns if c != "date"]
    out = out.groupby("date", as_index=False)[value_cols].mean()
    return out, {"time": time_col, "mapping": mapping}


def canonicalise_pa_hourly_speed(df: pd.DataFrame, info: HapiDatasetInfo | dict) -> tuple[pd.DataFrame, dict]:
    time_col = find_time_column(df, info)
    speed_col = resolve_bulk_speed(df, info)
    out = pd.DataFrame({
        "date": pd.to_datetime(df[time_col], errors="coerce", utc=True).dt.tz_localize(None).dt.normalize(),
        "solar_wind_speed_km_s": pd.to_numeric(df[speed_col], errors="coerce"),
    })
    out.loc[out.solar_wind_speed_km_s <= 0, "solar_wind_speed_km_s"] = np.nan
    out = out.dropna(subset=["date"])
    out = out.groupby("date", as_index=False).solar_wind_speed_km_s.mean()
    return out, {"time": time_col, "speed": speed_col, "daily_reducer": "arithmetic mean of finite positive hourly values"}


def assemble_thesis17(cpi_daily: pd.DataFrame, pa_daily_speed: pd.DataFrame | None = None,
                      require_all_particles: bool = True) -> pd.DataFrame:
    out = cpi_daily.copy()
    if pa_daily_speed is not None:
        out = out.merge(pa_daily_speed[["date", "solar_wind_speed_km_s"]], on="date", how="left")
    if "solar_wind_speed_km_s" not in out:
        out["solar_wind_speed_km_s"] = np.nan
    if require_all_particles:
        missing_particles = [c for c in PARTICLE_COLUMNS if c not in out.columns]
        if missing_particles:
            raise ValueError(f"Missing thesis CPI channels after mapping: {missing_particles}")
    for c in THESIS_COLUMNS[1:]:
        if c not in out:
            out[c] = np.nan
    out["time_decimal_year"] = thesis_decimal_year(out["date"])
    return out[[*THESIS_COLUMNS, "date"]].sort_values("date").reset_index(drop=True)



def cpi_core_report(df: pd.DataFrame, spacecraft: str, source_notes: dict | None = None) -> dict:
    """Report the particle-only Phase-A backbone without implying full thesis17 assembly."""
    date = pd.to_datetime(df["date"], errors="coerce")
    fields = {}
    for c in PARTICLE_COLUMNS:
        s = pd.to_numeric(df[c], errors="coerce") if c in df else pd.Series(dtype=float)
        fields[c] = {
            "finite": int(np.isfinite(s).sum()),
            "missing": int(s.isna().sum()),
            "zeros": int((s == 0).sum()),
        }
    report = {
        "spacecraft": spacecraft,
        "product": "CPI daily particle backbone",
        "rows": int(len(df)),
        "start_date": None if date.isna().all() else str(date.min().date()),
        "stop_date": None if date.isna().all() else str(date.max().date()),
        "fields": fields,
        "phase_a_core_ready": all(c in df.columns for c in PARTICLE_COLUMNS),
        "full_thesis17_ready": all(c in df.columns for c in [
            *PARTICLE_COLUMNS, "x_au", "y_au", "z_au", "solar_wind_speed_km_s"
        ]),
        "historical_checks": {
            "p10_reference_rows_from_dissertation": 7042 if spacecraft.lower() == "p10" else None,
            "row_count_matches_p10_reference": bool(len(df) == 7042) if spacecraft.lower() == "p10" else None,
        },
    }
    if source_notes:
        report["sources"] = source_notes
    return report

def thesis17_report(df: pd.DataFrame, spacecraft: str, source_notes: dict | None = None) -> dict:
    date = pd.to_datetime(df["date"], errors="coerce")
    fields = {}
    for c in THESIS_COLUMNS[1:]:
        s = pd.to_numeric(df[c], errors="coerce") if c in df else pd.Series(dtype=float)
        fields[c] = {
            "finite": int(np.isfinite(s).sum()),
            "missing": int(s.isna().sum()),
            "zeros": int((s == 0).sum()),
        }
    report = {
        "spacecraft": spacecraft,
        "rows": int(len(df)),
        "start_date": None if date.isna().all() else str(date.min().date()),
        "stop_date": None if date.isna().all() else str(date.max().date()),
        "fields": fields,
        "historical_checks": {
            "p10_reference_rows_from_dissertation": 7042 if spacecraft.lower() == "p10" else None,
            "row_count_matches_p10_reference": bool(len(df) == 7042) if spacecraft.lower() == "p10" else None,
        },
    }
    if source_notes:
        report["sources"] = source_notes
    return report


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
