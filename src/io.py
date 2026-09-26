from __future__ import annotations

from pathlib import Path
import io
import numpy as np
import pandas as pd

from .schema import THESIS_COLUMNS, NASA_FIELD_ALIASES, P11_AFFECTED_COLUMNS


def decimal_year_1998_to_datetime(values: pd.Series) -> pd.Series:
    """Convert the dissertation's YY.fraction timestamps to calendar dates.

    The dissertation explicitly gives 72.1721 as the 63rd day of 1972 and states
    that the fractional part is based on day/365.  Replication mode therefore
    deliberately uses 365 even in leap years.  Day 1 maps to Jan 1, hence the
    ``-1`` when converting ordinal day to a zero-based timedelta.
    """
    x = pd.to_numeric(values, errors="coerce")
    yy = np.floor(x).astype("Int64")
    frac = x - yy.astype(float)
    year = yy.where(yy >= 100, 1900 + yy)
    ordinal = np.rint(frac * 365.0).astype("Int64")
    ordinal = ordinal.clip(lower=1, upper=365)
    base = pd.to_datetime(year.astype(str) + "-01-01", errors="coerce")
    return base + pd.to_timedelta(ordinal - 1, unit="D")


def _read_text_source(source) -> str:
    if hasattr(source, "read"):
        raw = source.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return raw
    return Path(source).read_text(encoding="utf-8", errors="replace")


def load_thesis_dat(source) -> pd.DataFrame:
    """Load the original 17-column p10cpi.dat / p11cpi.dat layout."""
    raw = _read_text_source(source)
    df = pd.read_csv(io.StringIO(raw), sep=r"\s+", header=None, comment="#")
    if df.shape[1] != 17:
        raise ValueError(f"Expected 17 columns from the thesis schema; found {df.shape[1]}.")
    df.columns = THESIS_COLUMNS
    df["date"] = decimal_year_1998_to_datetime(df["time_decimal_year"])
    return df


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    for c in df.columns:
        key = str(c).strip().upper().replace(" ", "_")
        if key in NASA_FIELD_ALIASES:
            rename[c] = NASA_FIELD_ALIASES[key]
    return df.rename(columns=rename)


def _derive_date(df: pd.DataFrame) -> pd.Series:
    # Flexible handling for OMNI/SPDF exports and ordinary CSVs.
    candidates = ["date", "datetime", "time", "epoch", "timestamp"]
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c in lower:
            return pd.to_datetime(df[lower[c]], errors="coerce")

    if "time_decimal_year" in df.columns:
        return decimal_year_1998_to_datetime(df["time_decimal_year"])

    year_col = next((c for c in df.columns if str(c).strip().upper() in {"YEAR", "YYYY"}), None)
    doy_col = next((c for c in df.columns if str(c).strip().upper() in {"DOY", "DAY", "DAY_OF_YEAR"}), None)
    if year_col is not None and doy_col is not None:
        y = pd.to_numeric(df[year_col], errors="coerce").astype("Int64")
        d = pd.to_numeric(df[doy_col], errors="coerce").astype("Int64")
        return pd.to_datetime(y.astype(str) + d.astype(str).str.zfill(3), format="%Y%j", errors="coerce")

    raise ValueError("Could not infer dates. Provide date/datetime, time_decimal_year, or YEAR+DOY columns.")


def load_csv(source) -> pd.DataFrame:
    df = pd.read_csv(source, comment="#")
    df = _normalise_columns(df)
    df["date"] = _derive_date(df)
    return df


def load_auto(source, filename: str | None = None) -> pd.DataFrame:
    """Load thesis whitespace data or a NASA/ordinary CSV with minimal guessing."""
    name = (filename or getattr(source, "name", "") or "").lower()
    if name.endswith(".csv"):
        return load_csv(source)
    # For .dat/.txt first try the exact 17-column historical layout.
    try:
        return load_thesis_dat(source)
    except Exception:
        if hasattr(source, "seek"):
            source.seek(0)
        return load_csv(source)


def apply_p11_cpi_quality_mask(df: pd.DataFrame) -> pd.DataFrame:
    """Mask thesis columns 6–13 from day 239 of 1980 onward, as documented."""
    out = df.copy()
    cutoff = pd.Timestamp("1980-08-26")  # day 239 in leap-year 1980
    mask = out["date"] >= cutoff
    for c in P11_AFFECTED_COLUMNS:
        if c in out.columns:
            out.loc[mask, c] = np.nan
    return out


def data_quality_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for c in [x for x in THESIS_COLUMNS[1:] if x in df.columns]:
        s = pd.to_numeric(df[c], errors="coerce")
        rows.append({
            "field": c,
            "rows": len(s),
            "finite": int(np.isfinite(s).sum()),
            "missing": int(s.isna().sum()),
            "zeros": int((s == 0).sum()),
        })
    return pd.DataFrame(rows)
