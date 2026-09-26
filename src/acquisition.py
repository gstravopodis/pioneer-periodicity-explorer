from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import io
import json
import re
import time
from collections import Counter
from typing import Iterable

import numpy as np
import pandas as pd
import requests

from .provenance import PDS_SOURCES, HAPI_BASE_URL, SPDF_CPI_DAILY
from .schema import PARTICLE_COLUMNS


def download_file(url: str, destination: Path, timeout: int = 120) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with destination.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return destination




def download_spdf_cpi_daily(spacecraft: str, destination: str | Path, timeout: int = 180) -> Path:
    """Download the legacy SPDF 24-hour CPI ASCII file.

    The HAPI conversion of the Pioneer daily CPI collection currently fails on
    legacy two-digit year fields (e.g. ``unable to parse time: 72Z``).  The
    underlying SPDF daily ASCII product remains the most direct source for a
    historical replication, so we support it explicitly rather than modifying
    or fabricating dates client-side.
    """
    sc = spacecraft.lower()
    if sc not in SPDF_CPI_DAILY:
        raise ValueError(f"Unsupported spacecraft: {spacecraft}")
    return download_file(SPDF_CPI_DAILY[sc]["url"], Path(destination), timeout=timeout)


def _numeric_ascii_rows(path: str | Path) -> tuple[np.ndarray, dict]:
    """Read the dominant numeric record layout from a legacy ASCII file.

    Header/documentation lines are ignored.  Fortran D exponents are accepted.
    We select the modal numeric column count so a few numeric-looking header
    lines cannot silently alter the schema.
    """
    parsed: list[list[float]] = []
    widths: list[int] = []
    sample_raw: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line:
            continue
        toks = [t for t in re.split(r"[\s,]+", line) if t]
        vals: list[float] = []
        try:
            for t in toks:
                vals.append(float(t.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
        if len(vals) < 4:
            continue
        parsed.append(vals)
        widths.append(len(vals))
        if len(sample_raw) < 5:
            sample_raw.append(line)
    if not parsed:
        raise ValueError(f"No numeric data records found in {path}")
    mode_width, mode_count = Counter(widths).most_common(1)[0]
    rows = [r for r in parsed if len(r) == mode_width]
    if len(rows) < 10:
        raise ValueError(
            f"Could not identify a stable numeric record layout in {path}; "
            f"modal width={mode_width}, records={len(rows)}"
        )
    arr = np.asarray(rows, dtype=float)
    return arr, {
        "numeric_records_total": len(parsed),
        "modal_width": int(mode_width),
        "modal_records": len(rows),
        "width_histogram": {str(k): int(v) for k, v in sorted(Counter(widths).items())},
        "sample_numeric_lines": sample_raw,
    }


def _expand_legacy_year(y: np.ndarray) -> np.ndarray:
    y = np.rint(y).astype(int)
    out = y.copy()
    mask = (out >= 0) & (out < 100)
    # Pioneer coverage is 1972-1992/95, so the two-digit values are unambiguous.
    out[mask] += 1900
    return out


def _date_from_year_doy(year: np.ndarray, doy: np.ndarray) -> pd.Series:
    years = _expand_legacy_year(year)
    days = np.rint(doy).astype(int)
    if np.any((days < 1) | (days > 366)):
        bad = days[(days < 1) | (days > 366)][:5]
        raise ValueError(f"Invalid day-of-year values in legacy CPI data: {bad.tolist()}")
    text = pd.Series([f"{yy:04d}{dd:03d}" for yy, dd in zip(years, days)])
    return pd.to_datetime(text, format="%Y%j", errors="raise")


def _date_from_thesis_decimal_year(v: np.ndarray) -> pd.Series:
    yy = np.floor(v).astype(int)
    frac = v - yy
    years = _expand_legacy_year(yy)
    # The dissertation encoded YY + DOY/365, including leap years.
    doy = np.rint(frac * 365.0).astype(int)
    doy = np.clip(doy, 1, 366)
    return _date_from_year_doy(years, doy)


def parse_spdf_cpi_daily(path: str | Path, spacecraft: str) -> tuple[pd.DataFrame, dict]:
    """Parse legacy SPDF Pioneer CPI daily ASCII into canonical daily columns.

    Several archival layouts are supported deliberately and *only* when their
    structure is unambiguous:

    - YEAR, DOY, HOUR, 12 CPI channels               (15 columns)
      (the actual legacy ``p10cpi.h24`` layout observed in the archive)
    - YEAR, DOY, 12 CPI channels, X, Y, Z             (17 columns)
    - SCID, YEAR, DOY, 12 CPI channels, X, Y, Z       (18 columns)
    - either coordinate-bearing layout with a trailing solar-wind speed
    - thesis-style decimal-year + 12 channels + X,Y,Z [+ speed]

    The 15-column SPDF daily file is intentionally treated as the *particle
    backbone only*: it contains the 12 CPI channels but not the heliocentric
    coordinates or solar-wind speed that were later merged into the lost
    17-column thesis files.  Those extra fields must come from separate
    official products; they are never inferred from the particle file.

    The function never guesses an arbitrary column window. Unsupported widths
    fail with diagnostics so we can inspect the real archive format.
    """
    arr, diag = _numeric_ascii_rows(path)
    ncol = arr.shape[1]

    def looks_year_doy(year_col: np.ndarray, doy_col: np.ndarray) -> bool:
        y = year_col[np.isfinite(year_col)]
        d = doy_col[np.isfinite(doy_col)]
        if len(y) == 0 or len(d) == 0:
            return False
        y_ok = np.mean(((y >= 70) & (y <= 99)) | ((y >= 1970) & (y <= 2100))) > 0.95
        d_ok = np.mean((d >= 1) & (d <= 366) & (np.abs(d - np.rint(d)) < 1e-6)) > 0.95
        return bool(y_ok and d_ok)

    def looks_hour(hour_col: np.ndarray) -> bool:
        h = hour_col[np.isfinite(hour_col)]
        if len(h) == 0:
            return False
        return bool(np.mean((h >= 0) & (h <= 23) & (np.abs(h - np.rint(h)) < 1e-6)) > 0.95)

    start = None
    date: pd.Series
    layout: str
    if ncol >= 15 and looks_year_doy(arr[:, 0], arr[:, 1]):
        # The official SPDF 24-hour files use YY, DOY, HOUR followed by the
        # 12 CPI channels.  Keep the HOUR in diagnostics, but daily replication
        # uses the calendar date exactly as the 1998 MATLAB vector did.
        if ncol == 15 and looks_hour(arr[:, 2]):
            start = 3
            date = _date_from_year_doy(arr[:, 0], arr[:, 1])
            layout = "year_doy_hour_particles"
            diag["hour_values"] = sorted({int(x) for x in np.rint(arr[:, 2]).astype(int)})[:24]
        else:
            start = 2
            date = _date_from_year_doy(arr[:, 0], arr[:, 1])
            layout = "year_doy"
    elif ncol >= 18 and np.nanstd(arr[:, 0]) < 1e-9 and looks_year_doy(arr[:, 1], arr[:, 2]):
        start = 3
        date = _date_from_year_doy(arr[:, 1], arr[:, 2])
        layout = "scid_year_doy"
    else:
        # Thesis-style decimal year is typically 72.x ... 92.x.
        c0 = arr[:, 0]
        finite = c0[np.isfinite(c0)]
        dec_ok = len(finite) and np.mean(((finite >= 70) & (finite < 100)) | ((finite >= 1970) & (finite < 2100))) > 0.95
        if not dec_ok:
            raise ValueError(
                f"Unsupported SPDF CPI layout for {spacecraft}: {ncol} columns and no "
                "recognisable YEAR/DOY or decimal-year time field. "
                f"Diagnostics: {diag}"
            )
        start = 1
        date = _date_from_thesis_decimal_year(c0)
        layout = "decimal_year"

    remaining = ncol - start
    if remaining not in {12, 15, 16}:
        raise ValueError(
            f"Unsupported SPDF CPI layout for {spacecraft}: {ncol} columns, time layout "
            f"{layout}, leaving {remaining} value columns; expected 12 (particle-only), "
            f"15 (12 channels + XYZ), or 16 (+ solar-wind speed). Diagnostics: {diag}"
        )

    vals = arr[:, start:]
    out = pd.DataFrame({"date": pd.to_datetime(date).dt.normalize()})
    for i, col in enumerate(PARTICLE_COLUMNS):
        out[col] = vals[:, i]

    has_coordinates = remaining >= 15
    if has_coordinates:
        out["x_au"] = vals[:, 12]
        out["y_au"] = vals[:, 13]
        out["z_au"] = vals[:, 14]
    has_speed = remaining == 16
    if has_speed:
        out["solar_wind_speed_km_s"] = vals[:, 15]

    out = out.sort_values("date").drop_duplicates("date", keep="first").reset_index(drop=True)
    diag.update({
        "layout": layout,
        "value_columns_after_time": int(remaining),
        "contains_coordinates": bool(has_coordinates),
        "contains_solar_wind_speed": bool(has_speed),
        "particle_channels": list(PARTICLE_COLUMNS),
        "rows_after_date_deduplication": int(len(out)),
        "start_date": str(out.date.min().date()),
        "stop_date": str(out.date.max().date()),
    })
    return out, diag


def download_pds_manifests(out_dir: str | Path) -> list[Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved = []
    for key, meta in PDS_SOURCES.items():
        slug = "p10" if "10" in key else "p11"
        saved.append(download_file(meta["manifest"], out / f"{slug}_collection_1day.csv"))
    return saved


def parse_collection_manifest(path: str | Path) -> list[dict]:
    """Parse a PDS4 collection inventory defensively."""
    rows = []
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as f:
        for raw in csv.reader(f):
            if not raw:
                continue
            item = {"raw": raw}
            if len(raw) >= 2:
                item["status"] = raw[0].strip()
                item["lidvid"] = raw[1].strip()
            elif len(raw) == 1:
                item["lidvid"] = raw[0].strip()
            rows.append(item)
    return rows


@dataclass
class HapiDatasetInfo:
    dataset_id: str
    start_date: str | None
    stop_date: str | None
    parameters: list[dict]
    raw: dict

    @property
    def parameter_names(self) -> list[str]:
        return [p.get("name", "") for p in self.parameters]


class HapiRequestError(RuntimeError):
    """A HAPI HTTP/API failure with enough detail to diagnose the server."""

    def __init__(self, message: str, *, status_code: int | None = None,
                 url: str | None = None, response_text: str | None = None):
        super().__init__(message)
        self.status_code = status_code
        self.url = url
        self.response_text = response_text


class HapiClient:
    """Small, dependency-light client for the PDS/PPI HAPI server.

    PDS/PPI currently advertises HAPI 3.1. HAPI 3 renamed request parameters
    ``id`` -> ``dataset`` and ``time.min/time.max`` -> ``start/stop``. The
    specification says HAPI 3 servers should still accept the old names, but in
    practice some servers are less reliable through the compatibility path.
    Therefore this client prefers native HAPI-3 syntax and falls back to the
    older syntax only when needed.
    """

    def __init__(self, base_url: str = HAPI_BASE_URL, timeout: int = 120,
                 session: requests.Session | None = None, retries: int = 2):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.session = session or requests.Session()
        self.retries = max(0, int(retries))
        self._major_version: int | None = None

    @staticmethod
    def _body_preview(response: requests.Response, limit: int = 1200) -> str:
        try:
            body = response.text.strip().replace("\x00", "")
        except Exception:
            return ""
        return body[:limit]

    def _get(self, endpoint: str, params: dict | None = None) -> requests.Response:
        url = f"{self.base_url}/{endpoint}"
        last_exc: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                r = self.session.get(url, params=params, timeout=self.timeout)
            except requests.RequestException as exc:
                last_exc = exc
                if attempt < self.retries:
                    time.sleep(0.5 * (2 ** attempt))
                    continue
                raise HapiRequestError(
                    f"Network error while contacting PDS/PPI HAPI at {self.base_url}: {exc}",
                    url=url,
                ) from exc

            # Retry transient 502/503/504, but not deterministic 500 responses.
            if r.status_code in {502, 503, 504} and attempt < self.retries:
                time.sleep(0.5 * (2 ** attempt))
                continue
            return r
        raise HapiRequestError(f"HAPI request failed: {last_exc}", url=url)

    def capabilities(self) -> dict:
        r = self._get("capabilities")
        try:
            r.raise_for_status()
            payload = r.json()
        except requests.RequestException as exc:
            raise HapiRequestError(
                f"HAPI capabilities request failed with HTTP {r.status_code}.",
                status_code=r.status_code,
                url=r.url,
                response_text=self._body_preview(r),
            ) from exc
        except ValueError as exc:
            raise HapiRequestError(
                "HAPI capabilities response was not valid JSON.",
                status_code=r.status_code, url=r.url, response_text=self._body_preview(r),
            ) from exc
        version = str(payload.get("HAPI", ""))
        try:
            self._major_version = int(version.split(".", 1)[0])
        except Exception:
            self._major_version = None
        return payload

    def _prefer_v3(self) -> bool:
        if self._major_version is None:
            try:
                self.capabilities()
            except HapiRequestError:
                # If capabilities itself is unavailable, prefer current syntax.
                self._major_version = 3
        return (self._major_version or 3) >= 3

    def _json_request(self, endpoint: str, candidates: list[dict]) -> dict:
        errors: list[HapiRequestError] = []
        for params in candidates:
            r = self._get(endpoint, params=params)
            if r.ok:
                try:
                    payload = r.json()
                except ValueError as exc:
                    err = HapiRequestError(
                        f"HAPI {endpoint} response was not valid JSON.",
                        status_code=r.status_code, url=r.url,
                        response_text=self._body_preview(r),
                    )
                    errors.append(err)
                    continue
                status = payload.get("status", {})
                if status and status.get("code", 1200) != 1200:
                    err = HapiRequestError(
                        f"HAPI {endpoint} error {status}", status_code=r.status_code,
                        url=r.url, response_text=self._body_preview(r),
                    )
                    errors.append(err)
                    continue
                return payload

            errors.append(HapiRequestError(
                f"HAPI {endpoint} returned HTTP {r.status_code}.",
                status_code=r.status_code, url=r.url,
                response_text=self._body_preview(r),
            ))

        err = errors[-1]
        detail = f" URL: {err.url}" if err.url else ""
        body = f" Response: {err.response_text}" if err.response_text else ""
        raise HapiRequestError(f"{err}{detail}{body}", status_code=err.status_code,
                               url=err.url, response_text=err.response_text)

    def catalog(self) -> list[dict]:
        r = self._get("catalog")
        try:
            r.raise_for_status()
            payload = r.json()
        except requests.RequestException as exc:
            raise HapiRequestError(
                f"HAPI catalog request failed with HTTP {r.status_code}.",
                status_code=r.status_code, url=r.url, response_text=self._body_preview(r),
            ) from exc
        return payload.get("catalog", [])

    def info(self, dataset_id: str) -> HapiDatasetInfo:
        # Native HAPI 3 first; old HAPI 1/2 compatibility syntax second.
        if self._prefer_v3():
            candidates = [{"dataset": dataset_id}, {"id": dataset_id}]
        else:
            candidates = [{"id": dataset_id}, {"dataset": dataset_id}]
        raw = self._json_request("info", candidates)
        return HapiDatasetInfo(
            dataset_id=dataset_id,
            start_date=raw.get("startDate"),
            stop_date=raw.get("stopDate"),
            parameters=raw.get("parameters", []),
            raw=raw,
        )

    def _data_candidates(self, dataset_id: str, start: str, stop: str,
                         parameters: Iterable[str] | None = None) -> list[dict]:
        selected = list(parameters) if parameters else None
        v3 = {"dataset": dataset_id, "start": start, "stop": stop, "format": "csv"}
        v2 = {"id": dataset_id, "time.min": start, "time.max": stop, "format": "csv"}
        if selected:
            value = ",".join(selected)
            v3["parameters"] = value
            v2["parameters"] = value
        return [v3, v2] if self._prefer_v3() else [v2, v3]

    def data(self, dataset_id: str, start: str, stop: str,
             parameters: Iterable[str] | None = None,
             info: HapiDatasetInfo | None = None) -> pd.DataFrame:
        info = info or self.info(dataset_id)
        selected = list(parameters) if parameters else info.parameter_names
        if not selected:
            raise ValueError(f"HAPI info returned no parameters for {dataset_id}")

        errors: list[HapiRequestError] = []
        for q in self._data_candidates(dataset_id, start, stop, parameters):
            r = self._get("data", params=q)
            if r.ok:
                text = r.text.strip()
                if not text:
                    return pd.DataFrame(columns=selected)
                clean = "\n".join(line for line in text.splitlines() if line and not line.startswith("#"))
                try:
                    return pd.read_csv(io.StringIO(clean), header=None, names=selected)
                except Exception as exc:
                    raise HapiRequestError(
                        f"HAPI CSV parsing failed for {dataset_id} ({start} to {stop}): {exc}",
                        status_code=r.status_code, url=r.url,
                        response_text=clean[:1200],
                    ) from exc

            errors.append(HapiRequestError(
                f"HAPI data returned HTTP {r.status_code} for {dataset_id} "
                f"({start} to {stop}).",
                status_code=r.status_code, url=r.url,
                response_text=self._body_preview(r),
            ))

        err = errors[-1]
        attempts = "\n".join(
            f"  - HTTP {e.status_code}: {e.url}" +
            (f"\n    {e.response_text}" if e.response_text else "")
            for e in errors
        )
        raise HapiRequestError(
            "PDS/PPI HAPI rejected both native HAPI-3 and legacy request syntax.\n" + attempts,
            status_code=err.status_code, url=err.url, response_text=err.response_text,
        )

    def data_year_chunks(self, dataset_id: str, start: str, stop: str,
                         parameters: Iterable[str] | None = None) -> pd.DataFrame:
        """Fetch a long mission in year-sized requests and concatenate it."""
        info = self.info(dataset_id)
        s = pd.Timestamp(start)
        e = pd.Timestamp(stop)
        frames: list[pd.DataFrame] = []
        cur = s
        while cur < e:
            nxt = min(pd.Timestamp(year=cur.year + 1, month=1, day=1), e)
            start_iso = cur.strftime("%Y-%m-%dT%H:%M:%SZ")
            stop_iso = nxt.strftime("%Y-%m-%dT%H:%M:%SZ")
            part = self.data(
                dataset_id,
                start_iso,
                stop_iso,
                parameters=parameters,
                info=info,
            )
            if not part.empty:
                frames.append(part)
            cur = nxt
        if not frames:
            names = list(parameters) if parameters else info.parameter_names
            return pd.DataFrame(columns=names)
        return pd.concat(frames, ignore_index=True)


def save_hapi_info(info: HapiDatasetInfo, destination: str | Path) -> Path:
    p = Path(destination)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(info.raw, indent=2), encoding="utf-8")
    return p
