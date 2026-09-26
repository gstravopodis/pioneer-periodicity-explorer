"""Acquire official Pioneer products and rebuild the thesis-style 17-column files.

The preferred CPI path is the legacy NASA/SPDF 24-hour ASCII product.  This is
not a workaround of convenience: it is the archival product closest to the
input used in the 1998 dissertation.  The PDS/PPI HAPI metadata endpoint is
still used when available, but its data endpoint currently fails on the
Pioneer daily files because the backend tries to parse legacy two-digit years
(``unable to parse time: 72Z``).

Solar-wind speed is obtained from the CPI daily file when that field is present.
Otherwise the script uses the official PDS/PPI HAPI Plasma Analyzer hourly
product and collapses valid positive speeds to a daily arithmetic mean.

Examples
--------
python scripts/fetch_nasa.py --spacecraft p10
python scripts/fetch_nasa.py --spacecraft both
python scripts/fetch_nasa.py --spacecraft both --cpi-source spdf
python scripts/fetch_nasa.py --spacecraft p10 --cpi-source hapi

Outputs are written under data/nasa/{raw,metadata,processed}/ by default.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.acquisition import (
    HapiClient,
    HapiRequestError,
    download_spdf_cpi_daily,
    parse_spdf_cpi_daily,
    save_hapi_info,
)
from src.nasa_assembly import (
    canonicalise_cpi_daily,
    canonicalise_pa_hourly_speed,
    assemble_thesis17,
    thesis17_report,
    cpi_core_report,
    sha256_file,
)
from src.provenance import HAPI_DATASETS, SPDF_CPI_DAILY


def _safe_info(client: HapiClient, dataset_id: str, path: Path, label: str) -> dict | None:
    """Fetch HAPI metadata without making metadata availability a hard dependency."""
    try:
        print(f"{label} HAPI info")
        info = client.info(dataset_id)
        save_hapi_info(info, path)
        return info
    except Exception as exc:
        print(f"{label} HAPI metadata unavailable; continuing with archive source: {exc}")
        return None


def _acquire_cpi(
    sc: str,
    raw_dir: Path,
    meta_dir: Path,
    client: HapiClient,
    source: str,
    timeout: int,
) -> tuple[object, dict]:
    """Return canonical daily CPI dataframe and provenance block."""
    cpi_meta = HAPI_DATASETS[sc]["cpi_daily"]
    cpi_info = _safe_info(
        client,
        cpi_meta["id"],
        meta_dir / f"{sc}_cpi_daily_info.json",
        f"[{sc}] CPI daily:",
    )

    errors: list[str] = []

    if source in {"auto", "spdf"}:
        raw_path = raw_dir / f"{sc}_cpi_daily_spdf.h24"
        try:
            print(f"[{sc}] downloading legacy SPDF CPI daily ASCII")
            download_spdf_cpi_daily(sc, raw_path, timeout=timeout)
            cpi, parse_diag = parse_spdf_cpi_daily(raw_path, sc)
            provenance = {
                **cpi_meta,
                "transport": "NASA/SPDF legacy 24-hour ASCII",
                "url": SPDF_CPI_DAILY[sc]["url"],
                "archive_dir": SPDF_CPI_DAILY[sc]["archive_dir"],
                "note": SPDF_CPI_DAILY[sc]["note"],
                "raw_sha256": sha256_file(raw_path),
                "parser": parse_diag,
                "hapi_metadata_saved": cpi_info is not None,
            }
            print(
                f"[{sc}] SPDF CPI parsed: {len(cpi)} daily rows; "
                f"{parse_diag.get('start_date')} -> {parse_diag.get('stop_date')}; "
                f"layout={parse_diag.get('layout')}"
            )
            return cpi, provenance
        except Exception as exc:
            errors.append(f"SPDF CPI: {type(exc).__name__}: {exc}")
            if source == "spdf":
                raise RuntimeError(
                    f"[{sc}] Could not acquire the requested SPDF CPI source. "
                    f"{errors[-1]}"
                ) from exc
            print(f"[{sc}] SPDF CPI path failed; trying HAPI data as fallback.\n  {errors[-1]}")

    if source in {"auto", "hapi"}:
        if cpi_info is None:
            try:
                cpi_info = client.info(cpi_meta["id"])
                save_hapi_info(cpi_info, meta_dir / f"{sc}_cpi_daily_info.json")
            except Exception as exc:
                errors.append(f"HAPI CPI metadata: {type(exc).__name__}: {exc}")
                raise RuntimeError(
                    f"[{sc}] CPI acquisition failed from all permitted sources. "
                    + " | ".join(errors)
                ) from exc
        try:
            print(f"[{sc}] downloading CPI daily through HAPI {cpi_meta['coverage'][0]} -> {cpi_meta['coverage'][1]}")
            raw = client.data_year_chunks(cpi_meta["id"], *cpi_meta["coverage"])
            raw_path = raw_dir / f"{sc}_cpi_daily_hapi.csv"
            raw.to_csv(raw_path, index=False)
            cpi, mapping = canonicalise_cpi_daily(raw, cpi_info)
            provenance = {
                **cpi_meta,
                "transport": "PDS/PPI HAPI",
                "mapping": mapping,
                "raw_sha256": sha256_file(raw_path),
            }
            return cpi, provenance
        except Exception as exc:
            errors.append(f"HAPI CPI: {type(exc).__name__}: {exc}")
            raise RuntimeError(
                f"[{sc}] CPI acquisition failed from all permitted sources. "
                + " | ".join(errors)
            ) from exc

    raise ValueError(f"Unsupported CPI source: {source}")


def _acquire_pa_daily_speed(
    sc: str,
    raw_dir: Path,
    meta_dir: Path,
    client: HapiClient,
) -> tuple[object, dict]:
    pa_meta = HAPI_DATASETS[sc]["pa_hourly"]
    print(f"[{sc}] HAPI info: PA hourly")
    pa_info = client.info(pa_meta["id"])
    save_hapi_info(pa_info, meta_dir / f"{sc}_pa_hourly_info.json")
    print(f"[{sc}] downloading PA hourly {pa_meta['coverage'][0]} -> {pa_meta['coverage'][1]}")
    pa_raw = client.data_year_chunks(pa_meta["id"], *pa_meta["coverage"])
    pa_raw_path = raw_dir / f"{sc}_pa_hourly_hapi.csv"
    pa_raw.to_csv(pa_raw_path, index=False)
    pa, mapping = canonicalise_pa_hourly_speed(pa_raw, pa_info)
    provenance = {
        **pa_meta,
        "transport": "PDS/PPI HAPI",
        "mapping": mapping,
        "raw_sha256": sha256_file(pa_raw_path),
    }
    return pa, provenance


def fetch_one(sc: str, out: Path, client: HapiClient, cpi_source: str, timeout: int) -> None:
    meta_dir = out / "metadata"
    raw_dir = out / "raw"
    proc_dir = out / "processed"
    for d in (meta_dir, raw_dir, proc_dir):
        d.mkdir(parents=True, exist_ok=True)

    cpi, cpi_source_meta = _acquire_cpi(sc, raw_dir, meta_dir, client, cpi_source, timeout)

    # Phase A core does not require geometry or plasma speed.  Save the exact
    # daily particle backbone first so the historical ACF/CCF experiment can
    # proceed even when modern metadata/trajectory services are unavailable.
    core_path = proc_dir / f"{sc}_cpi_daily.csv"
    cpi.to_csv(core_path, index=False)
    core_report = cpi_core_report(cpi, sc, {"cpi": cpi_source_meta})
    core_report["processed"] = {
        "cpi_daily_csv": str(core_path),
        "cpi_daily_csv_sha256": sha256_file(core_path),
    }
    core_report_path = proc_dir / f"{sc}_cpi_report.json"
    core_report_path.write_text(json.dumps(core_report, indent=2), encoding="utf-8")
    print(
        f"[{sc}] Phase-A particle backbone ready: {len(cpi)} rows; "
        f"historical P10 row match={core_report['historical_checks'].get('row_count_matches_p10_reference')}"
    )

    # A full reconstruction is optional at this stage.  Only create a file
    # called thesis17 when all four non-particle inputs are genuinely present.
    required_extra = {"x_au", "y_au", "z_au", "solar_wind_speed_km_s"}
    if required_extra.issubset(cpi.columns):
        assembled = assemble_thesis17(cpi, None)
        pa_source_meta = {
            "transport": "not required",
            "note": "Coordinates and solar-wind speed were already present in the CPI source.",
        }
    else:
        missing = sorted(required_extra - set(cpi.columns))
        print(
            f"[{sc}] Full 17-column reconstruction deferred; missing separate official fields: "
            + ", ".join(missing)
        )
        print(
            f"[{sc}] This does NOT block the locked ~470 d ACF and ~340 d unshifted P10/P11 CCF replication."
        )
        return

    csv_path = proc_dir / f"{sc}_thesis17.csv"
    dat_path = proc_dir / f"{sc}cpi_reconstructed.dat"
    assembled.to_csv(csv_path, index=False)
    assembled.drop(columns=["date"]).to_csv(
        dat_path,
        sep=" ",
        header=False,
        index=False,
        na_rep="0",
        float_format="%.10g",
    )

    sources = {
        "cpi": cpi_source_meta,
        "pa": pa_source_meta,
        "assembly_note": (
            "Full thesis17 export is emitted only when particle channels, XYZ coordinates and "
            "solar-wind speed are all present from official sources. No missing field is fabricated."
        ),
    }
    report = thesis17_report(assembled, sc, sources)
    report["processed"] = {
        "csv": str(csv_path),
        "csv_sha256": sha256_file(csv_path),
        "dat": str(dat_path),
        "dat_sha256": sha256_file(dat_path),
    }
    report_path = proc_dir / f"{sc}_assembly_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spacecraft", choices=["p10", "p11", "both"], default="both")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "nasa")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument(
        "--cpi-source",
        choices=["auto", "spdf", "hapi"],
        default="auto",
        help=(
            "CPI transport. auto prefers the legacy NASA/SPDF daily ASCII and falls back to HAPI; "
            "spdf forces the historical archive; hapi forces the PDS/PPI HAPI data endpoint."
        ),
    )
    args = ap.parse_args()
    client = HapiClient(timeout=args.timeout)
    try:
        caps = client.capabilities()
        print(f"PDS/PPI HAPI version: {caps.get('HAPI', 'unknown')}; formats: {caps.get('outputFormats', [])}")
    except HapiRequestError as exc:
        print(f"PDS/PPI HAPI capabilities unavailable: {exc}. SPDF CPI acquisition can still proceed.")

    todo = ["p10", "p11"] if args.spacecraft == "both" else [args.spacecraft]
    for sc in todo:
        fetch_one(sc, args.out, client, args.cpi_source, args.timeout)


if __name__ == "__main__":
    main()
