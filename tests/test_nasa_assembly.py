import numpy as np
import pandas as pd

from src.acquisition import HapiDatasetInfo
from src.nasa_assembly import (
    resolve_cpi_channels,
    resolve_coordinates,
    resolve_bulk_speed,
    canonicalise_cpi_daily,
    canonicalise_pa_hourly_speed,
    assemble_thesis17,
    thesis_decimal_year,
    thesis17_report,
)
from src.schema import PARTICLE_COLUMNS, THESIS_COLUMNS


def _cpi_info():
    params = [{"name": "Time", "type": "isotime"}]
    names = ["RID2P", "RID2HE", "RID3P", "RID3HE", "RID4P", "RID4HE",
             "RID5P", "RID5HE", "RID5E1", "RID5E2", "RID7", "RID7ZG5"]
    params += [{"name": n, "type": "double", "description": n} for n in names]
    params += [
        {"name": "X", "type": "double", "units": "AU", "description": "heliocentric ecliptic Cartesian X"},
        {"name": "Y", "type": "double", "units": "AU", "description": "heliocentric ecliptic Cartesian Y"},
        {"name": "Z", "type": "double", "units": "AU", "description": "heliocentric ecliptic Cartesian Z"},
    ]
    return HapiDatasetInfo("cpi", None, None, params, {"parameters": params})


def test_cpi_alias_mapping_and_coordinates():
    info = _cpi_info()
    cols = [p["name"] for p in info.parameters]
    df = pd.DataFrame(columns=cols)
    m = resolve_cpi_channels(df, info)
    assert len(m) == 12
    assert m["p_11_20_mev"] == "RID2P"
    assert m["ions_zgt5_gt67_mev_n"] == "RID7ZG5"
    xyz = resolve_coordinates(df, info)
    assert xyz == {"x_au": "X", "y_au": "Y", "z_au": "Z"}


def test_pa_speed_mapping_and_daily_mean():
    params = [
        {"name": "Time", "type": "isotime"},
        {"name": "V", "type": "double", "units": "km/s", "description": "solar wind bulk velocity"},
        {"name": "N", "type": "double", "units": "#/cc", "description": "number density"},
    ]
    info = HapiDatasetInfo("pa", None, None, params, {"parameters": params})
    df = pd.DataFrame({
        "Time": ["1980-01-01T00:00:00Z", "1980-01-01T12:00:00Z", "1980-01-02T00:00:00Z"],
        "V": [300.0, 500.0, 450.0], "N": [5, 6, 7],
    })
    assert resolve_bulk_speed(df, info) == "V"
    daily, meta = canonicalise_pa_hourly_speed(df, info)
    assert daily.loc[0, "solar_wind_speed_km_s"] == 400.0
    assert "arithmetic mean" in meta["daily_reducer"]


def test_end_to_end_17_column_assembly():
    info = _cpi_info()
    n = 3
    raw = pd.DataFrame({"Time": pd.date_range("1980-01-01", periods=n, tz="UTC").astype(str)})
    for name in ["RID2P", "RID2HE", "RID3P", "RID3HE", "RID4P", "RID4HE",
                 "RID5P", "RID5HE", "RID5E1", "RID5E2", "RID7", "RID7ZG5"]:
        raw[name] = np.arange(1, n + 1, dtype=float)
    raw["X"] = [1, 2, 3]
    raw["Y"] = [0, 0, 0]
    raw["Z"] = [0, 0, 0]
    cpi, _ = canonicalise_cpi_daily(raw, info)
    pa = pd.DataFrame({"date": pd.date_range("1980-01-01", periods=n), "solar_wind_speed_km_s": [400, 410, 420]})
    assembled = assemble_thesis17(cpi, pa)
    assert list(assembled.columns[:-1]) == THESIS_COLUMNS
    assert assembled.shape == (3, 18)
    assert np.isclose(assembled.loc[0, "time_decimal_year"], 80 + 1/365)
    report = thesis17_report(assembled, "p10")
    assert report["rows"] == 3
    assert report["historical_checks"]["p10_reference_rows_from_dissertation"] == 7042


def test_thesis_decimal_year_uses_365_in_leap_year():
    x = thesis_decimal_year(pd.to_datetime(["1980-02-29"]))[0]
    assert np.isclose(x, 80 + 60/365)
