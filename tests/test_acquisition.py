from pathlib import Path
import pandas as pd

from src.acquisition import parse_collection_manifest, HapiClient


def test_manifest_parser(tmp_path: Path):
    p = tmp_path / "collection.csv"
    p.write_text("P,urn:nasa:pds:test:a::1.0\nS,urn:nasa:pds:test:b::1.0\n", encoding="utf-8")
    rows = parse_collection_manifest(p)
    assert len(rows) == 2
    assert rows[0]["status"] == "P"
    assert rows[1]["lidvid"].endswith("b::1.0")


class _FakeResponse:
    def __init__(self, payload=None, text="", status_code=200, url="https://fake"):
        self._payload = payload
        self.text = text
        self.status_code = status_code
        self.url = url

    @property
    def ok(self):
        return 200 <= self.status_code < 400

    def raise_for_status(self):
        if not self.ok:
            import requests
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class _FakeSession:
    def __init__(self, fail_v3_data=False):
        self.calls = []
        self.fail_v3_data = fail_v3_data

    def get(self, url, params=None, timeout=None, stream=False):
        self.calls.append((url, params))
        if url.endswith("/capabilities"):
            return _FakeResponse({"HAPI": "3.1", "status": {"code": 1200, "message": "OK"}}, url=url)
        if url.endswith("/info"):
            return _FakeResponse({
                "HAPI": "3.1",
                "status": {"code": 1200, "message": "OK"},
                "startDate": "1980-01-01T00:00:00Z",
                "stopDate": "1980-01-03T00:00:00Z",
                "parameters": [
                    {"name": "Time", "type": "isotime"},
                    {"name": "RID2P", "type": "double"},
                ],
            }, url=url)
        if url.endswith("/data"):
            if self.fail_v3_data and params and "dataset" in params:
                return _FakeResponse(text="internal error", status_code=500, url=url + "?v3")
            return _FakeResponse(text="1980-01-01T00:00:00Z,1.25\n1980-01-02T00:00:00Z,2.5\n", url=url)
        raise AssertionError(url)


def test_hapi_client_prefers_hapi3_parameter_names():
    s = _FakeSession()
    c = HapiClient(session=s)
    info = c.info("test")
    df = c.data("test", "1980-01-01T00:00:00Z", "1980-01-03T00:00:00Z", info=info)
    assert list(df.columns) == ["Time", "RID2P"]
    assert df.loc[1, "RID2P"] == 2.5
    info_call = [x for x in s.calls if x[0].endswith("/info")][0]
    assert info_call[1]["dataset"] == "test"
    data_call = [x for x in s.calls if x[0].endswith("/data")][0]
    assert data_call[1]["dataset"] == "test"
    assert data_call[1]["start"].startswith("1980-01-01")
    assert data_call[1]["stop"].startswith("1980-01-03")


def test_hapi_client_falls_back_to_legacy_data_syntax():
    s = _FakeSession(fail_v3_data=True)
    c = HapiClient(session=s)
    info = c.info("test")
    df = c.data("test", "1980-01-01T00:00:00Z", "1980-01-03T00:00:00Z", info=info)
    assert len(df) == 2
    data_calls = [x for x in s.calls if x[0].endswith("/data")]
    assert "dataset" in data_calls[0][1]
    assert "id" in data_calls[1][1]
    assert data_calls[1][1]["time.min"].startswith("1980-01-01")


def _write_numeric_rows(path: Path, rows):
    path.write_text("header ignored\n" + "\n".join(" ".join(str(x) for x in row) for row in rows) + "\n", encoding="utf-8")


def test_parse_spdf_year_doy_17_columns(tmp_path: Path):
    from src.acquisition import parse_spdf_cpi_daily
    rows = []
    # 2 time + 12 channels + XYZ = 17 numeric columns.
    for i in range(12):
        rows.append([72, 63 + i] + [100 + 10*j + i for j in range(12)] + [1.0+i/100, 2.0, 3.0])
    p = tmp_path / "p10cpi.h24"
    _write_numeric_rows(p, rows)
    df, diag = parse_spdf_cpi_daily(p, "p10")
    assert len(df) == 12
    assert str(df.loc[0, "date"].date()) == "1972-03-03"
    assert df.loc[0, "p_11_20_mev"] == 100
    assert df.loc[0, "ions_zgt5_gt67_mev_n"] == 210
    assert df.loc[0, "x_au"] == 1.0
    assert "solar_wind_speed_km_s" not in df.columns
    assert diag["layout"] == "year_doy"


def test_parse_spdf_scid_year_doy_18_columns(tmp_path: Path):
    from src.acquisition import parse_spdf_cpi_daily
    rows = []
    # SCID + 2 time + 12 channels + XYZ = 18 columns.
    for i in range(12):
        rows.append([10, 72, 63 + i] + [i + j/10 for j in range(12)] + [1.0, 0.0, -1.0])
    p = tmp_path / "p10_scid.h24"
    _write_numeric_rows(p, rows)
    df, diag = parse_spdf_cpi_daily(p, "p10")
    assert len(df) == 12
    assert diag["layout"] == "scid_year_doy"
    assert df.loc[3, "p_11_20_mev"] == 3.0


def test_parse_spdf_decimal_year_with_speed(tmp_path: Path):
    from src.acquisition import parse_spdf_cpi_daily
    rows = []
    # decimal year + 12 channels + XYZ + speed = 17 columns.
    for i in range(12):
        decimal_year = 72 + (63 + i) / 365.0
        rows.append([decimal_year] + [1.0 + j + i/100 for j in range(12)] + [1.2, -0.1, 0.05, 410+i])
    p = tmp_path / "thesis_style.dat"
    _write_numeric_rows(p, rows)
    df, diag = parse_spdf_cpi_daily(p, "p10")
    assert len(df) == 12
    assert diag["layout"] == "decimal_year"
    assert diag["contains_solar_wind_speed"] is True
    assert df.loc[0, "solar_wind_speed_km_s"] == 410


def test_parse_real_spdf_15col_year_doy_hour_particle_layout(tmp_path: Path):
    """Regression for the actual p10cpi.h24 record shape observed on 2026-09-25."""
    from src.acquisition import parse_spdf_cpi_daily

    lines = [
        "72 63 0 1.74E-05 1.91E-04 6.09E-04 9.40E-04 0.00E+00 3.48E-05 1.57E-04 2.96E-04 1.13E-03 6.61E-04 5.06E-01 1.90E-03",
        "72 64 0 8.73E-05 2.91E-04 7.42E-04 1.13E-03 2.91E-05 1.45E-05 2.62E-04 3.20E-04 1.27E-03 1.06E-03 4.79E-01 1.73E-03",
    ]
    # Need >=10 stable records for the conservative parser.
    p = tmp_path / "p10cpi.h24"
    p.write_text("\n".join(lines * 6) + "\n", encoding="utf-8")
    df, diag = parse_spdf_cpi_daily(p, "p10")
    # Duplicated dates are deliberately deduplicated after parsing.
    assert len(df) == 2
    assert str(df.loc[0, "date"].date()) == "1972-03-03"
    assert diag["layout"] == "year_doy_hour_particles"
    assert diag["value_columns_after_time"] == 12
    assert diag["contains_coordinates"] is False
    assert diag["contains_solar_wind_speed"] is False
    assert diag["hour_values"] == [0]
    assert abs(df.loc[0, "p_11_20_mev"] - 1.74e-05) < 1e-12
    assert abs(df.loc[0, "ions_zgt5_gt67_mev_n"] - 1.90e-03) < 1e-12
    assert "x_au" not in df.columns
