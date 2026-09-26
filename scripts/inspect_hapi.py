"""Inspect the official PDS/PPI HAPI parameter metadata before assembly."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.acquisition import HapiClient
from src.provenance import HAPI_DATASETS


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spacecraft", choices=["p10", "p11"], required=True)
    ap.add_argument("--dataset", choices=["cpi_daily", "cpi_hourly", "pa_hourly", "pa_trajectory"], default="cpi_daily")
    args = ap.parse_args()
    spec = HAPI_DATASETS[args.spacecraft][args.dataset]
    info = HapiClient().info(spec["id"])
    print(json.dumps({
        "id": spec["id"],
        "doi": spec["doi"],
        "startDate": info.start_date,
        "stopDate": info.stop_date,
        "parameters": info.parameters,
    }, indent=2))


if __name__ == "__main__":
    main()
