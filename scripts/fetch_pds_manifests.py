"""Fetch the official PDS4 daily CPI collection inventories for Pioneer 10/11.

Run locally on a machine with internet access:
    python scripts/fetch_pds_manifests.py
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.acquisition import download_pds_manifests, parse_collection_manifest

out = ROOT / "data" / "pds_manifests"
for path in download_pds_manifests(out):
    members = parse_collection_manifest(path)
    print(f"saved {path.name}: {len(members)} inventory rows")
