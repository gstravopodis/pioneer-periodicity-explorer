"""One-command NASA acquisition -> historical reconstruction -> Phase-A replication."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="p_11_20_mev")
    ap.add_argument("--skip-download", action="store_true", help="Use already assembled data/nasa/processed files")
    args = ap.parse_args()

    if not args.skip_download:
        try:
            subprocess.run([sys.executable, str(ROOT / "scripts" / "fetch_nasa.py"), "--spacecraft", "both"], check=True)
        except subprocess.CalledProcessError as exc:
            raise SystemExit(f"NASA acquisition failed (exit code {exc.returncode}). See the acquisition diagnostic immediately above.")

    # Core Phase A intentionally uses the particle-only daily backbone.
    # Geometry and plasma speed are required only for the later propagation-shift extension.
    p10 = ROOT / "data" / "nasa" / "processed" / "p10_cpi_daily.csv"
    p11 = ROOT / "data" / "nasa" / "processed" / "p11_cpi_daily.csv"
    if not p10.exists() or not p11.exists():
        raise SystemExit("Processed CPI daily files are missing; run scripts/fetch_nasa.py first.")

    subprocess.run([
        sys.executable, str(ROOT / "scripts" / "replicate.py"),
        "--p10", str(p10), "--p11", str(p11),
        "--channel", args.channel,
        "--out", str(ROOT / "results" / "phase_a" / args.channel),
    ], check=True)


if __name__ == "__main__":
    main()
