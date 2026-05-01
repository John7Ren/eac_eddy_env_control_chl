"""Step 02: Load WOA23 and calculate density/N2 background fields.

Run from the repository root:
    PYTHONPATH=src/python python scripts/python/02_prepare_background_fields.py --config config/paths_template.yml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from eac_eddy_chl.io import load_woa23_monthly, save_pickle


def read_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main(config_path: str | Path) -> None:
    cfg = read_config(config_path)
    paths = cfg["paths"]
    out_dir = Path(paths["processed_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading WOA23 monthly T/S/NO3 and calculating density/N2...")
    woa = load_woa23_monthly(paths["woa_root"])
    out_file = out_dir / "step02_woa_background_bundle.pkl"
    save_pickle(woa, out_file)
    print(f"Saved: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/paths_template.yml")
    args = parser.parse_args()
    main(args.config)
