"""Step 03: Check the region-mask file and print available region keys.

This script does not build the masks from scratch yet. It validates the saved
region mask dictionary that you already used in the paper workflow.

Run from the repository root:
    PYTHONPATH=src/python python scripts/python/03_build_region_masks.py --config config/paths_template.yml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml

from eac_eddy_chl.io import load_region_mask


def read_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main(config_path: str | Path) -> None:
    cfg = read_config(config_path)
    masks = load_region_mask(cfg["paths"]["region_mask_file"])
    print("Region-mask keys:")
    for key, mask in masks.items():
        arr = np.asarray(mask)
        print(f"  {key}: shape={arr.shape}, selected cells={np.sum(arr.astype(bool))}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/paths_template.yml")
    args = parser.parse_args()
    main(args.config)
