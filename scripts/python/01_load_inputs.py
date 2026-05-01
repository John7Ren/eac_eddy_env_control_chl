"""Step 01: Load eddy, MODIS, MLD, and WOA inputs and save a prepared input bundle.

Run from the repository root:
    PYTHONPATH=src/python python scripts/python/01_load_inputs.py --config config/paths_template.yml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from eac_eddy_chl.io import (
    load_eddy_pickles,
    load_mld_holte,
    load_modis_chl,
    load_modis_k490_zeu,
    load_region_mask,
    load_woa23_monthly,
    prepare_eddy_dicts,
    save_pickle,
)


def read_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main(config_path: str | Path) -> None:
    cfg = read_config(config_path)
    paths = cfg["paths"]
    out_dir = Path(paths["processed_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading CE/AE eddy pickles...")
    ce, ae = load_eddy_pickles(paths["data_root"])
    ce_subset, ae_subset = prepare_eddy_dicts(
        ce,
        ae,
        donut_file=paths.get("donut_file"),
        start_date=cfg.get("start_date", "2003-01-01"),
        end_date=cfg.get("end_date", "2022-12-31"),
    )

    print("Loading MODIS CHL...")
    chl = load_modis_chl(paths["chl_dir"])

    print("Loading MODIS K490 and calculating Zeu...")
    k490 = load_modis_k490_zeu(paths["k490_dir"])

    print("Loading MLD...")
    mld = load_mld_holte(paths["mld_file"])

    print("Loading region masks...")
    region_masks = load_region_mask(paths["region_mask_file"])

    bundle = {
        "ce_subset": ce_subset,
        "ae_subset": ae_subset,
        "chl": chl,
        "k490": k490,
        "mld": mld,
        "region_masks": region_masks,
    }
    out_file = out_dir / "step01_input_bundle.pkl"
    save_pickle(bundle, out_file)
    print(f"Saved: {out_file}")

    if cfg.get("load_woa_in_step01", False):
        print("Loading WOA23 and calculating density/N2. This can be slow...")
        woa = load_woa23_monthly(paths["woa_root"])
        save_pickle(woa, out_dir / "step01_woa_bundle.pkl")
        print(f"Saved: {out_dir / 'step01_woa_bundle.pkl'}")
    else:
        print("Skipped WOA23. Run 02_prepare_background_fields.py when ready.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/paths_template.yml")
    args = parser.parse_args()
    main(args.config)
