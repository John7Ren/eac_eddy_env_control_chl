"""Step 04: Build Paper 1 climatology dictionaries from prepared inputs.

This script connects the cleaned module functions to your existing workflow:
- background-map climatology for CHL, Zeu, MLD, N2, and nitracline depth;
- regional depth/environment climatologies;
- regional eddy-induced CHL-a anomaly climatologies.

Run from the repository root:
    PYTHONPATH=src/python python scripts/python/04_build_climatologies.py --config config/paths_template.yml
"""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from eac_eddy_chl.climatology import build_background_map_climatology, build_depth_climatology_for_region
from eac_eddy_chl.io import load_pickle, save_pickle


def read_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def main(config_path: str | Path) -> None:
    cfg = read_config(config_path)
    paths = cfg["paths"]
    out_dir = Path(paths["processed_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    input_bundle = load_pickle(out_dir / "step01_input_bundle.pkl")
    woa_bundle = load_pickle(out_dir / "step02_woa_background_bundle.pkl")

    ce_subset = input_bundle["ce_subset"]
    ae_subset = input_bundle["ae_subset"]
    region_masks = input_bundle["region_masks"]
    chl = input_bundle["chl"]
    k490 = input_bundle["k490"]
    mld = input_bundle["mld"]

    print("Building background-map climatology...")
    climt_bg_dict, climt2d_bg_dict = build_background_map_climatology(
        chl_ts=chl["chl_ts"],
        zeu_ts=k490["zeu_ts"],
        mld_holte=mld["mld_holte"],
        N_woa=woa_bundle["N_woa"],
        N2_woa=woa_bundle["N2_woa"],
        lat_woa=woa_bundle["lat_woa"],
        lon_woa=woa_bundle["lon_woa"],
        depth_N_woa=woa_bundle["depth_N_woa"],
        Lat_c=chl["Lat_c"],
        Lon_c=chl["Lon_c"],
        Lat_k=k490["Lat_k"],
        Lon_k=k490["Lon_k"],
    )

    print("Building regional climatologies...")
    regions_to_run = cfg.get("regions_to_run", ["ED", "Open", "Low-N"])
    depth_climatologies = {}
    eddy_climatologies = {}
    for region_name in regions_to_run:
        print(f"  {region_name}")
        depth_dict, eddy_dict = build_depth_climatology_for_region(
            region_name=region_name,
            region_masks=region_masks,
            ce_subset=ce_subset,
            ae_subset=ae_subset,
            chl_ts=chl["chl_ts"],
            zeu_ts=k490["zeu_ts"],
            N_woa=woa_bundle["N_woa"],
            N2_woa=woa_bundle["N2_woa"],
            mld_holte=mld["mld_holte"],
            time_woa=woa_bundle["time_woa"],
            depth_TS_woa=woa_bundle["depth_TS_woa"],
            depth_N_woa=woa_bundle["depth_N_woa"],
            lat_woa=woa_bundle["lat_woa"],
            lon_woa=woa_bundle["lon_woa"],
        )
        depth_climatologies[region_name] = depth_dict
        eddy_climatologies[region_name] = eddy_dict

    outputs = {
        "climt_bg_dict": climt_bg_dict,
        "climt2d_bg_dict": climt2d_bg_dict,
        "depth_climatologies": depth_climatologies,
        "eddy_climatologies": eddy_climatologies,
    }
    out_file = out_dir / "step04_climatology_outputs.pkl"
    save_pickle(outputs, out_file)
    print(f"Saved: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/paths_template.yml")
    args = parser.parse_args()
    main(args.config)
