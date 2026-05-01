# Next steps for the cleaned Paper 1 workflow

1. Copy this package into the root of your GitHub repo.
2. Copy `config/paths_template.yml` to `config/paths_local.yml` and edit paths.
3. Install the local package:

```bash
python -m pip install -e .
```

4. Run scripts from the repo root:

```bash
python scripts/python/01_load_inputs.py --config config/paths_local.yml
python scripts/python/02_prepare_background_fields.py --config config/paths_local.yml
python scripts/python/03_build_region_masks.py --config config/paths_local.yml
python scripts/python/04_build_climatologies.py --config config/paths_local.yml
```

Large outputs should go to your external data folder and then Zenodo, not GitHub.
