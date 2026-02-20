# planeflight_io

Python tools to create and read GEOS-Chem planeflight diagnostic input/output files, with a focus on comparing model output to aircraft campaign data.

> **Compatibility:** GEOS-Chem v13.0.0 and later.

---

## Overview

`planeflight_io.py` provides functions to:
- **Create** `Planeflight.dat.YYYYMMDD` input files that tell GEOS-Chem when, where, and what to sample along a flight track.
- **Read** GEOS-Chem `plane.log` output files into pandas DataFrames or xarray Datasets, individually or concatenated across multiple days.

---

## Installation

### Option 1: PYTHONPATH (no pip required)

Clone the repo and add it to your `PYTHONPATH`:

```bash
git clone https://github.com/your-org/planeflight_io.git
export PYTHONPATH="/path/to/planeflight_io:$PYTHONPATH"
```

Add the `export` line to your `~/.bashrc` or `~/.bash_profile` to make it permanent. Then `import planeflight_io as pln` works from anywhere.

### Option 2: pip editable install

```bash
git clone https://github.com/your-org/planeflight_io.git
cd planeflight_io
pip install -e .
```

After this, `import planeflight_io as pln` works from any Python environment that uses that pip installation — no path manipulation needed.

---

## Quick Start

```python
import planeflight_io as pln

# --- Create input files ---
pln.make_planeflight_inputs(
    savedir='/path/to/output/',
    gc_config='/path/to/geoschem_config.yml',
    datetimes=senex_time,   # pd.Series of pd.Timestamps (UTC)
    lat_arr=senex_lat,
    lon_arr=senex_lon,
    vert_arr=senex_pres,    # Pressure array in hPa
    vert_is_pres=True,
    tracers=['NO', 'O3', 'CO'],
    diags=['NOy', 'RO2'],
)

# --- Read a single plane.log ---
df, info_dict = pln.read_planelog(
    '/path/to/plane.log.20170116',
    spdb_yaml='/path/to/species_database.yml',
    config_yaml='/path/to/geoschem_config.yml',
)

# --- Concatenate multiple plane.log files into one NetCDF ---
ds = pln.read_and_concat_planelogs(
    '/path/to/planelog_dir/',
    spdb_yaml='/path/to/species_database.yml',
    config_yaml='/path/to/geoschem_config.yml',
    as_xarr=True,
    output_dir='/path/to/output/',
    output_file='all_PlaneLogs',
)
```

---

## Examples

Interactive Jupyter Notebooks and plain Python scripts are provided in the `examples/` folder.

| File | Description |
|---|---|
| `examples/Examples- Creating PlaneFlight Input Files.ipynb` | Step-by-step notebook for creating input files using SENEX campaign data |
| `examples/Examples- Reading In PlaneFlight Output Files.ipynb` | Step-by-step notebook for reading and plotting plane.log output |
| `examples/Examples- Creating PlaneFlight Input Files.py` | Same content as the notebook, as a plain Python script |
| `examples/Examples- Reading In PlaneFlight Output Files.py` | Same content as the notebook, as a plain Python script |
| `examples/working_examples/` | Real-world scripts used for various atmospheric campaigns |

To run the notebooks or example scripts, set the `path_to_examples` variable near the top of each file to point to the `examples/` folder in your local clone. No other edits are needed to run them with the shipped example data.

---

## Dependencies

- Python >= 3.8
- [`numpy`](https://numpy.org/)
- [`pandas`](https://pandas.pydata.org/)
- [`xarray`](http://xarray.pydata.org/en/stable/)
- [`pyyaml`](https://pyyaml.org/)

---

## GEOS-Chem Compatibility

This package is designed for use with **GEOS-Chem v13.0.0 and later**. The planeflight diagnostic and `geoschem_config.yml` format it relies on were introduced in v13. Earlier versions are not supported.

For more information on the planeflight diagnostic, see the [GEOS-Chem documentation](https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/planeflight.html).
