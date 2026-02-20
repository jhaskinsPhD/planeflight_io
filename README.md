# planeflight_io

Python tools to create and read GEOS-Chem planeflight diagnostic input/output files, with a focus on comparing model output to aircraft campaign data.

> **Compatibility:** GEOS-Chem v13.0.0 and later.

---

## Overview

`planeflight_io.py` provides functions to:
- **Create** `Planeflight.dat.YYYYMMDD` input files that tell GEOS-Chem when, where, and what to sample along a flight track.
- **Read** GEOS-Chem `plane.log` output files into pandas DataFrames or xarray Datasets, individually or concatenated across multiple days.

---

## What You Need

### Your campaign / observation data

`make_planeflight_inputs()` requires four arrays describing the observation points. Match these types, units, and conventions exactly — the function does not silently re-project or re-scale inputs:

| Argument | Required type | Units / convention |
|---|---|---|
| `datetimes` | `pd.Series` of `pd.Timestamp` | UTC (not local time) |
| `lat_arr` | 1-D `np.ndarray` | degrees North (−90 to 90) |
| `lon_arr` | 1-D `np.ndarray` | degrees East (−180 to 180, **not** 0–360) |
| `vert_arr` | 1-D `np.ndarray` | pressure in **hPa** (preferred) or altitude in **meters above ground** |

All four arrays must be the same length and NaN-free.

**`vert_is_pres` — pressure vs. altitude:** The planeflight diagnostic only natively supports altitude input for CCGG/tower-type observations; all aircraft data should use pressure (`vert_is_pres=True`). Using altitude for aircraft is technically possible — `planeflight_io` sets the TYPE string accordingly — but it is not advisable because of ambiguity between "above ground" and "above sea level" conventions. See [GH issue #320](https://github.com/geoschem/geos-chem/issues/320) for the full discussion.

### Your GEOS-Chem run files

Two YAML files from your GEOS-Chem run directory are required:

- **`geoschem_config.yml`** — used by `make_planeflight_inputs()` to:
  1. Determine your simulation type, so only compatible optional diagnostics are offered.
  2. Retrieve the full list of advected species, so `tracers='?ALL?'` works correctly.
  3. Map species names → tracer numbers, so outputs come in `mol/mol` (see below).
- **`species_database.yml`** — used when *reading* output (by `read_planelog()` and `read_and_concat_planelogs()`) to attach metadata — molecular weight, chemical formula, long name, and units — to each output column.

Both files live in your GEOS-Chem run directory (the directory where you ran or will run GEOS-Chem).

---

## Tracer Numbers vs. Tracer Names

The `use_tracer_names` argument in `make_planeflight_inputs()` controls the units of advected species in the `plane.log` output:

| `use_tracer_names` | What goes in the input file | Units in `plane.log` output |
|---|---|---|
| `False` (default, **recommended**) | Tracer numbers (e.g. `1`, `2`, `3`) | **mol/mol dry** — directly comparable to observations |
| `True` | Tracer names (e.g. `'NO'`, `'O3'`) | **molec/cm³** — requires conversion before comparing |

Using tracer **numbers** produces output in `mol/mol dry`, which can be compared directly to aircraft observations without any further conversion. Using tracer **names** is more human-readable in the input file, but GEOS-Chem will output advected species concentrations in `molec/cm³`. See [GH issue #796](https://github.com/geoschem/geos-chem/issues/796) for the full explanation of why this unit difference exists.

**Recommendation:** leave `use_tracer_names=False` (the default) unless you have a specific reason to use names. If you do use tracer names, pass `convert2_molmol=True` when reading output with `pln.read_and_concat_planelogs()`.

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
    savedir='/path/to/output/',               # Directory to save Planeflight.dat.YYYYMMDD files
    gc_config='/path/to/geoschem_config.yml', # Your run's geoschem_config.yml (see "What You Need")
    datetimes=senex_time,                     # pd.Series of pd.Timestamps — must be UTC
    lat_arr=senex_lat,                        # Latitudes in degrees North (−90 to 90)
    lon_arr=senex_lon,                        # Longitudes in degrees East (−180 to 180, not 0–360)
    vert_arr=senex_pres,                      # Pressure in hPa (preferred) or altitude in m
    vert_is_pres=True,                        # True = pressure input; False = altitude input
    tracers=['NO', 'O3', 'CO'],               # Specific tracers, or '?ALL?' for all advected species
    diags=['NOy', 'RO2'],                     # Optional diagnostics, or '?ALL?' for all compatible
)

# --- Read a single plane.log ---
df, info_dict = pln.read_planelog(
    '/path/to/plane.log.20170116',
    spdb_yaml='/path/to/species_database.yml',   # For column metadata (MW, formula, units)
    config_yaml='/path/to/geoschem_config.yml',  # To identify tracers vs. diagnostics
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
