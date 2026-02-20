# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:40:31 2021

@author: jhask
"""
import xarray as xr
import pandas as pd
import os
import numpy as np
import planeflight_io as pln

path_to_examples = os.path.dirname(os.path.abspath(__file__))

# ==============================================================================
# Background: what is the planeflight diagnostic?
# ==============================================================================
# The GEOS-Chem planeflight diagnostic samples the model at arbitrary locations
# and times during a simulation — following an aircraft track, ship cruise, or
# any set of observation points you define. Sampled values are written to
# plane.log files (one per day) that you can then read back in and compare to
# your observations.
#
# You define WHERE, WHEN, and WHAT to sample in Planeflight.dat.YYYYMMDD input
# files. This script creates those files using pln.make_planeflight_inputs().
#
# Full documentation on the diagnostic:
# https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/planeflight.html
#
# ==============================================================================
# What do you actually need to make planeflight input files?
# ==============================================================================
# Less than you might think. The only truly required inputs are:
#
#   ALWAYS REQUIRED:
#     - Your flight/observation data: time, lat, lon, pressure (or altitude)
#     - An explicit list of species you want sampled (e.g. ['NO', 'O3', 'CO'])
#     - Your simulation type string (e.g. 'fullchem') — only if you also want
#       optional diagnostics beyond tracers
#
#   OPTIONAL — but unlocks more features:
#     - geoschem_config.yml (pass as gc_config=...) lets you:
#         * Use tracers='?ALL?' to sample every advected species automatically
#         * Get output in mol/mol (instead of molec/cm3) via automatic
#           name-to-number mapping
#         * Validate your tracer names against the actual species in your run
#
# The examples below start with the simplest possible call (no config files at
# all) and progressively add capabilities as you have more information available.

# ==============================================================================
# Step 1: Load your flight data
# ==============================================================================
# make_planeflight_inputs() requires four arrays describing the observation
# points. These conventions must be matched exactly:
#
#   Argument    Required type              Units / convention
#   ----------  -------------------------  -------------------------------------------
#   datetimes   pd.Series of pd.Timestamp  UTC (not local time)
#   lat_arr     1-D array-like             degrees North (-90 to 90)
#   lon_arr     1-D array-like             degrees East (-180 to 180, NOT 0-360)
#   vert_arr    1-D array-like             pressure in hPa  OR  altitude in meters
#
# All four must be the same length and NaN-free. The function does not silently
# re-project or re-scale inputs — bad units or conventions will produce silently
# wrong GEOS-Chem samples.
#
# Pressure vs. altitude (vert_is_pres):
#   For aircraft data, always use pressure (vert_is_pres=True). Altitude input
#   is ambiguous between "above ground" and "above sea level" conventions and is
#   only natively supported for CCGG/tower-type observations. See:
#   https://github.com/geoschem/geos-chem/issues/320

senex_pth = path_to_examples + '/datafiles_for_examples/SENEX.nc'
ds = xr.open_dataset(senex_pth)

# Use only the first 2 days so we generate just 2 input files:
unq_dates = np.unique(ds.time.dt.date)
ds = ds.where(((ds.time.dt.date == unq_dates[0]) | (ds.time.dt.date == unq_dates[1])), drop=True)

# Times — pd.Series of pd.Timestamps in UTC:
senex_time = pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)
print(f'Time:  {type(senex_time[0])}')  # must be pandas Timestamp

# Lat/lon — 1-D arrays, lon must be in range -180 to 180:
senex_lat = ds.GpsLat.values
senex_lon = ds.GpsLon.values
print(f'Lat:   {senex_lat.dtype}, range [{senex_lat.min():.2f}, {senex_lat.max():.2f}] deg')
print(f'Lon:   {senex_lon.dtype}, range [{senex_lon.min():.2f}, {senex_lon.max():.2f}] deg')

# Pressure in hPa — confirm units before passing:
senex_pres = ds.StaticPrs.values
print(f'Pres:  {senex_pres.dtype}, units = {ds.StaticPrs.attrs["Units"]} (must be hPa / mb)')

# %%
# ==============================================================================
# Example 1: The minimum viable call — no geoschem_config.yml needed
# ==============================================================================
# You don't need your GEOS-Chem run directory to create planeflight input files.
# This example shows the absolute minimum: an explicit list of tracers, no
# optional diagnostics, no config file.
#
# When is this useful?
#   - You haven't made a rundir yet and only need a few known tracers to be outputted. 
#
# One consequence: 
#    Without passing a gc_config.yml file, tracer names are written to 
#   the input file instead of tracer numbers. This causes GEOS-Chem to output advected
#   species concentrations in molec/cm3 rather than mol/mol. W
#          Why? See: https://github.com/geoschem/geos-chem/issues/796  
#   But, that's fine — you can convert at read time using convert2_molmol=True in
#   pln.read_and_concat_planelogs(). Example 3 shows how to avoid this entirely
#   if you do have a gc_config.yml available.

ex1_dir = path_to_examples + '/example1/'
if not os.path.isdir(ex1_dir):
    os.mkdir(ex1_dir)

pln.make_planeflight_inputs(
    savedir=ex1_dir,
    gc_config=None,              # No config file needed for a basic explicit-tracer call
    datetimes=senex_time,
    lat_arr=senex_lat,
    lon_arr=senex_lon,
    vert_arr=senex_pres,
    vert_is_pres=True,
    tracers=['NO', 'O3', 'CO'],  # Explicit list — required when gc_config=None
    username='me',
    overwrite=True,
)

# %%
# ==============================================================================
# Example 2: Adding optional diagnostics — still no geoschem_config.yml needed
# ==============================================================================
# Beyond advected tracers, the planeflight diagnostic can output a variety of
# optional diagnostics: meteorological fields, aerosol optical depths, chemical
# family concentrations, and more. These are requested via the `diags` argument.
#
# Crucially, optional diagnostics do NOT require gc_config — they're organised
# by simulation type, not by species list. You only need to tell the function
# what kind of simulation you're running (simtype='fullchem').
#
# Use pln.get_compatible_input_diags() to explore what's available. Valid
# collection names include: 
#    'aer_uptake', 'aodb', 'aodc', 'aq_aer', 'chem_fams', 'defaults', 
#    'gmao_ice', 'gmao_met', 'hg', 'htep', 'isor', 'tomas'.

# See every optional diagnostic available for a fullchem simulation:
all_diags = pln.get_compatible_input_diags(simtype='fullchem', display=True)

# Or filter by collection to get just the ones you care about.
met_and_fam = pln.get_compatible_input_diags(
    simtype='fullchem',
    these_collections=['gmao_met', 'chem_fams'],
    display=True,
)

ex2_dir = path_to_examples + '/example2/'
if not os.path.isdir(ex2_dir):
    os.mkdir(ex2_dir)

pln.make_planeflight_inputs(
    savedir=ex2_dir,
    gc_config=None,              # Still no config file needed
    datetimes=senex_time,
    lat_arr=senex_lat,
    lon_arr=senex_lon,
    vert_arr=senex_pres,
    vert_is_pres=True,
    tracers=['NO', 'O3', 'CO'],
    diags=['NOy', 'RO2'],        # Optional diagnostics — only need simtype, not gc_config
    simtype='fullchem',          # Required when requesting optional diagnostics
    username='me',
    overwrite=True,
)

# %%
# ==============================================================================
# Example 3: With geoschem_config.yml — mol/mol output and species validation
# ==============================================================================
# If you have your geoschem_config.yml, passing it as gc_config unlocks three
# things that aren't possible without it:
#
#   1. Tracer numbers instead of names — the function reads the geoschem_config.yml 
#      file to map each species name to its tracer number and writes numbers to the 
#      input file. GEOS-Chem then outputs advected species in mol/mol dry, which is
#      directly comparable to observations with no further conversion needed.
#      (Without gc_config, names are written and output is in molec/cm3.)
#
#   2. tracers='?ALL?' wildcard — The config contains the full list of advected
#      species in your run. Passing '?ALL?' samples all of them automatically,
#      so you don't need to pass an explicit list of every tracer you want to output.
#      See Example 4.
#
#   3. Species validation — the function checks every tracer you request against
#      the config's species list, catching typos before they silently produce
#      empty columns in your plane.log output.
#
# If you don't have the geoschem_config.yml available, stick with Examples 1 & 2 and pass
# convert2_molmol=True when reading your output back in.

gc_config = path_to_examples + '/datafiles_for_examples/geoschem_config.yml'

ex3_dir = path_to_examples + '/example3/'
if not os.path.isdir(ex3_dir):
    os.mkdir(ex3_dir)

pln.make_planeflight_inputs(
    savedir=ex3_dir,
    gc_config=gc_config,         # Unlocks mol/mol output, species validation, and '?ALL?'
    datetimes=senex_time,
    lat_arr=senex_lat,
    lon_arr=senex_lon,
    vert_arr=senex_pres,
    vert_is_pres=True,
    tracers=['NO', 'O3', 'CO'],
    diags=['NOy', 'RO2'],
    username='me',
    overwrite=True,
    use_tracer_names=False,      # Default: write tracer numbers → output in mol/mol
)
# Open the files from this example and Example 1 side by side — the tracer
# entries will look different (numbers vs. names), and when you run GEOS-Chem
# the concentrations will be in different units as a result.

# %%
# ==============================================================================
# Example 4: Wildcards with geoschem_config.yml — everything, minus exclusions
# ==============================================================================
# With gc_config available you can use '?ALL?' to request every advected species
# and/or every compatible optional diagnostic at once. This is the most
# comprehensive option — useful when you want full output for an exploratory run
# or when you don't want to maintain a long explicit list.
#
# Combining '?ALL?' with *_minus lists lets you say "everything except..." which
# is often more concise than a long include list when you only want to drop a
# handful of irrelevant variables. Here we exclude some halogen tracers (not
# relevant to our study) and aerosol/ice diagnostics (to keep file size down).
#
# Note: use_tracer_names=True here so you can open these files alongside those
# from Example 3 and directly compare the tracer-number vs. tracer-name format.

tracers_minus = ['ClNO2', 'Cl2', 'ClO', 'HOCl', 'HCl', 'BrCl']  # exclude halogen tracers

diags_minus = [
    "AODC_SULF", "AODC_BLKC", "AODC_ORGC", "AODC_SALA", "AODC_SALC",
    "AODC_DUST", "AODB_SULF", "AODB_BLKC", "AODB_ORGC", "AODB_SALA", "AODB_SALC",
    "AODB_DUST", "GMAO_ICE00", "GMAO_ICE10", "GMAO_ICE20",
    "GMAO_ICE30", "GMAO_ICE40", "GMAO_ICE50", "GMAO_ICE60", "GMAO_ICE70",
    "GMAO_ICE80", "GMAO_ICE90",
]

ex4_dir = path_to_examples + '/example4/'
if not os.path.isdir(ex4_dir):
    os.mkdir(ex4_dir)

pln.make_planeflight_inputs(
    savedir=ex4_dir,
    gc_config=gc_config,
    datetimes=senex_time,
    lat_arr=senex_lat,
    lon_arr=senex_lon,
    vert_arr=senex_pres,
    vert_is_pres=True,
    tracers='?ALL?',              # Sample every advected species in the run
    tracers_minus=tracers_minus,  # ...except these
    diags='?ALL?',                # Sample every compatible optional diagnostic
    diags_minus=diags_minus,      # ...except these
    username='me',
    overwrite=True,
    use_tracer_names=True,        # Write names → output in molec/cm3 (compare to Example 3!)
)
