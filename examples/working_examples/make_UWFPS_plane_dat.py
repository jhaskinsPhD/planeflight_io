#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:12:05 2024

@author: u6044586
"""
import pandas as pd
import numpy as np
import xarray as xr
import planeflight_io as pln

# Set path to UWPFS merged 1 min campaign data:
uwfps_path = ('/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data'
              '/Raw_Data/UWFPS_2017/data/UWFPS_2017_1min_Merged.nc')

# Load in UWFPS 2017 data:
ds = xr.load_dataset(uwfps_path)

# Set path to the input file we want to get a list of tracers from!
inp_file = ('/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs'
            '/gc_2x25_merra2_fullchem_base/geoschem_config_run_jan_to_apr_2017.yml')

# Set path to the directory where Planedat output files will be saved:
savedir = ('/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data'
           '/GC_Input_Files/UWFPS_2017/PlaneDat_Inputs/2024_09_09-SHUYING/')

# Make sure pressure data is formatted correctly / should be in units of hPa.
print(ds.StaticPres.Units)
ds.StaticPres[ds.StaticPres < 0] = np.nan  # Don't allow negatives.

# Get time (in UTC) to not be an index and not be weird.
time0 = pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)

# Get required data arrays:
lat0 = ds.Lat.values
lon0 = ds.Lon.values
pres0 = ds.StaticPres.values

# Drop any points where any lat/lon/pressure/time is a NaN:
lat, lon, time, pres = pln.remove_nan_rows(lat0, lon0, time0, pres0)

# Make UWFPS planeflight input files:
pln.make_planeflight_inputs(savedir=savedir,
                            gc_config=inp_file,
                            datetimes=time,
                            lat_arr=lat,
                            lon_arr=lon,
                            vert_arr=pres,
                            vert_is_pres=True,
                            tracers='?ALL?',
                            diags='?ALL?',
                            username='jdh',
                            overwrite=True,
                            use_tracer_names=True)
