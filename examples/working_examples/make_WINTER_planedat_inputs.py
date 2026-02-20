#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 11:26:51 2025

@author: u6044586
"""

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

# Set path to WINTER merged 1 min campaign data: 
winter_path='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Merged_Data/l3_final_merges/old/WINTER.nc'

# Load in WINTER 2015 data: 
ds=xr.load_dataset(winter_path) 

# # Get time (in UTC) to not be an index and be pandas._libs.tslibs.timestamps.Timestamp
time0= pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)

# Get Static pressure in hPa: 
pres0=ds.RAF_PSXC.values
    
# Get plane lat & plane long: 
lat0=ds.Lat.values
lon0=ds.Lon.values

# Drop any points where any lat/lon/pressure/time is a nan: 
lat,lon,time,pres= pln.remove_nan_rows(lat0,lon0,time0,pres0)

# Set path to the input file we want to get a list of tracers from! 
gc_config='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base/geoschem_config_run_jan_to_apr_2017.yml'

# Set path to the directory where Planedat output files will be saved:
savedir ='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/GC_Input_Files/WINTER_2015/PlaneDat_Inputs/'

# Make WINTER planeflight input files: 
pln.make_planeflight_inputs(savedir=savedir,                 # Where input files will be saved 
                            gc_config=gc_config,             # Path to your run's gc_config.yml file 
                            datetimes=time,                  # Times to sample model                
                            lat_arr=lat,                     # Lats to sample model
                            lon_arr=lon,                     # Lons to sample model 
                            vert_arr=pres,                   # Pressure/altitude array 
                            vert_is_pres=True,               # Tell it the vert coord is pressure!  
                            tracers='?ALL?',                 # Wildcard to sample all tracesr
                            diags='?ALL?',                   # Wildcard to sample all optional diagnosticss
                            username= 'jdh',                 # Who made this file 
                            overwrite=True,                  # Overwrite files if they already exist.
                            use_tracer_names=True)          # Use tracer #s so outputs in mol/mol   
