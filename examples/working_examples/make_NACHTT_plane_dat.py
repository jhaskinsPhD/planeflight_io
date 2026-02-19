#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 19:02:53 2024

@author: u6044586
"""
import sys 
import pandas as pd 
import numpy as np 
import xarray as xr 

# Add Jessica's library to make planeflight input files:
sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/planeflight_io/')
import planeflight_io as pln

# Set path to NAHCTT merged 1 min campaign data: 
nachtt_path='/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/Raw_Data/NACHTT_2011/data/elevator/NACHTT_2011_1min_Merged.nc'

# Load in NACHTT data: 
ds=xr.load_dataset(nachtt_path) 

# Set path to the input file we want to get a list of tracers from! 
inp_file='/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem/GC_RunDirs/gc_2x25_merra2_fullchem/geoschem_config_run1.yml'

# Set path to the directory where Planedat output files will be saved:
savedir = '/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data/GC_Input_Files/NACHTT_2011/PlaneDat_Inputs/2024_09_19-JOEY/'

# Get time (in UTC) to not be an index and not be weird. 
time= pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)

# Create smaller xarray with ONLY the data we need, drop any NaNs from it! 
nachtt=  xr.Dataset({ 'lat':xr.DataArray( data= ds.Lat.values , dims=['obs']),
                     'lon':xr.DataArray( data= ds.Lon.values , dims=['obs']),
                    'alt':xr.DataArray( data= ds.Tower_Height_m.values , dims=['obs']),
                    'time':xr.DataArray( data= time, dims=['obs'])}).dropna(dim='obs')

# Drop all vars if any is nan: 
nachtt= nachtt.dropna(dim='obs', how='any')

# Need to re-do time config after dropping nan stuff: 
time= pd.to_datetime(nachtt.time.values).to_series().reset_index(drop=True)

# Make NACHTT planeflight input files (using altitude since this is tower data):
pln.make_planeflight_inputs(savedir=savedir,
                            gc_config=inp_file,
                            datetimes=time,
                            lat_arr=nachtt.lat.values,
                            lon_arr=nachtt.lon.values,
                            vert_arr=nachtt.alt.values,  # in m
                            vert_is_pres=False,          # altitude, not pressure
                            tracers='?ALL?',
                            diags='?ALL?',
                            username='jdh',
                            overwrite=True,
                            use_tracer_names=True)
