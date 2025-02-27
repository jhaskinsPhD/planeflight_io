#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 13:25:51 2024

@author: u6044586
"""
import xarray as xr 
import sys 
import os 

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/osbpack_io.py')
from obspack_io import * 

 
ds0=xr.open_dataset('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/data/data/processed/obspack_ch4.20060607.nc')
ds1=xr.open_dataset('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jhask/GEOSChem/GC_RunDirs/v14.2.3_SOAS2013/Obspack_IO/Inputs/obspack_input.20130602.nc') 


has=['obs','calendar_components','time','time_components','latitude','longitude','altitude','obspack_id', 'CT_sampling_strategy']
for var in has: 
    if var in list(ds0.data_vars): 
        print('OBSPACK GOOD: (', var,')\n')
        print('Var:', ds0[var])
        print('Values:', ds0[var].values[0])
        print('Type of var:', type(ds0[var].values))
        print('Type of var values:', type(ds0[var].values[0]))
        print('') 
    if var in list(ds1.data_vars): 
        print('OBSPACK MINE: (', var,')\n')
        print('Var:', ds1[var])
        print('Values:', ds1[var].values[0])
        print('Type of var:', type(ds1[var].values))
        print('Type of var values:', type(ds1[var].values[0]))
        print('') 
        
              
# home='/uufs/chpc.utah.edu/common/home/haskins-group1/users/jhask/GEOSChem/GC_RunDirs/v14.2.3_SOAS2013/Obspack_IO/Inputs/'         
# filenames = write_obspack_inputs_ground('SOAS', 
#                                       lat=32.903281,  # lat of the sample site in Deg N
#                                       lon=-87.249942,  # lat of the sample site in Deg E
#                                       alt=139, # height above ground level to sample at in meters
#                                       datestart='20130601 00:00:00', # When I want to start sampling 
#                                       dateend='20130715 00:00:00',   # When I want to end sampling 
#                                       samplefreq=1800, # How frequently the model should sample  in s
#                                       sample_stragety=4,  # Set it to give me an avg hourly output when it samples.
#                                       outpath=home) 
