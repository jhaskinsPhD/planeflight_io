# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 12:04:04 2022

@author: jhask
"""

# -*- coding: utf-8 -*-
"""
Script with examples of reading in ObsPack  output files from GEOS-Chem. 

@author: Dr. Jessica D. Haskins
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt 
import os 
import sys 

# Update 'path_to_this' to point to wherever you have put 'gcpy_campaigns' on your computer. 
path_to_this= 'C:\\Users\\jhask\\OneDrive\\Documents\\Python\\gcpy_extended\\gcpy_campaigns'
sys.path.insert(0, path_to_this) 

import obspack_io as obs

# The path to the example folder under path_to_this(shouldn't need to update). 
path_to_examples= path_to_this+'\\examples'

# =============================================================================
# ############### Read in & Concatenate Output ObsPack files  #################
# =============================================================================  
# Path to folder containing example GEOS-Chem ObsPack Output 
filepath = path_to_examples+'/datafiles_for_examples/'

# Concatonate the output ObsPack netcdf files into a single netcdf file. 
# You only need to do once after you get GEOS-Chem output... 
# Set concat=True to create the concat'd data and to False to load it in other times. 
concat= True

if concat is True: # Concat all obspack files... 
    # Remove the old concated file if it exists... 
    if os.path.isfile(filepath+'/obspack_output/all_Obspacks.nc' ): os.remove(filepath+'/all_Obspacks.nc')
    
    # NOTE: If you store the concated .nc file in the same place as the individuals, 
    # and then try to re-concat them without deleting that file first... you will get an error 
    # because it gets confused and tries to read that concate'd file and concat it with the indvs! 
    # So always store the concated file somewhere else or make sure you delete it before you regenerate it! 
    
    fileconcat = obs.read_and_concat_output_files(filepath+'/obspack_output/', outdir=filepath, outfile='all_ObsPacks')
    
else: # load in the concat'ed file. 
     fileconcat=filepath+'/all_ObsPacks.nc' 
    
# Open the concatenated obspack file with your GEOS-Chem output... 
ds=xr.open_dataset(fileconcat) 

# Plot the outputted O3 from ObsPack
tracer='O3' 
print('Original Units of ',tracer,': ', ds[tracer].units ,'* 1e9 = ppbv') 

plt.plot(ds.time, ds[tracer]*1e9, label='ObsPack', linewidth=2) 
plt.title(tracer+' (ppbv)')
plt.legend()
plt.tight_layout()
plt.show()
            