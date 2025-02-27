# -*- coding: utf-8 -*-
"""
Script with examples of reading in ObsPack & Planelight.dat output files
from GEOS-Chem. 

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

import planeflight_io as pln

# The path to the example folder under path_to_this(shouldn't need to update). 
path_to_examples= path_to_this+'\\examples'

# =============================================================================
# ############### Read in & Concatenate output plane.log files  #################
# =============================================================================  
# Path to folder containing example GEOS-Chem Planeflight Output 
filepath = path_to_examples+'/datafiles_for_examples/'

# =================================================================
#  OPTION 1: Read in a single plane.log files as a pandas dataframe, 
# =================================================================
pdat= pln.read_planelog(filepath+'planeflight_output/plane.log.20130610')

# =================================================================
#  OPTION 2: Read in several plane.log files, concatenate data from 
# into a single NetCDF file and take care of all the output if it came out in 
# molec cm-3 instead of mol mol-1 
# =================================================================

# Concatonate the output ObsPack netcdf files into a single netcdf file. 
# You only need to do once after you get GEOS-Chem output... 
# Set concat=True to create the concat'd data and to False to load it in other times. 
concat= False

if concat is True: # Concat all obspack files... 
    # Remove the old concated file if it exists... 
    if os.path.isfile(filepath+'all_PlaneLogs.nc4' ): os.remove(filepath+'all_PlaneLogs.nc4')
    
    # ~v13.3.4 Planeflight outputs are in molec cm-3 if you used tracer names (rather than tracer ID nums)
    # in your planeflight input files, so you may need to convert them  into mol mol-1.
    # This functionarliy is built into planeflight_io. If you DID user tracer#s 
    # in your planelog input file, don't use the # convert2molmol keyword! 
    fileconcat=pln.planelog_read_and_concat(filepath+'/planeflight_output/', 
                                            outdir=filepath,
                                            outfile='all_PlaneLogs',convert2_molmol=True) 
else:  # load in the concat'ed file. 
     fileconcat=filepath+'/all_PlaneLogs.nc4' 

# Open the concatenated obspack file with your GEOS-Chem output... 
ds=xr.open_dataset(fileconcat) 

# Plot the outputted O3 from plane.log
tracer='O3' 
print('Original Units of ',tracer,':  mol mol-1 * 1e9 = ppbv') 

plt.plot(ds.time, ds[tracer]*1e9, label='ObsPack', linewidth=2) 
plt.title(tracer+' (ppbv)')
plt.legend()
plt.tight_layout()
plt.show()
            
