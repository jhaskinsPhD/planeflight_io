# -*- coding: utf-8 -*-
"""
Script with examples of reading in plane.log output files from GEOS-Chem.

@author: Dr. Jessica D. Haskins
"""

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt 
import os 
import sys 

# Update 'path_to_this' to point to wherever you have put 'planeflight_io' on your computer.
path_to_this= '/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/planeflight_io/'
sys.path.insert(0, path_to_this)

import planeflight_io as pln

# The path to the example folder under path_to_this (shouldn't need to update).
path_to_examples= path_to_this+'/examples'

# =============================================================================
# ############### Read in & Concatenate output plane.log files  #################
# =============================================================================
# Path to folder containing example GEOS-Chem Planeflight Output
filepath = path_to_examples+'/datafiles_for_examples/'

# Paths to the config files required by read_planelog() and read_and_concat_planelogs():
spdb_yaml = filepath + 'species_database.yml'
config_yaml = filepath + 'geoschem_config.yml'

# =================================================================
#  OPTION 1: Read in a single plane.log file as a pandas dataframe.
# =================================================================
# outputs_nums/ contains plane.log files made from input files that used tracer
# numbers (rather than tracer names), so concentrations are already in mol/mol.
df, info_dict = pln.read_planelog(filepath+'outputs_nums/plane.log.20170116',
                                  spdb_yaml=spdb_yaml,
                                  config_yaml=config_yaml)

# =================================================================
#  OPTION 2: Read in several plane.log files, concatenate data from
# all of them into a single NetCDF file. If outputs came from input
# files that used tracer names (rather than tracer numbers), set
# convert2_molmol=True to convert molec/cm3 output to mol/mol.
# =================================================================

# Concatenate the output plane.log files into a single file.
# You only need to do this once after you get GEOS-Chem output.
# Set concat=True to create the concat'd data and False to load it in other times.
concat= False

if concat is True:
    # Remove the old concatenated file if it exists...
    if os.path.isfile(filepath+'all_PlaneLogs.nc'): os.remove(filepath+'all_PlaneLogs.nc')

    # outputs_names/ contains plane.log files made from input files that used tracer
    # names, so concentrations are in molec/cm3 and need conversion to mol/mol.
    ds = pln.read_and_concat_planelogs(filepath+'outputs_names/',
                                       spdb_yaml=spdb_yaml,
                                       config_yaml=config_yaml,
                                       as_xarr=True,
                                       convert2_molmol=True,
                                       output_dir=filepath,
                                       output_file='all_PlaneLogs',
                                       overwrite=True)
else:  # Load the previously saved concatenated file.
    ds = xr.open_dataset(filepath+'all_PlaneLogs.nc')

# Plot the outputted O3 from plane.log
tracer='O3'
print('Original Units of ',tracer,':  mol mol-1 * 1e9 = ppbv')

plt.plot(ds.time_UTC, ds[tracer]*1e9, label='Model', linewidth=2)
plt.title(tracer+' (ppbv)')
plt.legend()
plt.tight_layout()
plt.show()
            
