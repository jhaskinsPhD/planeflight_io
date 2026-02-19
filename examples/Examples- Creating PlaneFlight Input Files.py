# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:40:31 2021

@author: jhask
"""
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import os 
import sys 
import numpy as np

# Update 'path_to_this' to point to wherever you have put 'planeflight_io' on your computer.
path_to_this= '/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/planeflight_io/'
sys.path.insert(0, path_to_this)
import planeflight_io as pln

# Define path to examples folder: 
path_to_examples= path_to_this+'/examples'

#==============================================================================
# Step #0:  Visit Read-The-Docs to lear more about the Planeflight diagnostic: 
#==============================================================================
# https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/planeflight.html

# In this example, we'll be using the function pln.make_planeflight_inputs() to create input 
# plane.dat files for GEOS-Chem here. Read about the function's inputs/outputs here: 
"""
make_planeflight_inputs(savedir: str,
                        gc_config:str,
                        datetimes,
                        lat_arr,
                        lon_arr,
                        vert_arr,
                        vert_is_pres:bool,
                        tracers,
                        tracers_minus: list =[],
                        diags=[],
                        diags_minus: list =[],
                        username: str = 'user',
                        overwrite: bool = False,
                        use_tracer_names:bool=False)

INPUTS:
-------
Note: For all inputs where "ARRAY" is accepted, vars can be 1-D lists, np.ndarray,
      pd.Series, or xr.dataarray. Rounding/ # of accepted decimals and input data types are
      all automatically converted to what GEOS-Chem requires.

    (1) savedir      - STRING containing path to directory in which to save
                      the output planeflight.dat files at.

    (2) gc_config   - STRING containing path to your geoschem_config.yml file
                      used to read in/define valid names of transported tracers
                      and to determine simulation type, used to determine all valid
                      additional / optional planeflight diagnostics.

    (3) datetimes    - PANDAS SERIES of datetimes in UTC stored as Timestamps at which
                      to sample the model.Type(datetimes[0]) should return:
                      <class 'pandas._libs.tslibs.timestamps.Timestamp'>

                      To create input in the correct format do:
                            date_range = pd.date_range(start='2017-01-01', end='2017-01-03', freq='60s')
                            datetimes=pd.Series(date_range)

    (4) lat_arr      - ARRAY of latitudes at which to sample the model. (range: -90 to 90 deg)

    (5) lon_arr      - ARRAY of longitudes at which to sample the model (range: -180 to 180 deg)

    (6) vert_arr     - EITHER an ARRAY of pressures (hPa) OR altitudes above
                      the ground (meters) at which to sample the model. See
                      https://github.com/geoschem/geos-chem/issues/320
                      for discussion on whether input altitudes should be
                      "above ground" or "above sea level".

    (7) vert_is_pres - BOOLEAN indicating if "vert_arr" containined pressures or not.
                      When TRUE,  values are assumed to be pressures (hPa).
                      When FALSE, values are assumed to be altitudes (meters).

    (8) tracers     -  Either (1) an ARRAY of specific advected tracers you want to sample
                          OR (2) a STRING equal to '?ALL?' to sample all advected species
                      listed in your geoschem_config.yml file.

    (9) tracers_minus - (OPTIONAL) ARRAY containing STRINGS with all advected species
                      tracers you don't want to include (only relevant if you passed
                      tracers='?ALL?').

    (10) diags       -  (OPTIONAL) Either (1) an ARRAY containing STRINGS of any additional
                      diagnostics to sample from model (beyond tracers) OR (2) a
                      STRING equal to '?ALL?' to sample all available additional
                      diagnostics compatiable with your simulation type. Default is
                      to include the grid-box indexes planflight pulled from &
                      Pres/Temp/RH at center of grid box (e.g. ['GMAO_IIEV',
                      'GMAO_JJEV', 'GMAO_LLEV', 'GMAO_PRES','GMAO_RELH', 'GMAO_TEMP']).

    (11) diags_minus - (OPTIONAL) ARRAY containing STRINGS with all additional
                      diagnostics you don't want to include (only relevant
                      if you passed diags='?ALL?').

    (12) username   - (OPTIONAL) STRING contaiing name of user who created files.
                      This gets listed in header of resulting planedat input files.

    (13) overwrite  - (OPTIONAL) BOOLEAN of whether to overwrite any existing files
                      at 'savedir' with this name or not. If FALSE, & any files
                      under 'savedir' would be overwritten, a new sub-directory
                      under 'savedir' named "NEW_YYYYMMDD_HHMMSS" is created to
                      hold all the new output files. If TRUE, only files
                      with conflicting names under 'savedir' are overwritten.
                      Default is set to FALSE (not to overwrite files).

    (14) use_tracer_names -(OPTIONAL) BOOL indicating if you want to write the file
                      with tracer names rather than with tracer numbers. Default is
                      set to FALSE (to write file using tracer numbers) since outputs
                      like this have units of mol/mol dry. If set to true, input
                      files will be written with tracer names instead, which will
                      result in the output files have units of molec/cm3. While
                      using tracer names is more readable, it adds an extra step of
                      potential error to compare directly to observations. Thus,
                      writing files with tracer numbers is reccommended. This code
                      can accomodate either though when reading in output files.

OUTPUT:
------
    (1) One file for each day listed in 'datetimes' named "Planeflight.dat.YYYYMMDD"
        written to "savedir" that can be passed to GEOS-Chem as input files for
        the planeflight diagnostic. By default, outputted files convert tracer names
        to tracer numbers so that output files have concentrations in units of mol/mol
        rather than molec/cm3 (which occurs if tracer names are used instead).
"""

#==============================================================================
# Step #1:  Open your flight data and curate the inputs for make_planeflight_inputs(): 
#==============================================================================
# This example ships with some merged data from the SENEX campaign: 
senex_pth= path_to_examples+'/datafiles_for_examples/SENEX.nc'
ds=xr.open_dataset(senex_pth)

# Select only the first 2 days so we just generate 2 input files... 
unq_dates=np.unique(ds.time.dt.date)
ds = ds.where(((ds.time.dt.date == unq_dates[0]) | (ds.time.dt.date == unq_dates[1])), drop=True)

# 1.1 Curate array to pass as (required) arg for "datetimes"  
#    Pull out time (in UTC). GEOS-Chem needs to know when to sample the model. 
#    Then, convert it to a ** PANDAS SERIES of datetimes stored as Timestamps** 
senex_time=pd.to_datetime(ds.time.values).to_series().reset_index(drop=True)

# Check that the type of "time[0]" is <class 'pandas._libs.tslibs.timestamps.Timestamp'>
print(type(senex_time[0]))

# 1.2 Curate arrays to pass as (required) arg for "lat_arr" and "lon_arr" 
#     Pull out latitude & longitude. GEOS-Chem needs to know where the plane is. 
senex_lat= ds.GpsLat.values
senex_lon=ds.GpsLon.values

# Check that the type of these variables is a 1-D  array of one of the following types: 
#       a 1-D list, np.ndarray, pd.Series, or xr.dataarray
print(type(senex_lat), type(senex_lon))

# 1.3 Curate Pressure/Altitude data. 
#    Next you need to figure out if you're going to pass GEOS-Chem pressure or 
#    altitude values to determine what vertical level to sample the model at. 
#    It's generally more accurate if you use pressure inputs rather than altitude. 
#    You'll wnt to use whats typically called the "static presure" field measured
#    during aircraft campaigns. The SENEX data does indeed have this field! 

senex_pres=ds.StaticPrs 

# Check that the type of these variables is a 1-D  array of one of the following types: 
#       a 1-D list, np.ndarray, pd.Series, or xr.dataarray
print(type(senex_pres))

# Then if using pressure make sure its in hPa / mbar.  For alt(meters!) 
print(ds.StaticPrs.attrs['Units'])

# NOTE: Planeflight only supports input of altitude for CCGG or tower data.
# All other observations must specify pressure in hPa. The "TYPE" string is 
# automatically set in this code to be for "tower or ccgg data" if you input 
# altitudes rather than pressures. So, this code will allow you to "hack" your 
# way around that limitation if you really want to... though not advisable. 
# So, use pressure inputs whenever possible if trying to sample aircraft! 

#==============================================================================
# Step #2:  Figure out what vars you want GEOS-Chem to output. 
#==============================================================================
# GEOS-Chem can sample any advected species or a variety of optional additional 
# diagnostics within planeflight. To see all additional diagnostics compatible w
# with your simulation you can call this function to get a list of all optional diagnostic 
# "collections" that are compatible with your simulation type: 
diags= pln.get_compatible_input_diags(simtype='fullchem',display=True)

# To retrieve only specific collections of diagnostics rather than all compatible ones,
# pass a list of collection names to the 'these_collections' argument. Valid collection
# names include: 'aer_uptake', 'aodb', 'aodc', 'aq_aer', 'chem_fams', 'defaults',
# 'gmao_ice', 'gmao_met', 'hg', 'htep', 'isor', 'tomas'.
met_diags = pln.get_compatible_input_diags(simtype='fullchem',
                                           these_collections=['gmao_met','chem_fams'],
                                           display=True)

# The following examples show a few different ways to create the files using
# either a list of species to output,  a wildcard to sample all advected 
# species in your simulation, and how to select various optional diagnostics too! 

# The script uses your gc_config file to determine your simtype, get a list of 
# all advected species in your simulation to ensure what optional diagnositcs 
# you request are compatible with your simulation type and to allow you 
# to request all species with a wildcard. So, you must tell it where your config file
# is to use make_planeflight_inputs(). 
gc_config=path_to_examples+'/datafiles_for_examples/geoschem_config.yml'

# %%
# =============================================================================
#  Example # 1:  Make a planeflight.dat file for all SENEX flights. 
                 # Sample just NO, O3, CO fro advected species and 
                 # The optional diagnostics 'NOy', and 'RO2'
# =============================================================================
# Define & make directory to hold example 1 output if it doesn't exist: 
ex1_dir= path_to_examples+'/example1/'
if not os.path.isdir(ex1_dir): os.mkdir(ex1_dir) 
    
pln.make_planeflight_inputs(savedir=path_to_examples+'/example1/', # Where input files will be saved 
                            gc_config=gc_config,                   # Path to your run's gc_config.yml file 
                            datetimes=senex_time,                  # Times to sample model                
                            lat_arr=senex_lat,                     # Lats to sample model
                            lon_arr=senex_lon,                     # Lons to sample model 
                            vert_arr=senex_pres,                   # Pressure/altitude array 
                            vert_is_pres=True,                     # Tell it the vert coord is pressure!  
                            tracers=['NO', 'O3', 'CO'],  # List of tracers to sample
                            diags=['NOy','RO2'],         # List of optional diagnostics to sample 
                            username= 'me',                        # Who made this file 
                            overwrite=True,                        # Overwrite files if they already exist.
                            use_tracer_names=False)                # Use tracer #s so outputs in mol/mol    

# Try re-running this example & compare the input files generated 
# after setting "overwrite" to False (so you retain each example file): 
#    
#      - with use_tracer_names set to True 
#          Note: Except for the names in the list of things to smaple, the input 
#                files will be similar, but these two files would generate advected 
#                species conc ouputs from GEOS-Chem in very dif units! 
#      - without passing the diags keyword 
#      - by passing a different username

# %%
# =============================================================================
#  Example #2:  Make a planeflight.dat file for all SENEX flights. 
#              Sample all tracers in your input file and all optional diagnostics
# =============================================================================
# Define & make directory to hold example 2 output if it doesn't exist: 
ex2_dir= path_to_examples+'/example2/'
if not os.path.isdir(ex2_dir): os.mkdir(ex2_dir) 

pln.make_planeflight_inputs(savedir=path_to_examples+'/example2/', # Where input files will be saved 
                            gc_config=gc_config,                   # Path to your run's gc_config.yml file 
                            datetimes=senex_time,                  # Times to sample model                
                            lat_arr=senex_lat,                     # Lats to sample model
                            lon_arr=senex_lon,                     # Lons to sample model 
                            vert_arr=senex_pres,                   # Pressure/altitude array 
                            vert_is_pres=True,                     # Tell it the vert coord is pressure!  
                            tracers='?ALL?',     # Use wildcard' to request all advected species
                            diags='?ALL?',       # Use wildcard' to request all compatiable optional diagnostics
                            username= 'me',                        # Who made this file 
                            overwrite=True,                        # Overwrite files if they already exist.
                            use_tracer_names=False)                # Use tracer #s so outputs in mol/mol 

# Compare these generated input files to those you generated in the previous example. 

# %%
# # =============================================================================
# # Example #3:  Make a planeflight.dat file for all SENEX flights. 
# #              Sample all tracers in your input file except those in tracers_minus
#                and all diagnostics except those in diags_minus 
# # =============================================================================
# Define & make directory to hold example 3 output if it doesn't exist: 
ex3_dir= path_to_examples+'/example3/'
if not os.path.isdir(ex3_dir): os.mkdir(ex3_dir) 

# List of tracers to NOT include... 
tracers_minus=['ClNO2','Cl2','ClO','HOCl','HCl','BrCl'] 

# List of diagnostics to NOT include... 
diags_minus=["AODC_SULF","AODC_BLKC","AODC_ORGC","AODC_SALA","AODC_SALC",
"AODC_DUST","AODB_SULF","AODB_BLKC","AODB_ORGC","AODB_SALA","AODB_SALC",
"AODB_DUST","GMAO_ICE00","GMAO_ICE10","GMAO_ICE20", 
"GMAO_ICE30", "GMAO_ICE40","GMAO_ICE50", "GMAO_ICE60","GMAO_ICE70", "GMAO_ICE80","GMAO_ICE90"]

pln.make_planeflight_inputs(savedir=path_to_examples+'/example3/', # Where input files will be saved 
                            gc_config=gc_config,                   # Path to your run's gc_config.yml file 
                            datetimes=senex_time,                  # Times to sample model                
                            lat_arr=senex_lat,                     # Lats to sample model
                            lon_arr=senex_lon,                     # Lons to sample model 
                            vert_arr=senex_pres,                   # Pressure/altitude array 
                            vert_is_pres=True,                     # Tell it the vert coord is pressure!  
                            tracers='?ALL?',                       # Use wildcard' to request all advected species
                            tracers_minus=tracers_minus,     # Exclude tracers in this list though... 
                            diags='?ALL?',                         # Use wildcard' to request all compatiable optional diagnostics
                            diags_minus=diags_minus,         # Exclude these diagnostics though.... 
                            username= 'me',                        # Who made this file 
                            overwrite=True,                        # Overwrite files if they already exist.
                            use_tracer_names=True)           # Output files will have concentrations of Advected species in molec/cm3 

