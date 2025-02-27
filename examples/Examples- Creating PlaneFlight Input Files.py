# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:40:31 2021

@author: jhask
"""
import xarray as xr
import pandas as pd
import matplotlib as plt 
import os 
import sys 

# Update 'path_to_this' to point to wherever you have put 'gcpy_campaigns' on your computer. 
path_to_this= 'C:\\Users\\jhask\\OneDrive\\Documents\\Python\\gcpy_extended\\gcpy_campaigns'
sys.path.insert(0, path_to_this)

import obspack_io as obs
import planeflight_io as pln

# ======   Open the SENEX data, and make sure its formatted OK! ================
path_to_examples= path_to_this+'\\examples'

# Open SENEX merged data. 
senex=path_to_examples+'\\datafiles_for_examples\\SENEX.nc'
sn=xr.open_dataset(senex) 

# Convert times to a ** pandas datetime array ** 
time=pd.to_datetime(sn.time.values).to_series().reset_index(drop=True)
 

# =============================================================================
#  Example # 1:  Make a planeflight.dat file for all SENEX flights. Sample just NO, O3, CO.
# =============================================================================

tf= os.path.isdir(path_to_examples+'\\example1\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example1\\') # Make directory to hold output.
    
pln.make_planeflightdat_files(outpath=path_to_examples+'\\example1\\',
                            datetimes=time,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            typestr='aSENEX', # Don't make typestrs begin with "S" or errors! 
                            tracers=['NO', 'O3', 'CO'],
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'])
# =============================================================================
#  Example #2:  Make a planeflight.dat file for all SENEX flights. 
#              Sample all tracers in your input file and all optional diagnostics
# =============================================================================

tf= os.path.isdir(path_to_examples+'\\example2\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example2\\') # Make directory to hold output.
pln.make_planeflightdat_files(outpath=path_to_examples+'\\example2\\',
                            input_file=path_to_examples+'\\datafiles_for_examples\\input.geos.standard', 
                            typestr='aSENEX', # Don't make typestrs begin with "S" or errors! 
                            datetimes=time,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'],
                            print_diag_options= False)

# =============================================================================
# Example #3:  Make a planeflight.dat file for all SENEX flights. 
#               Sample all tracers in your input file except those in "minus". 
# =============================================================================

tf= os.path.isdir(path_to_examples+'\\example3\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example3\\') # Make directory to hold output.

minus=["AODC_SULF","AODC_BLKC","AODC_ORGC","AODC_SALA","AODC_SALC",
"AODC_DUST","AODB_SULF","AODB_BLKC","AODB_ORGC","AODB_SALA","AODB_SALC",
"AODB_DUST","HG2_FRACG" ,"HG2_FRACP","GMAO_ICE00","GMAO_ICE10","GMAO_ICE20", 
"GMAO_ICE30", "GMAO_ICE40","GMAO_ICE50", "GMAO_ICE60","GMAO_ICE70", "GMAO_ICE80","GMAO_ICE90"]

pln.make_planeflightdat_files(outpath=path_to_examples+'\\example3\\',
                            input_file=path_to_examples+'\\datafiles_for_examples\\input.geos.standard', 
                            typestr='aSENEX',
                            datetimes=time,
                            lat_arr=sn.GpsLat.values, 
                            lon_arr=sn.GpsLon.values, 
                            pres_arr=sn.StaticPrs.values,
                            username= 'jdh',
                            drop_dupes = False,
                            diags = ['all'],
                            diags_minus=minus,
                            print_diag_options= False)
    
