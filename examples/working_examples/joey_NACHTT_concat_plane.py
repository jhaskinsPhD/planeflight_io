#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 14:33:58 2024

@author: u6044586
"""
import planeflight_io as pln

################################################################################
# Concatenate the plane log outputs
################################################################################

# # Path where output plane.log files ive:
# dirpath=('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
#          '/GC_RunDirs/gc_2x25_nacht2011_base/OutputDir/Plane_Logs/')

# # Path to your GEOS-Chem run's species database yaml file
# spdb_yaml=('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
#            '/GC_RunDirs/gc_2x25_nacht2011_base/species_database.yml')
# config_yaml=('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
#              '/GC_RunDirs/gc_2x25_nacht2011_base/geoschem_config.yml')


dirpath = ('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
           '/GC_RunDirs/gc_2x25_PLYACLNO2/OutputDir/Plane_Logs/')
config_yaml = ('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
               '/GC_RunDirs/gc_2x25_PLYACLNO2/geoschem_config.yml')
spdb_yaml = ('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem'
             '/GC_RunDirs/gc_2x25_PLYACLNO2/species_database.yml')


# Read in all plane.log files, concatenate them into a single X-Array netcdf file
# that will be stored at the dirpath which will include the attributes telling you
# what each variable is, and what units it is in
# and because we used tracer names in the plane input files, we need to set
# convert2_molmol = True so our output is in units of mol/mol.

model = pln.read_and_concat_planelogs(
    dirpath,
    spdb_yaml,
    config_yaml,
    as_xarr=False,
    convert2_molmol=True,
    overwrite=False,
    output_file='/uufs/chpc.utah.edu/common/home/u6044586/nachtt.xlsx')

################################################################################
# Doing model measurement comaprisons of GEOS-Chem at NACHTT Site
################################################################################
# Open concatenated plane.log data saved as xarray dataset at:
# model_pth=('/uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/GEOSChem/GC_RunDirs'
#            '/gc_2x25_nacht2011_base/OutputDir/Plane_Logs/planelog_concat_20110217_20110314.nc')
# model=xr.open_dataset(model_pth)

# # Open NACHTT Campaign Data:
# data_pth=('/uufs/chpc.utah.edu/common/home/haskins-group1/data/Campaign_Data'
#           '/Raw_Data/NACHTT_2011/data/elevator/NACHTT_2011_1min_Merged.nc')
# data=xr.open_dataset(data_pth)

# # Remove data that wasn't used to sample the model b/c of nans. Needed to
# # get data same len as model outputs:
# data=data.where((np.isnan(data['Lat'])==False) &
#               (np.isnan(data['Lon'])==False) &
#               (np.isnan(data['Tower_Height_m'])==False) &
#               (np.isnan(data['time'])==False), drop=True)

# # Scatter plot compairison to make sure units OK:
# plt.scatter(data['ClNO2_pptv'],model['ClNO2']*1e12,color='r', s=5)
