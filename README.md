# gcpy_campaigns
This repository contains Python 3 scripts used to prep inputs files and read output files for GEOS-Chem with a particular focus on functionality related to comparing model output to campaign data.

Note: ** Functionality designed to be compatible with GEOS-Chem v13.0.0 and onwards **

## obspack_io.py
Functions in osbpack_io.py allow you to create ObsPack input files that are formatted correctly to be read in by GEOS-Chem. Functionality is included to make files for sampling at ground sites or via a plane (or other mobile platform). Also contains functions that read in GEOS-Chem Obspack output files, with functionality to read in individual files, and concatenate them into a single netcdf file, also transforming output time into a python friendly datetime. 

## planeflight_io.py
Similarly, functions in planeflight_io.py allow you to create Planeflight.dat input files that are formatted correctly to be read in by GEOS-Chem. Functionality is included to make files with Pressure or Altitude inputs. Additionally, functions allow you to read in GEOS-Chem output plane.log files as pandas dataframes either individually or concatenated into a single file. 

## Requirements:
***(Universal dependencies:)*** 
* Python 3
* `numpy`: https://numpy.org/
* `pandas`: https://pandas.pydata.org/

***(For planeflight_io only)***
* `astropy.io` : https://docs.astropy.org/en/stable/index.html

***(For obspack_io only)***
* `xarray`: http://xarray.pydata.org/en/stable/


## Examples:
Examples are provided in the subfolder for making input files for various campaigns as well as examples of how to read in data using this library. 
