# planeflight_io
This repository contains Python 3 scripts used to prep input files and read output files for the GEOS-Chem planeflight diagnostic, with a focus on functionality related to comparing model output to campaign data.

Note: ** Functionality designed to be compatible with GEOS-Chem v13.0.0 and onwards **

## planeflight_io.py
Functions in planeflight_io.py allow you to create Planeflight.dat input files that are formatted correctly to be read in by GEOS-Chem. Functionality is included to make files with Pressure or Altitude inputs. Additionally, functions allow you to read in GEOS-Chem output plane.log files as pandas dataframes or xarray datasets, either individually or concatenated into a single file.

## Requirements:
* Python 3
* `numpy`: https://numpy.org/
* `pandas`: https://pandas.pydata.org/
* `xarray`: http://xarray.pydata.org/en/stable/
* `pyyaml`: https://pyyaml.org/


## Examples:
Examples are provided in the subfolder for making input files for various campaigns as well as examples of how to read in data using this library.
