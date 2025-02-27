"""
Scripts to make readable Obspack input files for GEOS-Chem for planes or ground sites
and read output from GEOS-Chem, concatonate files into a single .nc file for easy use.

Created on Sun Feb 14 12:05:25 2021.

@author: Dr. Jessica D.Haskins
"""
import pandas as pd
import numpy as np
import xarray as xr
import os 
import datetime
import glob
import sys

import campaign_utils as ut
from shutil import copyfile

#############################################################################
########################   Input Handling    ###############################
#############################################################################

def write_obspack_inputs_flights(sitename, lat , lon , alt,
                               time, sample_stragety: int, outpath: str, 
                               quiet:bool = False):
    """Create ObsPack input files for GEOS Chem v13.0.0. for a Flight.
    
    # =========================================================================
    #                             INPUTS
    # =========================================================================
    # sitename = string with sitename used in obspack ID, (e.g.'SOAS')
    # lat      = np array with latitudes of plane
    # lon      = np array with longitudes of plane (** MUST be between -180E, 180W!)
    # alt      = np array with plane height in meters above sea level
    # time     = pandas datetime series in UTC with times you want samples at.
    #
    # Sample Stragety = Integer specifiing how the model will average to
    #                   the timebase given (calc'd by  datestart,dateend,
    #                   and freq). Only Valid option are:
    #                             1  =   4-hour avg;
    #                             2  =   1-hour avg;
    #                             3  =   90-min avg;
    #                             4  =   instantaneous
    #
    # outfile_path  = string containing absolute path where the netcdf files
    #                 will be written.
    #
    # quiet = boolean of whether you'd like to print progress to screen.
    #
    # =========================================================================
    #                             OUTPUTS
    # =========================================================================
    #
    #  out= Returns list of netcdf files written to outpath for each day in
    #       the time period that  can be used as GEOS Chem inputs.
    #
    #  NOTE: The output files MUST be named with YYYYMMDD in order to be
    #        read in by GEOS-Chem properly. So do NOT change the names.
    #        The convention is as follows:
    #        filename = 'obspack_input.YYYYMMDD.nc'
    """
    # Create dataframe with Y M D H M S as columns with integer values, 
    # convert to numpy array. Contains entire campaign.
    df=pd.DataFrame() 
    df['Y'] = time.dt.strftime("%Y").astype(int)
    df['M'] =time.dt.strftime("%m").astype(int)
    df['D'] = time.dt.strftime("%d").astype(int)
    df['h'] =time.dt.strftime("%H").astype(int)
    df['m'] = time.dt.strftime("%M").astype(int)
    df['s'] = time.dt.strftime("%S").astype(int)
    
    # Time components in this format is the old format for time ObsPack
    # used as input in GEOS-Chem versions prior to 13.3.4 
    time_components = df.to_numpy() 
    
    # Get the Date only, no HMS so we can select points from a specifc day
    date_only=time.dt.strftime("%Y%m%d") # List of all days, for full time.
    indv_days= np.unique(date_only) # List of unique days 
    
    # Set the epoch that Obspack wants the difference in seconds from. 
    epoch=datetime.datetime(1970, 1,1,0,0,0)

    # Now get the total number of seconds between when you want samples and the epoch.     
    dts0=[ str((datetime.datetime(df.loc[i,'Y'],df.loc[i,'M'],df.loc[i,'D'],df.loc[i,'h'],df.loc[i,'m'],df.loc[i,'s']) - epoch).total_seconds()) for i in range(0,len(df)) ]

    # Make sure it stores times as integers and not datetime64[ns] objects .. 
    dts=np.array([np.int32(val.split('.0')[0]) for val in dts0],dtype=int)

    # Check to make sure the user passed values for lat/ long/ alt / time are valid. 
    if (any(abs(lat)> 90)) or (np.isfinite(lat).any()  ==False or (any(np.isnan(lat))==True)):
        sys.exit('Latitudes passed are either > 90 or contain non -Finite values')
    if (any(abs(lon)> 180)) or (any(np.isfinite(lon))  ==False) or (any(np.isnan(lon))==True):
        sys.exit('Longitudes passed are either > 180 or contain non -Finite values')
    if (any(alt < 0)) or (any(np.isfinite(alt))  ==False) or (any(np.isnan(alt))==True):
        sys.exit('Altitude passed are < 0 or contain non -Finite values')
    if (len(dts[np.isnan(dts)])> 0):
        print(dts[np.isnan(dts)])
        sys.exit('Times passed contain non -Finite values')
        
    # Create unique obspack array of ID strings: 
    n = np.arange(0, len(date_only)).astype(str)  # Create a uniq obs #, same len as time 
    ids= np.array([(sitename+ str(stri)) for stri in n]).astype('|S200') # Make IDs- same len as time
 
    str200= n[0].astype('|S200') # Make IDs- same len as time
    # ==================   Begin making Indv Files:  ==========================
    all_files = list()  # empty list to contain netcdf filenames generated.
    
    # Loop over each day you need an ObsPack File for. Grab time, lat lon, alt, ids, 
    # only on this date. Create the xarray and then save the .nc file. 
    for i in range(0, len(indv_days)): 
        # Get index of obs points on this day from larger arrays. 
        inds= np.where(date_only == indv_days[i])[0] 
        
        # Actual data you're putting inside the NetCDF file. 
        these_times= dts[inds].astype(int)
        these_ids=  ids[inds]

        # If all are same length then great, just use the index selected from time. 
        if all([len(v)==len(dts) for v in [lat, lon, alt]]):
            these_lats= lat[inds].astype(np.float32)
            these_lons= lon[inds].astype(np.float32)
            these_alts= alt[inds].astype(np.float32)
        else: 
            # If lat/long/alt not same len as time, then, assume you'll sample 
            # all these sites at each point in time. 
            these_lats=lat.astype(np.float32)
            these_lons=lon.astype(np.float32)
            these_alts=alt.astype(np.float32)
        
        # Round to OK # of digits! 
        these_lats=np.round_(these_lats,decimals=3) 
        these_lons=np.round_(these_lons,decimals=3) 
        these_alts=np.round_(these_alts,decimals=2) 
        
        # ====================================================================
        # =========    Make Big XArray with everything required.      ========
        # ====================================================================
        # Format of vars & attributes comes from GC example shown here:
        #   http://wiki.seas.harvard.edu/geos-chem/index.php/ObsPack_diagnostic
        
        ds = xr.Dataset({ # Need var "obs" which is just an indexing dimension
            'obs': xr.DataArray(
                data= np.arange(1, len(these_times)+1).astype(np.int64) ,
                dims=['obs'],
                attrs={
                    "long_name": "obs",
                    "_Storage":"chunked",
                    "_ChunkSizes": 1024,
                    "_Endianness":"little",
                }).astype(np.int64),
            # Need var "time" which gives the time you want to sample at.
            'time': xr.DataArray(
                data= these_times,
                dims=['obs'],
                attrs={
                    "units": "Seconds since 1970-01-01 00:00:00 UTC",
                    "_FillValue": -999999999,
                    "long_name": "Seconds since 1970-01-01 00:00:00 UTC",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(int),
            
            # Duplicate lat as many times as you collect obs points. Station
            # isnt' moving.
            'latitude': xr.DataArray(
                data=these_lats.astype(np.float32),
                dims=['obs'],
                attrs={
                    "units": "degrees_north",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "Sample latitude",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # Duplicate lon as many times as you collect obs points.
            'longitude': xr.DataArray(
                data=these_lons.astype(np.float32),
                dims=['obs'],
                attrs={
                    "units": "degrees_east",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "Sample longitude",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # Duplicate alt as many times as you collect obs points.
            'altitude': xr.DataArray(
                data=these_alts.astype(np.float32),
                dims=['obs'],
                attrs={
                    "units": "meters",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "sample altitude in meters above sea level",
                    "comment": "Altitude is surface elevation plus sample intake height in meters above sea level.",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # And now also pass our unique IDs (for each sample point).
            # Dimension is obs #, 200 len string.
            'obspack_id': xr.DataArray(
                data=these_ids,
                dims=['obs','string_of_200chars'],
                attrs={
                    "long_name": "Unique ObsPack observation id",
                    "comment": "Unique observation id string that includes obs_id, dataset_id and obspack_num.",
                    "_Storage":"chunked",
                    "_ChunkSizes": [1, 200],
                    "_DeflateLevel": 5,
                }),
            
            # Duplicate the int that indicates sampling stragety as many times
            # as you want to sample the model.
            'CT_sampling_strategy': xr.DataArray(
                data=np.full(len(these_times), sample_stragety, dtype=int),
                dims=['obs'],
                attrs={
                    "_FillValue": int(-9),
                    "long_name": "model sampling strategy",
                    "values": "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(int),
        })
        
        # ====================================================================
        # NOTE: FILES MUST BE NAMED LIKE THIS IN ORDER TO BE READ IN GEOS CHEM
        #       OBSPACK_MOD.F90
        filename = 'obspack_input.'+ str(indv_days[i])+ '.nc'
        # ====================================================================
    
        # Convert our Xarray Data set to a netcdf file w/ this name at outpath
        ds.to_netcdf(outpath + filename)
        
        # Append name of this file to our list of file names that we return
        # to the user.
        all_files.append(outpath + filename)
        
        blank='========================================================================\n'
        if quiet is False: 
            if i ==0: print(blank, 'Obspack files saved at : ', outpath, blank)
            print(filename)
            if i ==len(indv_days): print(blank, 'Done!', blank)
        
    return all_files
    

def write_obspack_inputs_ground(sitename: str, lat: int, lon: int, alt: int,
                               datestart: str, dateend: str, samplefreq: int,
                               sample_stragety: int, outpath: str, quiet: bool =False):
    """Create ObsPack input files for GEOS Chem v13.0.0. for a single gound site.
    
    # =========================================================================
    #                             INPUTS
    # =========================================================================
    sitename = string with sitename used in obspack ID, (e.g.'SOAS-Tower')
    
    lat = integer with latitude of obs site in degrees north (e.g. 35.45)
   
    lon = integer with longitude of obs site in degrees east
    
    alt = integer with elevation plus sample intake height in meters above sea level
    
    datestart = string format YYYYMMDD HH:MM:SS' indicating when you need
                  obsPack netcdf  file to begin
                  
    dateend = string format YYYYMMDD HH:MM:SS' indicating when you need
                  obsPack netcdf  file to stop
                  
    freq = integer for # of seconds you'd like to step betweeen datestart,
            dateend (e.g. 3600 for hourly steps)
            
    Sample Stragety = Integer specifiing how the model will average to
                      the timebase given (calc'd by  datestart,dateend,
                      and freq). Only Valid option are:
                                1  =   4-hour avg;
                                2  =   1-hour avg;
                                3  =   90-min avg;
                                4  =   instantaneous
    
    outfile_path  = string containing absolute path where the netcdf files
                    will be written.
    
    quiet = boolean of whether you'd like to print progress to screen.
    
    # =========================================================================
    #                             OUTPUTS
    # =========================================================================
    
      out= Returns list of netcdf files written to outpath for each day in
          the time period that  can be used as GEOS Chem inputs.
    
    #  NOTE: The output files MUST be named with YYMMDD in order to be
    #        read in by GEOS-Chem properly. So do NOT change the names.
    #        The convention is as follows:
    #        filename = 'obspack_input.YYYYMMDD.nc'
    
    # =========================================================================
    #                              Example
    # =========================================================================
    
    SOAS Centerville Ground collection site was at lat=32.903281,
    lon=-87.249942, alt=125. Let's sample the model there ever hour
    (outputting the hourly average).This code snippiet will make ALL the
    .nc files we need for ObsPack input for the length of our run (6/1/2013-
    7/15/2013) and save them to my desktop.
    
    import obspack_io as obs 
    import xarray
    
    # SOAS ground site at  lat=32.903281, lon= -87.249942, alt=123m elevation + 2m sample height 
    filenames = obs.write_obspack_inputs_ground('SOAS', lat=32.903281, 
                                          lon=-87.249942, alt=125,
                                          datestart='20130601 00:00:00',
                                          dateend='20130715 00:00:00', 
                                          samplefreq=3600,
                                          sample_stragety=2, 
                                          outpath='~\\SOAS\\ObsPack_Inputs\\')
    
     # Open the first file with xarray, not allowing xarray to changing var types
     # as you read in the files, and print info about what's in it! 
     dat = xr.open_dataset(filename[0], mask_and_scale=False, decode_times=False)
     print(dat)
    
    """
    # Check to make sure passed values are valid. 
    if (abs(lat)> 90) or  (np.isfinite(lat) is False):
        sys.exit('Latitudes for '+sitename+'are either > 90 or contain non -Finite values')
        
    if (abs(lon)> 180) or  (np.isfinite(lon)  is False):
        sys.exit('Longitudes for '+sitename+' are either > 180 or contain non -Finite values')
        
    if (alt < 0) or  (np.isfinite(alt)  is False):
        sys.exit('Altitudes for '+sitename+' are < 0 or contain non -Finite values')
    
    # Round to OK # of digits! 
    lat=np.round_(lat,decimals=3) 
    lon=np.round_(lon,decimals=3) 
    alt=np.round_(alt,decimals=2) 
        
    # Figure out # of days between start & Enddate (e.g. # of NC files we need)
    start = datetime.datetime.strptime(datestart, '%Y%m%d %H:%M:%S')
    endd = datetime.datetime.strptime(dateend, '%Y%m%d %H:%M:%S')
    delta = endd - start
    
    if (delta.seconds != 0):
        days = delta.days + 1  # add a day if there's some leftover time
    else:
        days = delta.days  # otherwise stick with the round number.
    
    starts = pd.date_range(datestart, dateend, freq='1D').to_series()
    starts2Use = starts[:days]  # only use the # of days we decided on.
    ends2Use = starts2Use + datetime.timedelta(days=1) - datetime.timedelta(seconds=1)
    
    # Don't do the whole day if the inputted stop time is shorter.
    if ends2Use.iloc[-1] > endd: ends2Use.iloc[-1] = endd
    
    # ===================================================================
    all_files = list()  # empty list to contain netcdf filenames generated.
    
    # Loop over # of days we need a file for and make a netcdf!
    for t in range(0, len(starts2Use)):
        
        # Get datetimes in our day at sample freq
        sdates = pd.date_range(starts2Use.iloc[t], ends2Use.iloc[t], freq=str(samplefreq) + 's')
        
        # Set the epoch that Obspack wants the difference from... 
        epoch=datetime.datetime(1970, 1,1,0,0,0)
        
        # Now get the total number of seconds between when you want samples and the epoch. 
        dts0= (sdates-epoch).total_seconds().map(str)
        
        # Make sure it stores times as ints and not datetime64[ns] objects .. 
        dts=np.array([np.int32(val.split('.0')[0]) for val in dts0],dtype=int)
        
        df=pd.DataFrame() 
        df['Y'] = sdates.strftime("%Y").astype(np.int64)
        df['M'] =sdates.strftime("%m").astype(np.int64)
        df['D'] = sdates.strftime("%d").astype(np.int64)
        df['h'] =sdates.strftime("%H").astype(np.int64)
        df['m'] = sdates.strftime("%M").astype(np.int64)
        df['s'] = sdates.strftime("%S").astype(np.int64) 
        
        # Time components in this format is the old format for time ObsPack
        # used as input in GEOS-Chem versions prior to 13.3.4 
        time_components = np.round(df.to_numpy()).astype(np.float64)
        
        # =====================================================================
        # Create unique obspack array of ID strings:
        n = np.arange(0, len(dts)).astype(str)  # Create a uniq obs #
        
        # Create a prefix containing the sitename, start/stop dates
        prefix = sitename + '_from_' + \
            datestart.split(' ')[0] + '_to_' + dateend.split(' ')[0] + '_n'
            
        # Make obspack IDs. ... old method before needed new dim for string200! 
        ids = np.array([(prefix + stri) for stri in n]).astype('|S200')
        
        # ====================================================================
        # =========    Make Big XArray with everything required.      ========
        # ====================================================================
        # Format of vars & attributes comes from GC example shown here:
        #   http://wiki.seas.harvard.edu/geos-chem/index.php/ObsPack_diagnostic
        
        ds = xr.Dataset({
             # Need var "obs" which is just an indexing dimension
            'obs': xr.DataArray(
                data= np.arange(1, len(dts)+1).astype(int) ,
                dims=['obs'],
                attrs={
                    "long_name": "Sample latitude",
                    "_Storage":"chunked",
                    "_ChunkSizes": 1024,
                    "_Endianness":"little",
                }).astype(int),
            # 'calendar_components': xr.DataArray(
            #     data= np.arange(0,6).astype(np.int64) ,
            #     dims=['calendar_components'],
            #     ).astype(np.int64),
            # Need var "time" which gives the time you want to sample at.
            'time': xr.DataArray(
                data= dts.astype(int),
                dims=['obs'],
                attrs={
                    "units": "Seconds since 1970-01-01 00:00:00 UTC",
                    "_FillValue": -999999999,
                    "long_name": "Seconds since 1970-01-01 00:00:00 UTC",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(int),
            
            # Duplicate lat as many times as you collect obs points. Station
            # isnt' moving.
            'latitude': xr.DataArray(
                data=np.full(len(dts), lat, dtype=np.float32),
                dims=['obs'],
                attrs={
                    "units": "degrees_north",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "Sample latitude",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # Duplicate lon as many times as you collect obs points.
            'longitude': xr.DataArray(
                data=np.full(len(dts), lon, dtype=np.float32),
                dims=['obs'],
                attrs={
                    "units": "degrees_east",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "Sample longitude",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # Duplicate alt as many times as you collect obs points.
            'altitude': xr.DataArray(
                data=np.full(len(dts), alt, dtype=np.float32),
                dims=['obs'],
                attrs={
                    "units": "meters",
                    "_FillValue": np.float32(-1.e+34),
                    "long_name": "sample altitude in meters above sea level",
                    "comment": "Altitude is surface elevation plus sample intake height in meters above sea level.",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(np.float32),
            
            # And now also pass our unique IDs (for each sample point).
            # Dimension is obs #, 200 len string.
            'obspack_id': xr.DataArray(
                data=ids,
                dims=['obs'],
                attrs={
                    "long_name": "Unique ObsPack observation id",
                    "comment": "Unique observation id string that includes obs_id, dataset_id and obspack_num.",
                    "_Storage":"chunked",
                    "_ChunkSizes": 1,
                    "_DeflateLevel": 5,
                }),
            
            # Duplicate the int that indicates sampling stragety as many times
            # as you want to sample the model.
            'CT_sampling_strategy': xr.DataArray(
                data=np.full(len(dts), sample_stragety, dtype=int),
                dims=['obs'],
                attrs={
                    "_FillValue": int(-9),
                    "long_name": "model sampling strategy",
                    "values": "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous",
                    "_Storage":"chunked",
                    "_ChunkSizes": 778,
                    "_DeflateLevel": 5,
                    "_Endianness":"little",
                }).astype(int),
        })
        # ====================================================================
        # NOTE: FILES MUST BE NAMED LIKE THIS IN ORDER TO BE READ IN GEOS CHEM!
        filename = 'obspack_input.'+ str(starts2Use.iloc[t]).split(' ')[0].replace('-', '') + '.nc'
        # ====================================================================
    
        # Convert our Xarray Data set to a netcdf file w/ this name at outpath
        ds.to_netcdf(outpath + filename)
        ds.close()
        
        # Append name of this file to our list of file names that we return
        # to the user.
        all_files.append(outpath + filename)
        
        blank='========================================================================\n'
        if quiet is False: 
            if t ==0: print(blank, 'Obspack files saved at : ', outpath, blank)
            print(filename)
            if t ==len(starts2Use): print(blank, 'Done!', blank)
        
    return all_files


def combine_common_obspack_inputs(folder_1: str , folder_2:str , outpath: str, 
                            copy_not_common : bool = False):
    """Look in two folders for ObsPack Files and combine common date files.
    
    This function is useful if you want to sample the model on the same date 
    at multiple ground sites. For example, maybe  I want to 
    sample GEOS-Chem at the 2011 BEACHON-RomBAS site and at the ALC 2011 site 
    which took place in different grid boxes on the same days. 
    I used the write_obsPack_inputs_ground() function to make ObsPack files for each campaign
    individually, and then this function to combine the files for the same dates (in time order). 
    
    You could seperate the output by either using the lat/lon of the ground sites 
    or probably, easier, by filtering by the obspack_ID string (which is unique
    for those from BEACHON vs. ALC)
    
    Args: 
    ----
        folder_1: string path to a folder containing ObsPack files. 
                    (outpath for make Obspack Inputs)
        folder_2: string path to a folder containing ObsPack files. 
                    (outpath for make Obspack Inputs)
        outpath: string to path you want to store the combined ObsPack files in.
        
        copy_not_common: Boolean if you want to copy files that the two folders
            don't have in common to this outpath. So have a folder of all 
            ObsPack files from Folder 1, folder 2, with the commons combined 
            and the not commons just copied over from where they came from. Default is False. 
            
    Output: 
    ------
        files_combo= list of files that were combined and saved at outpath.
    """
    # Get all file names in these source folders: 
    filenames_1=glob.glob(folder_1+"*.nc")
    filenames_2=glob.glob(folder_2+"*.nc")
    
    # Split the path of the filename, and compare just the dates of the files: 
    dates1= [filenames_1[i].split('obspack_input.')[1] for i in range(0,len(filenames_1))]
    dates2= [filenames_2[i].split('obspack_input.')[1] for i in range(0,len(filenames_2))]
    
    # Get list of Common dates we have 2 ObsPack files for: 
    common = list(set(dates1) & set(dates2))
    common = ['obspack_input.' + sub for sub in common] 
    
    # Loop over each file we have 2 ObsPack files for: 
    for i in range(0, len(common)):  
        # Full path of each file for the same date: 
        indv_filename1= folder_1+ common[i]
        indv_filename2= folder_2+ common[i]
        
        # Open each data set: 
        one = xr.open_dataset(indv_filename1, mask_and_scale=False, decode_times=False)
        two = xr.open_dataset(indv_filename2, mask_and_scale=False, decode_times=False)
        
        # Combine the two data sets: 
        both= xr.concat([one, two], dim='obs')
        
        # Get the indicies that would sort them in time: 
        ind= np.argsort(both.time.values.astype(int)) 
        
        # Re-Index "Both" so that everything is in time order.
        both.time.values = both.time.values[ind]
        both.latitude.values = both.latitude.values[ind]
        both.longitude.values = both.longitude.values[ind]
        both.altitude.values = both.altitude.values[ind]
        both.obspack_id.values = both.obspack_id.values[ind]
        both.CT_sampling_strategy.values = both.CT_sampling_strategy.values[ind]
        
        # Convert combined DataSet to a netcdf file w/ this name at outpath
        both.to_netcdf(outpath+ common[i])
        
    print('Common files saved at: ', outpath) # Let 'em know.
    
    if copy_not_common is True: 
        # Get list of dates we don't have 2 ObsPack files for: 
        not_common = list(set(dates1) ^ set(dates2))
        not_common = ['obspack_input.' + sub for sub in not_common] 
        for i in range(0, len(not_common)):  
            # Two possible paths to this file in each folder: 
            indv_filename1= folder_1+ not_common[i]
            indv_filename2= folder_2+ not_common[i]
            
            # Copy it over to the new outpath from its source when you find it: 
            if indv_filename1 in filenames_1: 
                copyfile(indv_filename1, outpath+not_common[i])
                
            if indv_filename2 in filenames_2: 
                copyfile(indv_filename2, outpath+not_common[i])
                
        print('Uncommon files copied over to: ', outpath)
                
    return common

#############################################################################
########################   Output Handling    ###############################
#############################################################################


def sep_obspack_output(path_to_dir:str , split_id_at: str = '_from_', 
                            file_patterns: str = ['ObsPack', '.nc4'], 
                            quiet:bool =False, debug: bool = False):
    
    """Parse all obspack files in a directory and seperate their contents based on 
    their obspack_ids. Useful if you sampled the model at multiple different 
    places in one Obspack Input file but want the output seperated. 
    
    Args:
    -----
        path_to_dir = String with filepath to folder containing GEOS Chem output 
                      Obspack netcdf files
                      
        split_id_at (optional) = String you'd like to split the ObsPack ID strings
                on where whatever preceeds this string is the "campaign" name prefix... 
                Default is to split ID's on'_from_' which is the naming convention used in 
                write_obspack_inputs() to careate obspack_id strings. 
            
        file_patterns (optional) = List of things that much be in a filename in order 
                to actually concatenate them. (Useful for filtering by date,
                file extension, etc).Default is to look for .nc4 files that contain
                "ObsPack" as they would as GEOS-Chem output. 
                
        debug (optional) = Boolean of whether to print metrics that can help 
                determine whether this seperation is being done correctly or not. 
                Default is set to False. 
                
        
    Output: 
    -------
    
        New ObsPack ouput files contained in folders named 
        by the Obspack-ID they came from at path_to_dir
        
    """
    # Look for all the ObsPack netCDF files in the directory given 
    file_list = ut._find_files_in_dir(path_to_dir, file_patterns)
    
    # Get all inv names of the files in the directory.  
    names=[file.split(path_to_dir)[1] for file in file_list]
    
    save_dirs=list([])# To hold output folders that you make.. 
    for f,file in enumerate(file_list): # Loop over all files in the directory. 
        ds=xr.open_dataset(file, mask_and_scale=False, decode_times=False)
        
        # Split the obspack IDs to get the sitename / prefixes only. 
        prefs=[idn.decode().split(split_id_at)[0] for idn in ds.obspack_id.values]
        
        camps= np.unique(prefs) # Get list of unique campaign prefixes. 
        
        # Notify User that only one campaign was found in this particular file... 
        if len(camps) <=1: 
            if quiet is False: print("sep_obspack_output.py ... Found only: ", camps, 'in: ', names[f])
            if debug is True: 
                print("All ObsPack IDs contained in this file:")
                [print(str(ids)) for ids in ds.obspack_id.values]
        
        for n in range(0,len(camps)): 
            ds_n=ds.copy()  # Make a new copy of the original ds for each new campaign. 
            
            # Get the IDs of all the Obspack vals that don't have this campaign in their name
            n_ids=[ij for ij, ids in enumerate(ds.obspack_id.values) if camps[n] not in ids.decode()]
           
            ds_n['indx'] = ("obs", np.arange(0,len(ds.obspack_id.values)))
    
            # Tell xarray that we want "indx" to be a coordinate.
            ds_n = ds_n.set_coords('indx') # because obs contains repeated values if both sampled at same time. 
            
            # And tell it to replace obs # with time as the preferred dimension.
            ds_n = ds_n.swap_dims({"obs": "indx"})
    
            # drop all the rows from the data set if they're from a different campaign. 
            ds_n=ds_n.drop_sel(indx=n_ids)
            
            #Switch the index back to obs and drop the "index". 
            ds_n = ds_n.swap_dims({"indx": "obs"}).drop("indx")
            
            # Specify the path where this campaign's data will be saved. 
            fld= os.path.join(os.path.dirname(path_to_dir), camps[n]+'//')
            tf= os.path.isdir(fld) # Check if folder exists... 
            if not tf: os.mkdir(fld) # Make directory to hold output from this campaign only
            save_dirs.append(camps[n]+'//') 
            out= os.path.join(fld+names[f])#  Join the output folder and filename together. 
            ds_n.to_netcdf(out) # And save the data from this campaign only. 
            ds_n.close()
            
        ds.close()
            
    return save_dirs


def _switch_obs_2_time_dim(ds, seperate: bool = False):
    """Function to create a single time variable that is the midpoint of the 
    ObsPack averaging interval, and make it the xarray coordinate. """ 
    #print(ds.data_vars)
            
    # Get the midpoint of the average pulled from the model:
    midpoint = pd.to_datetime((ds.averaging_interval_start.values +ds.averaging_interval.values/2))
    
    # Make the time midpoint a new variable in the dataset.
    t = midpoint.to_series().reset_index(drop=True)
    ds['time'] = ("obs", t)
    
    # Tell xarray that we want time to be a coordinate.
    ds = ds.set_coords('time')
    
    # And tell it to replace Obs # with time as the preferred dimension.
    ds = ds.swap_dims({"obs": "time"})
    
    return ds

    
def read_and_concat_output_files(path_to_dir: str, outdir: str = None,
                   outfile: str = None, match_pattern= ['ObsPack', '.nc4']):
    """
    Concatonate output ObsPack files into a single file from a directory.
    # If you only want to open a single file, try xarray.open_dataset().
    
    Args
    ----
        path_to_dir = String with filepath to folder containing GEOS Chem output 
                      Obspack netcdf files
    
        outdir(optional) = # String path to where concatonated file will be saved 
    
        outfilename (optional) = string name of output file, no extension needed.
        
        match_pattern (optional) = List of things that much be in a filename in order 
                to actually concatenate them. (Useful for filtering by date,
                file extension, etc).Default is to look for .nc4 files that contain
                "ObsPack" as they would as GEOS-Chem output. 
    """
    # If outdir not set, then set it to the same as the input file path.
    outdir = path_to_dir if outdir is None else outdir
    
    # If outfilename not set, then set it to be concat_ObsPack.nc
    outfile = 'concat_ObsPack' if outfile is None else outfile
    
    # Get a list of variables that GCPy should not read.
    # These are mostly variables introduced into GCHP with the MAPL v1.0.0
    # update.  These variables contain either repeated or non-standard
    # dimensions that can cause problems in xarray when combining datasets.
    skip_vars = ["anchor","ncontact","orientation","contacts","cubed_sphere"]
    
    # Look for all the ObsPack netCDF files in the path
    file_list = ut._find_files_in_dir(path_to_dir, match_pattern)
    
    # Return a single xarray Dataset containing data from all files
    # Specify a pre-processing function that makes time midpoint the coordinate.
    ds = xr.open_mfdataset(file_list, drop_variables=skip_vars,
                                 preprocess=_switch_obs_2_time_dim)
    ds.close()
    # Specify the path and filename for the concatenated data
    out = str(os.path.join(outdir, outfile+'.nc'))
    
    # Write concatenated data to a netCDF file
    ds.to_netcdf(out)
    
    print('File saved to:' + out) 
    
    return out 
