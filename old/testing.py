#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 11:59:22 2024

@author: u6044586
"""
import pandas as pd 
import datetime 
import numpy as np 
import xarray as xr

datestart='20130601 00:00:00' # When I want to start sampling 
dateend='20130602 00:00:00'  # When I want to end sampling 
samplefreq= 1800 

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

# Get datetimes in our day at sample freq
sdates = pd.date_range(starts2Use.iloc[0], ends2Use.iloc[0], freq=str(samplefreq) + 's')

# Set the epoch that Obspack wants the difference from... 
epoch=datetime.datetime(1970, 1,1,0,0,0)

# Now get the total number of seconds between when you want samples and the epoch. 
dts0= (sdates-epoch).total_seconds().map(str)

# Make sure it stores times as ints and not datetime64[ns] objects .. 
dts=np.array([np.int32(val.split('.0')[0]) for val in dts0],dtype=int)

# =====================================================================
# Create unique obspack array of ID strings:
n = np.arange(0, len(dts)).astype(str)  # Create a uniq obs #



for stri_i, stri in enumerate(n[0:6]):
    st='SOAS'+stri
    id_i=np.array([char for char in st]).astype('|S200') 
    
    if stri_i==0: 
        ids= np.array(id_i,dtype='|S200')
    else: 
        ids=np.vstack((ids,id_i),dtype='|S200')

ids=ids.astype('|S200')
#str200= n[0].astype('|S200') 




ds = xr.Dataset({ # Need var "obs" which is just an indexing dimension
    'obs': xr.DataArray(
        data= np.arange(1, 7).astype(int) ,
        dims=['obs'],
        attrs={
            "long_name": "obs",
            "_Storage":"chunked",
            "_ChunkSizes": 1024,
            "_Endianness":"little",
        }).astype(int),
    'string_of_200chars': xr.DataArray(
        data= ids[0],
        dims=['string_of_200chars']
        ),
    'obspack_id': xr.DataArray(
        data=ids,
        dims=['obs','string_of_200chars'],
        attrs={
            "long_name": "Unique ObsPack observation id",
            "comment": "Unique observation id string that includes obs_id, dataset_id and obspack_num.",
            "_Storage":"chunked",
            "_ChunkSizes": [1, 200],
            "_DeflateLevel": 5,
        }),
    })


    
    
    