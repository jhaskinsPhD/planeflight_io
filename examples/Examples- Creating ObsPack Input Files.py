# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 10:45:25 2022

@author: jhask
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 10:40:31 2021

@author: jhask
"""
cccccccccccccccc

# ======   Open the SENEX data, and make sure its formatted OK!=============

# Open SENEX merged data. 
senex=path_to_examples+'\\datafiles_for_examples\\SENEX.nc'
sn=xr.open_dataset(senex) 
 
# Create a smaller dataset from the larger one, with ONLY the stuff that Obspack needs,
# And drop any rows that have a NaN for any of these vars in them, from the dataset. 
sen =  xr.Dataset({'alt':xr.DataArray( data= sn.GpsAlt.values , dims=['obs']),  
                   'lat':xr.DataArray( data= sn.GpsLat.values , dims=['obs']), 
                   'lon':xr.DataArray( data= sn.GpsLon.values , dims=['obs']), 
                   'time':xr.DataArray(data= sn.time.values ,  dims=['obs'])}).dropna(dim='obs') 

# Convert times to a ** pandas datetime array ** 
time= pd.to_datetime(sen.time.values).to_series().reset_index(drop=True) 

# =============================================================================
# ########### EXAMPLE #4: Sample SENEX Flights Using ObsPack   ################
# =============================================================================
#  Make Obspack Files to sample the model along all of SENEX's flight tracks. 

tf= os.path.isdir(path_to_examples+'\\example4\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example4\\') # Make directory to hold output.
    
# Make ObsPack files for all flights during SENEX 
# Pass lat, long, alt and time as arrays, set sample stragety=4 for instantantoues sampleling of the model. 
filenames = obs.write_obspack_inputs_flights('SENEX',lat=sen.lat.values,  
                                                     lon=sen.lon.values, 
                                                     alt=sen.alt.values, 
                                                     time=time,
                                                     sample_stragety=4, 
                                                     outpath=path_to_examples+'\\example_4\\') 
 
# ============================================================================= 
# ######## EXAMPLE #5: Sample SOAS Centreville Ground Site UsSing ObsPack ######
# ============================================================================= 
# Create ObsPack files to sample the model at the the Centreville, Alabama ground 
# site during SOAS sampling the model there every hour of the campaign between 
# 6/1/2013 and 7/15/2013. 

tf= os.path.isdir(path_to_examples+'\\example5\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example5\\') # Make directory to hold output.
    
filenames = obs.write_obspack_inputs_ground('SOAS-Ground', 
                                      lat=32.903281,  # lat of the sample site 
                                      lon=-87.249942,  # lat of the sample site 
                                      alt=2, # height above ground level to sample at 
                                      datestart='20130601 00:00:00', # When I want to start sampling 
                                      dateend='20130715 00:00:00',   # When I want to end sampling 
                                      samplefreq=3600, # How frequently the model should sample  in s
                                      sample_stragety=2,  # Set it to give me an avg hourly output when it samples.
                                      outpath=path_to_examples+'\\example_5\\') 
 
# ============================================================================= 
#  ####### EXAMPLE #6: Sample SENEX and SOAS simultaneously using ObsPack  ####
# ============================================================================= 
# Suppose you want to sample GEOS-Chem at the SOAS Centerville site and  
# along the SENEX flight path during June- July of 2013. Some of the SENEX flights 
# took place on the same days as the SOAS continous ground monitoring, so 
# in example #4 and #5 above, we have a few Obspack files that are for the same date. 
# GEOS-Chem will not allow you to have multiple Obspacks files for the same date, 
# as input, so you need to combine the sampling data in those together to sample the 
# model at both. We achieve this by, FIRST, sampling at each individually, as we did 
# in the prior two examples, and then looping through out files, finding when they 
# are on the same dates, and combining that input. We can seperate them after 
# GEOS-Chem outputs stuff using the Obspack ID (they're different for SOAS/SENEX). 
 
tf= os.path.isdir(path_to_examples+'\\example6\\') # Check if folder exists... 
if not tf: os.mkdir(path_to_examples+'\\example6\\') # Make directory to hold output.
     
# Folders to each of their indvidual ObsPack Files  
senex_folder = path_to_examples+'\\example_4\\' 
soas_folder = path_to_examples+'\\example_5\\' 
combo_folder= path_to_examples+'\\example_6\\'
 
# Tell it to copy over the not-common ObspackFiles (and combo the common ones!) 
commons= obs.combine_common_obspack_inputs(soas_folder, senex_folder, 
                                     outpath=combo_folder,  
                                     copy_not_common= True) 

# The files you'd pass to GEOS-Chem as input to sample SENEX AND SOAS, are now in combo_folder.
# Some files tell it to sample just SOAS, others tell it to sample both. We have a 
# function to seperate the output when it's sampled like this too! 

# Now, open a combined" file we generated  and take a peek for sanity!  
t=xr.open_dataset(combo_folder+commons[0]) 
 
lat= t.latitude.values 
lon= t.longitude.values 
alt= t.altitude.values 
tm= t.time.values 
obspack_id= str(t.obspack_id.values) 
sample= t.CT_sampling_strategy.values

