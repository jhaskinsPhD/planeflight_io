#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Internal functions not intended for direct user usage called within main functions 
in planeflight_io.py 

Created on Thu Feb 27 21:13:59 2025

@author: Dr. Jessica D. Haskins
GitHub: @jhaskinsPhD
Email: jessica.haskins@utah.edu 
"""

import re
import os
import sys
import warnings
import yaml
import numpy as np
import pandas as pd 
import xarray as xr 

###############################################################################
# INTERNAL HELPER FUNCTIONS- FOR PLANEFLIGHT_IO 
###############################################################################
def _check_simtype(simtype:str, blank_allowed:bool=False): 
    """Internal function to check that argument passed for "simtype" is allowed/ 
    recognized (case insensitive).  
    
    Inputs: 
    ------ 
         (1) simtype - STRING containing user's simulation type.(e.g. any option 
                       assigned for 'sim_name' in createRunDir.sh
         
         (2) blank_allowed - BOOL indicating if blank input for simtype should 
                             result in ValueError or not. If TRUE, will bypass 
                             error. If FALSE (Default) will throw error if blank. 
                            
    """
    # Check that 'simtype' is one of the allowed values (e.g. any option assigned for 
    # 'sim_name' in createRunDir.sh or blank (returns all not broken!))
    allowed=['fullchem','aerosol','carbon','Hg','POPs','tagO3','TransportTracers',
             'metals','CH4','CO2','tagCO']
    
    # Add "blank" simtype to allowed list if blank_allowed is set to True. 
    if blank_allowed==True:  allowed.append('') 
        
    # Case insensitive check in syptpe is in allowed. 
    if simtype.lower() not in [name.lower() for name in allowed]: 
        raise ValueError("Invalid input for 'simtype'='"+simtype+"'"+ 
                          'Valid options include (case insensitive): \n\t'+ 
                          ','.join(allowed)).with_traceback(sys.exc_info()[2])
        
    return 


def _parse_gc_config(config_yaml: str):
    """Internal function to parse the geoschem_config.yml files and return
    the simulation type and a list of all advected species. Primarily called from 
    within make_planeflight_inputs() when a geoschem_config/yml file is passed to 
    create list of tracers that planeflight will output and also to determine
    what other planeflights diagnostics are compatible with user's simulation type.
    
    INPUTS:
    -------
        (1) config_yaml -  A STRING containing path to geoschem_config.yaml file. 
        
    OUTPUTS: 
    ------- 
        (1) adv_species - A LIST of all transported species 
        
        (2) simtype - A STRING containing the simulation type (e.g. 'fullchem') 
    """
    
    # Check that the geoschem_config.yaml files exists, otherwise throw an error: 
    if not os.path.isfile(config_yaml):
        raise FileNotFoundError('The following file could not be found:\n\t'+
                                config_yaml).with_traceback(sys.exc_info()[2])
    # Open the geoschem_config.yml file and load contents as dict: 
    with open(config_yaml, 'r') as f:
        gc_config=yaml.load(f, Loader=yaml.FullLoader)
    
    # Extract a list of all transported species from the geoschem_config.yaml dict:  
    adv_species=gc_config['operations']['transport']['transported_species']
    
    # Extract the simulation type from the geoschem_config.yaml dict: 
    simtype=gc_config['simulation']['name']
    
    return adv_species, simtype


def _read_spdb_yml(spdb_yaml):
    """Function to read in the geoschem species database yaml file"""
    # Check that the species database file exists, otherwise throw an error: 
    if not os.path.isfile(spdb_yaml):
        raise FileNotFoundError('The following file could not be found:\n\t'+
                                spdb_yaml).with_traceback(sys.exc_info()[2])
    
    # Load the species database as a dict from the species_database.yml file: 
    with open(spdb_yaml, 'r') as f:
        spdb=yaml.load(f, Loader=yaml.FullLoader)
        
    return spdb


def _read_planeflight_diags_yml(inputs_only:bool=False, no_collections:bool=False): 
    """Internal function to read in the planeflight_diagnostics.yml file 
    as a dictionary which contains "collections" of all vars outside of 
    advected species that you can "request" to output from planeflight.
    
    Inputs: 
    -------
    
        (1) inputs_only - BOOL indicating if you ONLY want to return diagnostics 
                          valid for creating planeflight files. Default is set to FALSE (e.g. 
                          returns valid input diags to create files & the standard 
                          output collection (which planeflight adds to all output files) 
                          
        (2) no_collections  - BOOL indicating if you want to return diags organized by
                        "collection" or not. Default is FALSE (to include them
                        arranged by collection) in a nested dict format where first
                        key is collection name. Setting to TRUE removes the nested
                        dict organizing them by collection so diagnostic name is
                        first key.
    
    Outputs: 
    ------- 
       (2) diags - DICTIONARY with contents of planeflight_diagnostics.yml
       
    NOTE: All info in planeflight_diagnostics.yml read in here is compiled from 
          within planeflight_mod.F90 and from the Planefligth Diagnostic ReadTheDocs page: 
          https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/planeflight.html
    
    """
    # Get path to planeflight_diagnostics.yml (stored in same dir as this file) 
    diag_yaml = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'planeflight_diagnostics.yml')
    
    # Check that the planeflight_diagnostics.yml files exists at this path:  
    if not os.path.isfile(diag_yaml):
        raise FileNotFoundError('Could not find required file, "planeflight_diagnostics.yml".' 
                                'This function expected to find it at: '+
                                diag_yaml).with_traceback(sys.exc_info()[2])
        
    # Open the planeflight_diagnostics.yml containing all optional planeflight 
    # diagnostics including names/dictionaries. 
    with open(diag_yaml, 'r') as f:
        diags=yaml.load(f, Loader=yaml.FullLoader)
        
    # Remove "std_out" collection if inputs_only set to TRUE. 
    if inputs_only is True: 
        #This collection is  included in planeflight.yml so we can assign attrs 
        # for those vars if _convert_to_xarr() is called. But, if asked to retrieve 
        # only valid INPUT diagnotisic options, we will remove it 
        # (since you can't request these vars in planeflight input files). 
        diags.pop('std_output') 
        
    if no_collections is False:
        # Standard/default output has all diags arranged by collections. return that.
        return diags
    else: 
        # Remove collections nested dict outputs and just return the diags themselves. 
        diag_dict=dict({})
        for collection in list(diags.keys()): 
            for key in diags[collection]['Diagnostics']: 
                diag_dict[key]=diags[collection]['Diagnostics'][key]
                
        return diag_dict 


def _get_tracer_name_num_mapping(config_yaml, keys_are_nums:bool=True):
    """Internal helper function to create dictionary associating tracer name with 
    tracer index numbers. Created by assuming order of species listed in the 
    geoschem_config.yml file are the order of the tracer numbers (which appears to 
    line up with doing this isntead by parsing the geoschem output log file...)
    
    (1) config_yaml -  A STRING containing path to geoschem_config.yml file. 
    
    (2) keys_are_nums - BOOLEAN indicating if you want the output dictionary to 
                           have tracer numbers as keys (TRUE) or as values (FALSE).
                           If set to TRUE, output dict, name_map, will be of the form: 
                              name_map['TRA_052']='ClNO2', 
                           If set to FALSE output dict, name_map, will be of the form:  
                              name_map['ClNO2']='TRA_052'. 
                           Default is set to TRUE (output dict keys are tracer #s). 
    OUTPUT: 
    -------
        
        (1) name_map     - DICTIONARY with key/value pairs of tracer names (ClNO2) & 
                           tracer numbers (e.g. TRA_052). 
    """
    # Get list of all advected species 
    adv_species, _ = _parse_gc_config(config_yaml)
    
    # Create dict to map advected species names to tracer numbers... 
    name_map=dict({});
    for i,tra_name in enumerate(adv_species): 
        tra_num= f'TRA_{i+1:03}'
        if keys_are_nums is True: 
            name_map[tra_num]=tra_name
        else: 
            name_map[tra_name]=tra_num
        
    return name_map


def _display_diags(diags:dict, collections:list): 
    """Internal function to parse "collections" within dictionary "diags" 
    and create a single string that can be printed to screen displaying info 
    about all the diagnostic collections within input dict, "diags". 
    """
    
    lines=[] # Empty list to hold lines. 
    
    # Loop over each collection & pull out info/format, append to lines: 
    for coll in collections: 
        lines.append("-"*50)
        lines.append(f"{diags[coll]['Collection']} Collection")
        lines.append("-"*50)

        lines.append('Diagnostics:')
        for key in list(diags[coll]['Diagnostics'].keys()):
            lines.append(f"\t{key:<10} = {diags[coll]['Diagnostics'][key]['FULLNAME']}")
        if 'Notes' in list(diags[coll].keys()): 
            lines.append(f"Notes:\n\t{diags[coll]['Notes']}")
    
    # Join all lines together with line breaks & output single string:     
    info= '\n'.join(lines)
    
    return info  


def _check_str_arr_inputs(item, item_nm:str, allowed_vals, allowed_nm:str,
                          wildcard_str:bool=False, debug:bool=False ): 
    """ Internal helper function to check that input for "item" is EITHER a 
    string wildcard  matching '?ALL?' (if wildcard_str) is set to TRUE, OR is
    array like, with all array elements as type(string), and with all elements
    are within "allowed_vars". Returns if an error is found & why the error happened,
    Function only called within make_planeflight_inputs() to help with user 
    input validation & has not other valid usage. 
    
    """ 
    # Initialize stuff. 
    is_error=True;err_reason='unknown'; give=allowed_vals; 
    
    if type(item)==str and wildcard_str is True: 
        if debug==True: print("Item is string.")
        # Indicate error if type(item) is string & doesn't match wildcard: 
        if item.lower()!='?all?': 
            is_error=True
            err_reason= 'If "'+item_nm+'" is passed as a STR, only valid input'+\
                        'is "?ALL?".You passed '+item_nm+'="'+item+'".'
            give=item
            if debug==True: print("Item is string and not allowed wildcard.")
        else: 
            if debug==True:  print("Item is string and is an allowed wildcard.")
            is_error=False; err_reason=''; give=allowed_vals
         
    # If input is array like... recognizes lists, np.ndarrays, pd.Series, or xr.dataarrays
    elif ((type(item)==list) or (type(item)==np.ndarray) or (type(item)==pd.core.series.Series)\
        or (type(item)== xr.core.dataarray.DataArray)):
        if debug==True: print('Item is iterable.')
        
        # Check that all elements of array like item are strings, if not throw warning: 
        if not all(isinstance(elem, (str, np.str_)) for elem in item): 
            is_error= True
            err_reason= 'All elements of array like input for "'+ item_nm+\
                '" are not type(str) or type(np.str_), but must be.'
            give=item
            if debug==True: print('Item is iterable, but not all elements are str.')
        else: 
            if debug==True: print('Item is ierable and all are string. ')
            # Check that all elements listed in item are indeed in allowed_vals: 
            if any(elem not in allowed_vals for elem in item):
                if debug==True: print('Item is iterable of strs but has some unallowed values.')
                baddies=','.join([elem for elem in item if elem not in allowed_vals])
                # Warn about unrecognized items and continue with just the recognized ones:
                warnings.warn('The following items in input for "'+ item_nm+
                              '" could not be found in '+allowed_nm+f':\n\t{baddies}.\n'+
                              'These variables will not be requested in created planeflight input files.',
                              UserWarning, stacklevel=2)
                ok_list=[elem for elem in item if elem in allowed_vals]
                is_error=False; err_reason=''; give=ok_list;
            else: 
                if debug==True: print('Item is ierable, all are string, and all allowed. ')
                is_error=False; err_reason=''; give=item
    else: 
        # If "item" is not str or array like, return type error: 
        if debug==True: print('Item was not identified and vals assigned were unknown.')
        is_error= True
        err_reason='Input for "'+item_nm+'" is type='+str(type(item))+', and is not allowed.'
        give=item 
        
    return is_error, err_reason, give


def _read_planelog_to_df(filename:str): 
    """Internal function used to read in a single planedat output file and convert it 
    to a pandas dataframe. Doesn't do any renaming/unit conversions (handled elsewhere)
    but does ensure all columns are converted to appropriate data type. 
    
    INPUTS:
    ------ 
        (1) filename - STRING containing absolute path to an output plane.log file. 
    
    OUTPUTS:
    --------
        (2) df - pandas DATAFRAME with time_UTC as index and data for all outputs as 
                 columns in same units as all data in original file. 
    """
    # Check that input file exists, otherwise throw an error: 
    if not os.path.isfile(filename):
        raise FileNotFoundError('The following input plane.log file passed to'+ 
                                'read_planelog() as "filename" could not be found:\n\t'+
                                 filename).with_traceback(sys.exc_info()[2])
        
    #==========================================================================
    #      Read in header/data and deal with potential data split across lines: 
    #==========================================================================
    # Open file and read in all lines to a list.
    with open(filename, 'r') as fin:
        lines = fin.readlines()
    
    # Make an empty list to hold lists containing data for all headers from each line: 
    ln_list=list([]);
    
    # Loop over each line & seperate data from each header contained on it
    for i in range(0, len(lines)): 
        
        # plane.log is an ASCII delimited file, which doesn't play nice with 
        # most read-in functions. We assume ' ' is delimiter betwen header/data 
        # items, and remove blank spaces/new line chars. So, this extracts only 
        # the actual data from inbetween all those spaces on the line, in a list: 
        ln_i=[s.replace('\n','') for s in lines[i].split(' ') if len(s.replace('\n','')) > 0]
        
        # On first & second lines, define expected odd/even line length: 
        if i==0: odd_len=len(ln_i)
        if i==1: even_len=len(ln_i) 
        
        # If you try to save a bunch of diagnostis and all tracers, output files 
        # my have data/headers for a single time index spread over 2 lines. This affects 
        # how we should read in the file, so decide if lines are split or not after 
        # checking the first 2 lines & comparing their lengths. 
        if i==2: split_lines=True if odd_len!=even_len else False 
    
        # If data is not split across multiple lines are not split, then all lists 
        # of data on a line should be the same length. But, if data for 1 time 
        # index is split across multople lines, then odd # lines should be same 
        # len as line #1 and even # lines should be same len as line #2. So, 
        # here, we check to make sure the list of elements in all lines are the 
        # expected length. If this is not true, we could pull data into wrong 
        # dataframe column or have columns with dif lengths! So, throw an error 
        # if some inconsistency is found. 
        if ((i > 1) and (i%2==0)): 
            if len(ln_i)!=odd_len:
                raise ValueError('Inconsistent length found on line #'+str(i+1)+ 
                                  '\n  Expected len='+str(odd_len)+
                                  '\n    Actual len= '+str(len(ln_i))).with_traceback(sys.exc_info()[2])
        elif ((i>1) and (i%2==1)): 
            if len(ln_i)!=even_len:
                raise ValueError('Inconsistent length found on line #'+str(i+1)+ 
                                  '\n  Expected len='+str(even_len)+
                                  '\n    Actual len= '+str(len(ln_i))).with_traceback(sys.exc_info()[2])
    
        # Keep the list of all items on the line in a nested list with data from all lines:   
        ln_list.append(ln_i)
    
    # Construct the full header/data for each time index so lists come out identical 
    # to make a dataframe regardless of if data was split across multiple lines in files. 
    if split_lines==False:
        # If data is not split across multiple lines, the full header is only on line #1 
        header= ln_list[0]
        
        # And all data for each time index is on the same line. 
        data=ln_list[1:]
        
    else: 
        # If lines are split, header spills over onto line #2 as well, so combine 
        # the first two lines to get the full header listing. 
        header= ln_list[0]+ln_list[1]
        
        # If data is on split lines, then one line of data should be combo of every 2 lines
        data=[ln_list[i]+ln_list[i+1] for i in range(2, len(ln_list),2)]
        
        # Check that we have the same # of elements for all data lines: 
        if len(np.unique([len(data_i) for data_i in data])) >1: 
            raise ValueError('Issue with concatenating data across split lines.'+
                              'At this point, all lines in var, "data" should have same len'+
                              ', but they do not.').with_traceback(sys.exc_info()[2])
    
    # Replace '-' char in headers (pandas doesn't like them in column names) 
    header=[name.replace('-','_') for name in header]
    
    #==========================================================================
    #                    DUMP DATA INTO A PANDAS DATAFRAME: 
    #==========================================================================
    # At this point, data from each line is in the nested list "data"/"headers" 
    # and artifacts from how data is split across lines is all delt with. 
    
    # Create a pandas df from the data with columns matching plane.log headers: 
    df=pd.DataFrame(columns=header)
    
    # Loop over all columns in the df & fill each with an array for all time indexes: 
    for i, col in enumerate(df.columns):
    
        # Pull column(i) from line n and add to array we'll assign to the df column. 
        arr=[data[n][i] for n in range(0,len(data))] # all data are still strings. 
        
        # Convert list of string data in "arr" into a np.array with relevant d-type. 
        if col in ['TYPE', 'YYYYMMDD','HHMM', 'TIME_LT']: 
            # Keep "TYPE" and time columns as string arrays. 
            arr=np.array(arr,dtype=str)
        elif col in ['P_I','I_IND','J_IND','T_IND','GMAO_IIEV','GMAO_JJEV','GMAO_LLEV']:
            # Make index vars ints. But, to convert str to np.int64, 
            # you first have to convert to np.float to not get an error...  
            arr=np.array(arr,dtype=np.float64).astype(np.int64)
        else: 
            # All other vars are fine as np.float64 (just to be safe!). 
            arr=np.array(arr,dtype=np.float64)
        
        # Assign column to array with that specific data type.  
        df[col]=arr 
    
    # Concatenate 'YYYYMMDD' and 'HHMM' columns to datetime object &put in its own df: 
    time_utc = pd.DataFrame({'time_UTC': pd.to_datetime(df['YYYYMMDD']+df['HHMM'], 
                                                        format='%Y%m%d%H%M')})
    
    # Use pd.concat to join this "time_utc" df to the original DataFrame 
    # This method avoids a df fragmentation warning... 
    df = pd.concat([df, time_utc], axis=1)
    
    # Set the 'time_UTC' column as the df's index
    df.set_index('time_UTC', inplace=True)
    
    # Drop the 'YYYYMMDD' and 'HHMM' columns after. 
    df.drop(columns=['YYYYMMDD','HHMM'], inplace=True) 
    
    return df 


def _build_output_meta_dict(df:pd.core.frame.DataFrame, config_yaml:str,spdb_yaml:str ): 
    """Internal function to parse the geoschem_config.yml, species_database.yml, 
    and planeflight_diagnostics.yml files and return a dictionary with meta data
    on species long_name, units, etc. for all columns in output planeflight data. 
    
    INPUTS:
    -------
        (1) df          -  A pandas DATAFRAME with all planeflight outputs.
        
        (2) config_yaml -  A STRING containing path to geoschem_config.yml file. 
        
        (3) spdb_yaml   -  A STRING containing path to species_database.yml file. 
        
    OUTPUTS: 
    ------- 
        (1) info_dict - A nested DICT with 1st key equal to advected species/diagnostic
                        name containing keys 'long_name', 'units' (and 'MW_g' and 'formula' 
                        if these inside species database entry for advected species'
    """
    
    # Read in species database yaml
    spdb= _read_spdb_yml(spdb_yaml)
    
    # Read in all planeflight diags (not including advected species) as dictionary. 
    diags= _read_planeflight_diags_yml(inputs_only=False, no_collections=True)

    # Get dictionary mapping tracer numbers to tracer names...     
    nums2names=_get_tracer_name_num_mapping(config_yaml, keys_are_nums=True)

    # Initialize dictionary to hold info about all vars outputted from planeflight. 
    info_dict=dict({}) 
        
    # Loop over all vars outputted in plane.log:
    for col in df.columns:

        # Skip metadata columns added during concatenation (handled after loop).
        if col in ['PlaneLog_File', 'File_Index']:
            continue

        # Figure out if this is a column containing an advected species
        # concentration, optional diagnostic or is unrecognized...
        if col in list(nums2names.keys()) + list(nums2names.values()):
            is_conc=True;
        elif col in list(diags.keys()):
            is_conc=False;
        else:
            warnings.warn(f"Can't determine if column '{col}' is an advected species "
                          "or optional diagnostic. Skipping metadata for this column.",
                          UserWarning, stacklevel=2)
            continue

        # Initialize dict to hold info about this specific tracer:
        tracer_info=dict({'Long_name':'???','Units':'???'})
        
        # Pull info from species database if its an advected species... 
        if is_conc is True: 
            # Check if column name is advected species tracer NUMBER. 
            if re.match('^TRA_\d{3,4}$', col): 
                
                if col in list(nums2names.keys()):
                    # Lookup species name of this tracer number. 
                    name=nums2names[col]# Name of key in info_dict is species NAME, not tracer #... 
                    
                    # Dump info about var's name, units, formula, molecular weight into dict. 
                    if name in list(spdb.keys()):
                        for key in list(spdb[name].keys()):
                            if key in ['MW_g', 'Formula']:
                                tracer_info[key]=spdb[name][key]
                            elif key =='FullName': 
                                tracer_info['Long_name']='Concentration of '+spdb[name]['FullName']
                    else:
                        warnings.warn(f"Couldn't find info about advected species: {name}",
                                      UserWarning, stacklevel=2)
                else:
                    warnings.warn(f"Entry for tracer number: {col} not found in "
                                  "tracer name/num mapping built from config file.",
                                  UserWarning, stacklevel=2)
                    name=col # Name of key in info_dict is column name/ tracer number...
                        
                # When tracer numbers are used in planeflight, output units of advected 
                # species are are mol mol-1 dry (v/v). So, add that info to dict here: 
                tracer_info['Units']='mol/mol'
                
            elif col in list(nums2names.values()):
                name=col; # Name of key in info_dict is column name.
                
                # Dump info about var's name, units, formula, molecular weight into dict. 
                for key in list(spdb[col].keys()):
                    if key in ['MW_g', 'Formula']:
                        tracer_info[key]=spdb[name][key]
                    elif key =='FullName': 
                        tracer_info['Long_name']='Concentration of '+spdb[name]['FullName']

                # When tracer names are used in planeflight, output units of advected 
                # species are are molec/cm3. So, add that info to dict here: 
                tracer_info['Units']='molec/cm3'
                
        else: 
            name=col; # Name of key in info_dict is column name. 
            
            # Otherwise, pull info out o the optional diagnotistic dictionary: 
            if col in list(diags.keys()): 
                for key in list(diags[col].keys()): 
                    if key =='FULLNAME': 
                        tracer_info['Long_name']=diags[col][key]
                    elif key=='UNITS':
                        tracer_info['Units']=diags[col][key]
            else:
                warnings.warn(f"Couldn't find info about optional diagnostic named: {col}",
                              UserWarning, stacklevel=2)
                    
        # Add info about this tracer to the dictionary: 
        info_dict[name]=tracer_info

    # Add items for indexing when concatenating: 
    info_dict['PlaneLog_File']=dict({'Long_name': 'Original path to raw plane.log file this data came from.'})
    info_dict['File_Index']=dict({'Long_name': 'File index of all plane.log files concatenated this data came from.'})
    return info_dict


def _convert_moleccm3_to_mr(df:pd.core.frame.DataFrame, info_dict:dict):
    """Concentrations in plane.log files are outputted in units of molec cm-3 
    if **TRACER NAMES** rather than **TRACER NUMBERS** were used in the analogous 
    input files. This function uses the Pressure/Temperature outputs to convert 
    *ONLY CONCENTRATION* outputs from molec cm-3 into mol mol-1. Will skip converting
    vars that are not in units of molec cm-3. 
    
    INPUT: 
    ------ 
      (1) df         -  PANDAS DATAFRAME containing plane.log output.
      
      (2) info_dict  -  A nested DICT with 1st key equal to advected species/diagnostic
                        name containing keys 'long_name', 'units' (and 'MW_g' and 'formula' 
                        if these inside species database entry for advected species). 
              
    OUTPUT: 
    ------ 
      (1) df         -  Same as input, but with all species concentration columns
                        now expressed in units of mol/mol instead of molec/cm3.
                        
      (2) info_dict -  Same as input, but with new units that species concentration
                        columns are expressed in, if updated. 
    """
    #==========================================================================
    # Check to make sure Pressure/ Temperature variables exist in the output 
    # (since they're required to convert to mol/mol from molec/cm3):
    #==========================================================================
    if 'PRESS' in df.columns: 
        pres_var='PRESS' # Pressure in hPa for each flight track point
    elif 'GMAO_PRES' in df.columns:
        pres_var='GMAO_PRES' # Pressure at center of grid box in hPa
    else: 
        raise ValueError('No variable for pressure found in your output' +\
                         'files to use in converting concentrations from' +\
                         'molec cm-3 to mol mol-1. Check your output files to see'+\
                         'if "PRESS" or "GMAO_PRESS" are listed in the headers.'
                         ).with_traceback(sys.exc_info()[2]) 
            
    if 'GMAO_TEMP' not in df.columns:#
        raise ValueError('No variable for tempeature found in your output' +\
                         'files to use in converting concentrations from' +\
                         'molec cm-3 to mol mol-1. Check your output files to see'+\
                         'if "GMAO_TEMP" are listed in the headers.'
                         ).with_traceback(sys.exc_info()[2])
            
    #==========================================================================
    # Make Conversion factor for molec/cm3 to mol/mol using temp/pressure in file. 
    #==========================================================================
    # Get all vars in correct units for n= PV/RT to come out in mols: 
    T= df.GMAO_TEMP.values    # Temp in K 
    R=8.314462e4              # Gas Constant in cm3 hPa K-1 mol-1
    P=df[pres_var].values     # Pres in hPa
    Av=6.022e23               # Avagardo's # (molec/mol)

    # Define unit_conversion as a single multiplication factor:
    #                     (      R      T     1/P     1/ Av )
    #   X molec (CMPD)   {   cm3 hPa    K      1       mol  }   Y mols (CMPD)
    #  --------------- * {  ------- * ----*  ----*  ------- } =  -------------
    #    1 cm3 (air)     {    K mol     1     hPa     molec }     mols (air)
    #                                     ^^^
    #                    {       moleccm3_to_molmol        }
    moleccm3_to_molmol=(R*T)/(P*Av)
    
    #==========================================================================
    # Loop over all columns in df. Only convert cols with units of molec/cm3. 
    #==========================================================================
    for col in df.columns: 
        # Get info about the column's curernt units from info_dict: 
        if col not in ['PlaneLog_File','File_Index']: 
 
            col_units=info_dict[col]['Units']
            
            # Only do unit conversion /change stuff if units are molec/cm3.
            if col_units=='molec/cm3':
                
                # Do unit conversion: 
                df[col]=df[col]*moleccm3_to_molmol
                
                # Update info_dict with new units changed too. 
                info_dict[col]['Units']='mol/mol'
                
                # And add a note that this was done... 
                info_dict[col]['Note']="Converted raw output data for this var from molec/cm3"+\
                     "to mol/mol using Temperature=GMAO_TEMP and Pressure={pres_var} "+\
                     "within planeflight_io package."
        
    return df, info_dict


def _convert2_xr(df:pd.core.frame.DataFrame, info_dict:dict):
    """Internal function to convert output df to xarray datset with all 
    metadata about what units/what data is contained within attributes of
    all xarray datarrays inside the input dataset. This function only used 
    if the `output_xarray` option is  set to TRUE in read_planelog() or
    read_and_concat_planelogs(). 
    
    INPUT: 
    ------ 
    (1) df- Pandas dataframe with plane.log output 
          
    (2) info_dict - A nested DICT with 1st key equal to advected species/diagnostic
                    name containing keys 'long_name', 'units' (and 'MW_g' and 'formula' 
                    if these inside species database entry for advected species' 
         
    OUTPUT: 
    ------ 
    (1) ds - Xarray dataset with all columns in df now as data_arrays that each 
             include, at minimum, attributes for "long_name" and "units" (plus others)
             using info from the species database.yml & the planeflight_diagnostics.yaml. 
            
    """
    
    # Ensure 'time_UTC' is index of the df (messes up dims/coords of ds if not)  
    if ((df.index.name!='time_UTC') and ('time_UTC' in df.columns)):
        df = df.set_index('time_UTC')
            
    # Convert the concatenated df to an xarray dataset 
    ds = df.to_xarray() 
    
    # Ensure 'time_UTC' is set as a coordinate.
    ds = ds.set_coords('time_UTC')
    
    # Loop over all vars outputted in plane.log: 
    for var_i in list(ds.data_vars): 
        
        # Get meta data about this specific variable. 
        tracer_info=info_dict[var_i] 
        
        # Add all info from tracer info dict as attributes of this var in the output ds: 
        for key in list(tracer_info.keys()): 
            ds[var_i].attrs[key] = tracer_info[key]
                        
        # Delete "tracer_info" so no info about older tracers gets assigned to wrong var.
        del tracer_info 
    
    return ds


def _get_unique_filename(filepath: str, overwrite: bool, exts: list) -> str:
    """Function to generate a unique file path without an extension, ensuring 
       no existing files with the specified extensions are present in the directory. 
       This function appends an incrementing number to the base file path if
       necessary to avoid conflicts.
    
      INPUTS:
      -------
      (1) filepath : str
          The base file path (without extension) to check for uniqueness. This 
          should be a full path, such as '/dir/to/file'.
    
      (2) overwrite : bool
          Flag to determine behavior when a file already exists. If set to True,
          the function will return the provided filepath as-is, ignoring any 
          existing files. If False, the function will ensure that no files with 
          the specified extensions exist with the same base name.
    
      (3) exts : list
          A list of file extensions (e.g., ['.csv', '.nc']) to check against the
          base filepath. The function seeks to ensure that none of these extensions 
          lead to an existing file so pairs of files will have same base name. 
    
      Returns:
      --------
      (4) new_filepath : str
          A unique file path (without extension) that ensures no file with any 
          of the specified extensions exists. If overwrite is True, the original 
          filepath is returned unchanged.
      """
    if overwrite:
        return filepath
    else: 
        # Start with the original filepath
        new_filepath = filepath; i = 1
    
        # Loop until a unique filename is found
        while any(os.path.exists(f"{new_filepath}{ext}") for ext in exts):
            new_filepath = f"{filepath}_{i}"
            i += 1
    
        return new_filepath
   
    
def _save_outputs(df:pd.core.frame.DataFrame, info_dict:dict, savefile:str, 
                  as_xarr:bool=False, overwrite:bool=False):
    """Internal function to save output planelogs as df/yaml or as xarray ds.""" 
    
    if as_xarr is False: 
        # Get a unique filename pair for df in .csv & info_dict as .yml considering 
        # whether to overwrite any existing files or not... 
        savefile=_get_unique_filename(savefile, overwrite, exts=['.csv','.yaml'])
        
        # Save output pandas dataframe as a csv file: 
        df.to_csv(savefile+'.csv')
        
        # Save metadata dictionary to a YAML file
        with open(savefile+'.yaml', 'w') as file:
            yaml.dump(info_dict, file, default_flow_style=False)
            
        # Print output file paths to screen 
        print('plane.log output file: '+\
              '\n  (1) Data stored in a pandas dataframe saved at: \n\t'+savefile+'.csv'+\
              '\n  (2) Associated metadata in a dict with column names '+\
              'as keys saved at: \n\t'+savefile+'.yaml')
            
        return df, info_dict 
            
    else: # Convert output to xarray dataset instead... 
        
        # Function handles conversion to ds, attaches everything in info_dict as
        # attributes to individual data vars 
        ds=_convert2_xr(df, info_dict)
        
        # Get a unique filename considering whether to overwrite it or not if it already exists. 
        savefile=_get_unique_filename(savefile, overwrite, exts=['.nc'])
        
        # Save output xarray dataset as netcdf file
        ds.to_netcdf(savefile+'.nc') 
        
        # Print  output file path to screen & return ds
        print('Concatenated plane.log data saved as xarray dataset with'+\
              ' metadata attrs at: \n\t'+savefile+'.nc')
        
        return ds 
