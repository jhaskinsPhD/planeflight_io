"""
Script to make readable Planelight.dat.YYYYMMDD input files for GEOS-Chem
and read in outputted planeflight.log files in as pandas dataframes with 
option to concatenate files. 

Created on Sun Mar 21 14:21:14 2021

@author: Dr. Jessica D. Haskins
GitHub: @jhaskinsPhD
Email: jessica.haskins@utah.edu 
"""
import re
import os
import sys
import yaml
import numpy as np
import pandas as pd 
import xarray as xr 
from datetime import datetime


###############################################################################
# INTERNAL HELPER FUNCTIONS- NOT INTENDED TO BE USED OUTSIDE MAIN FUNCTIONS
###############################################################################
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
                          
        (2) no_collectrions - BOOL indicating if you want to return diags organized by 
                        "collection" or not. Default is FALSE (to include them 
                        arraged by collection) in a nested dict format where first 
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
        # Standard/defulat output has all diags arranged by collections. return that. 
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

        lines.append('Diagnostiscs:')
        for key in list(diags[coll]['Diagnostics'].keys()):
            lines.append(f"\t{key:<10} = {diags[coll]['Diagnostics'][key]['FULLNAME']}")
        if 'Notes' in list(diags[coll].keys()): 
            lines.append(f"Notes:\n\t{diags[coll]['Notes']}")
    
    # Join all lines together with line breaks & output single string:     
    info= '\n'.join(lines)
    
    return info  


def _check_str_arr_inputs(item, item_nm:str, allowed_vals, allowed_nm:str,
                          wildcard_str:bool=False ): 
    """ Internal helper function to check that input for "item" is EITHER a 
    string wildcard  matching '?ALL?' (if wildcard_str) is set to TRUE, OR is
    array like, with all array elements as type(string), and with all elements
    are within "allowed_vars". Returns if an error is found & why the error happened,
    Function only called within make_planeflight_inputs() to help with user 
    input validation & has not other valid usage. 
    
    CHANGE LOG: 
    ------------
    10/01/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
    """ 
     
    if type(item)==str and wildcard_str is True: 
        # Indicate error if type(item) is string & doesn't match wildcard: 
        if item.lower()!='?all?': 
            is_error=True
            err_reason= 'If "'+item_nm+'" is passed as a STR, only valid input'+\
                        'is "?ALL?".You passed '+item_nm+'="'+item+'".' 
            return is_error, err_reason, item
         
    # If input is array like... recognizes lists, np.ndarrays, pd.Series, or xr.dataarrays
    elif ((type(item)==list) or (type(item)==np.ndarray) or (type(item)==pd.core.series.Series)\
        or (type(item)== xr.core.dataarray.DataArray)):
        
        # Check that all elements of array like item are strings, if not throw warning: 
        if not all(isinstance(elem, (str, np.str_)) for elem in item): 
            is_error= True
            err_reason= 'All elements of array like input for "'+ item_nm+\
                '" are not type(str) or type(np.str_), but must be.'
            return is_error, err_reason, item
        else: 
            # Check that all elements listed in item are indeed in allowed_vals: 
            if any(elem not in allowed_vals for elem in item):
                
                # If any are found that aren't allowed, print & ask what to do next: 
                print('The following items in input for "'+ item_nm+\
                '" could not be found in '+allowed_nm+':\n') 
                [print(elem) for elem in item if elem not in allowed_vals]
                
                while True: 
                    user_input = input('Enter "1" to drop these items from the input list, "'+
                                       item_nm+'" & proceed.\n'+'Enter "2" to exit now.')
                    if user_input ==1: 
                        # Revise list to only include entries in allowed_vars
                        out_list=[elem for elem in item if elem in allowed_vals]
                        break
                    elif user_input == 2: 
                        exit(0)  
                    else: 
                        print('Invalid input. Please enter 1 to drop these items from the'+
                             'tracer list & proceed or enter 2 to exit.')
                return False, '', out_list      
    else: 
        # If "item" is not str or array like, return type error: 
        is_error= True
        err_reason='Input for "'+item_nm+'" is type='+type(item)+', and is not allowed.'
        return is_error, err_reason, item


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
    fin = open(filename, 'r')  
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
        
        # Figure out if this is a column containing an advected species 
        # concentration, optional diagnostic or is unrecognized...  
        if col in list(nums2names.keys()) + list(nums2names.values()): 
            is_conc=True;
        elif col in list(diags.keys()):
            is_conc=False;
        else: 
            if col not in ['PlaneLog_File','File_Index']: 
                raise Warning(f"Cant' determine if column {col} is an advected species" +\
                              "or optional diagnostic!!").with_traceback(sys.exc_info()[2]) 
            
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
                        for key in list(spdb[col].keys()):
                            if key in ['MW_g', 'Formula']:
                                tracer_info[key]=spdb[name][key]
                            elif key =='FullName': 
                                tracer_info['Long_name']='Concentration of '+spdb[name]['FullName']
                    else: 
                        raise Warning(f"Couldn't find info about advected species: {name}"\
                                      ).with_traceback(sys.exc_info()[2]) 
                else: 
                    raise Warning(f"Entry for tracer number: {col} not found in"+\
                                  "tracer name/num mapping built from config file."\
                                  ).with_traceback(sys.exc_info()[2]) 
                        
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
                raise Warning("Couldn't find info about column / optional diagnotic"+\
                              f"named:{col}" ).with_traceback(sys.exc_info()[2]) 
                    
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
        savefile=_get_unique_filename(savefile, overwrite, exts=['.csv','.yml'])
        
        # Save output pandas dataframe as a csv file: 
        df.to_csv(savefile+'.csv')
        
        # Save metadata dictionary to a YAML file
        with open(savefile+'.yaml', 'w') as file:
            yaml.dump(info_dict, file, default_flow_style=False)
            
        # Print output file paths to screen 
        print('plane.log output file: '+\
              '\n  (1) Data stored in a pandas dataframe saved at: \n\t'+savefile+'.csv'+\
              '\n  (2) Associated metadata in a dict with column names '+\
              'as keys saved at: \n\t'+savefile+'.yml') 
            
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
              'metadata attrs at: \n\t'+savefile+'.nc')
        
        return ds 
    
###############################################################################
# MAIN FUNCTIONS- INTENDED TO BE USED 
###############################################################################
def get_compatible_input_diags(simtype:str='',these_collections:list=[],
                                   display:bool=False): 
    """Function to load all optional planeflight diagnotics in an 
    output dicitonary as keys with info on diagnostic's FULLNAME/UNITS stored 
    in values that are compatiable with a user's simulation type & print results
    to screen if requested. 
    
    INPUT: 
    -------
    
    (1) simtype           - STRING containing user's simulation type. 
                            Used to determine what available planflight diagnostics are 
                            compatiable with your simulation type.
                    
    (2) these_collections - (OPTIONAL): LIST of STRINGS containing diagnostic 
                            collection names that you would like to retrieve. 
                            Useful if you only want to get specific collections 
                            of diags, not just ALL diags compatiable with simtype: 
                            Valid options to pass include: 
                              ['aer_uptake', 'aodb', 'aodc', 'aq_aer', 'broken', 
                               'chem_fams', 'defaults', 'gmao_ice', 'gmao_met', 
                               'hg', 'htep', 'tomas']
                    
    (2) display           - (OPTIONAL): BOOL indicating if compatible requested 
                            diagnotics should be also be printed to screen.     
    OUTPUT: 
    -------
    
    (1) diag_dict - DICTIONARY with keys corresponding to diagnostic names
                    used in planeflight.dat input files / plane.log output 
                    files with VALUES being a dictionary with "FULLNAME" and 
                    "UNITS" keys describing what the diagnostic is / its units.  
    """ 
    # Check that user passed an allowed/recognized simtype or throw warning. 
    _check_simtype(simtype, blank_allowed=False)

    # Read in all valid input planeflight diags (not including advected species) as dictionary. 
    diags= _read_planeflight_diags_yml(inputs_only=True)
        
    # Initialize dict holding strings of all available diagnostic collections: 
    ok_diags=list(diags.keys()) 
    
    ###########################################################################
    # Remove invalid collections from "ok_diags" based on user's simulation type
    ###########################################################################
    if simtype.lower()!='hg':    ok_diags.remove('hg') 
    if simtype.lower()!='tomas': ok_diags.remove('tomas')
    if simtype.lower()=='carbon': 
        for item in ['aer_uptake','htep','chem_fams']: ok_diags.remove(item) 
    
    # The "broken" collection includes diags that do not work with any simtype b/c
    # they need source code updates in GEOS-Chem... Only return if "all" is requested. 
    ok_diags.remove('broken') 

    ###########################################################################
    # If user asked for specific compatible diags by listing collections: 
    ###########################################################################
    if len(these_collections) > 0: 
        
        # Check that user only passed valid args for "collections" or throw error:   
        invalid= [item for item in these_collections if item.lower() not in list(diags.keys())]
        if len(invalid)>0: 
            raise ValueError('Invalid planeflight diagnostic collection name(s):\n\t'+
                             ','.join(invalid)+'\nValid collection names include: \n\t'+
                             ','.join(list(diags.keys()))).with_traceback(sys.exc_info()[2])
                        
        # Check for any requested collections that aren't compatiable with simtype: 
        not_compat=[item for item in these_collections if item not in ok_diags]
        if len(not_compat) >0:
            # If all requested collections aren't copmatiable, throw error: 
            if len(not_compat)==len(these_collections): 
                info=_display_diags(diags, ok_diags)
                raise ValueError('None of the requested collections are compatiable with your '+ 
                                 'simulation type. The following options are those compatiable:'+info )
                
            else: 
                # Otherwise just throw a warning we're dropping some... 
                raise Warning('The following requested collections are not compatiable '+ 
                              'with your simulation type and will NOT be included in the '+
                              'outputted diagnostic collection:\n\t'+','.join(not_compat)+
                              '\nUse pln.get_compatible_input_diags(simtype,display=True) to display '+
                              'all diagnostics compatiable with your simulation type.'
                              ).with_traceback(sys.exc_info()[2]) 
            
        # Define collections so it only includes compatiable collections: 
        collections=[item for item in these_collections if item in ok_diags]
            
    else: # If specific collections weren't requested, return all compatiable options: 
        collections=ok_diags 
    
    ###########################################################################
    # Parse the now curated collections list and build output diags dict: 
    ###########################################################################
    diag_dict=dict({})
    for c in collections: 
        for key in diags[c]['Diagnostics']: 
            diag_dict[key]=diags[c]['Diagnostics'][key]
    
    ###########################################################################
    # Print compatiable diags to screen if asked 
    ###########################################################################
    if display is True: 
        # Format info about diags in parse_this to print to screen: 
        info=_display_diags(diags, collections)
        print(info)
        
    return diag_dict


def make_planeflight_inputs(savedir: str, 
                            gc_config:str,
                            datetimes,
                            lat_arr, 
                            lon_arr, 
                            vert_arr,
                            vert_is_pres:bool, 
                            tracers,
                            diags,
                            diags_minus: list =[],
                            username: str = 'user',
                            overwrite: bool = False,
                            drop_dupes: bool = False):
    """Function to create planeflight.dat files in correct format for GEOS-Chem input.
    
    INPUT:
    -------
    Note: For all inputs where "ARRAY" is accepted, vars can be 1-D lists, np.ndarray,
          pd.Series, or xr.dataarray. 
    
       (1) savedir      - STRING containing path to directory in which to save
                          the output planeflight.dat files at. 
                      
       (2) gc_config   - STRING containing path to your geoschem_config.yml file 
                          used to read in/define valid names of transported tracers
                          and to determine simulation type, used to determine all valid 
                          additional optional planeflight diagnostics.
       
       (2) datetimes    - PANDAS SERIES of datetimes in UTC stored as Timestamps at which 
                          to sample the model.Type(datetimes[0]) should return: 
                          <class 'pandas._libs.tslibs.timestamps.Timestamp'>
                          
                          To create input in the correct format do: 
                               date_range = pd.date_range(start='2017-01-01', end='2017-01-03', freq='60s')
                               datetimes=pd.Series(date_range) 
       
       (3) lat_arr      - ARRAY of latitudes at which to sample the model. (range: -90 to 90 deg)
       
       (4) lon_arr      - ARRAY of longitudes at which to sample the model (range: -180 to 180 deg)
    
       (5) vert_arr     - EITHER an ARRAY of pressures (hPa) OR altitudes above 
                          the ground (meters) at which to sample the model. See 
                          https://github.com/geoschem/geos-chem/issues/320 
                          for discussion on whether input altitudes should be 
                          "above ground" or "above sea level".     
                         
       (6) vert_is_pres - BOOLEAN indicating if "vert_arr" containined pressures or not. 
                          When TRUE,  values are assumed to be pressures (hPa). 
                          When FALSE, values are assumed to be altitudes (meters).       
                          
       (7) tracers     -  Either (1) an ARRAY of specific advected tracers you want to sample 
                              OR (2) a STRING equal to '?ALL?' to sample all advected species 
                          listed in your geoschem_config.yml file. 
                          
       (8) diags       -  (OPTIONAL) Either (1) an ARRAY containing STRINGS of any additional 
                          diagnostics to sample from model (beyond tracers) OR (2) a 
                          STRING equal to '?ALL?' to sample all available additional 
                          diagnostics compatiable with your simulation type. Default is 
                          to include the grid-box indexes planflight pulled from & 
                          Pres/Temp/RH at center of grid box (e.g. ['GMAO_IIEV',
                          'GMAO_JJEV', 'GMAO_LLEV', 'GMAO_PRES','GMAO_RELH', 'GMAO_TEMP']).  
                          
       (9) diags_minus - (OPTIONAL) ARRAY containing STRINGS with all additional 
                         diagnostics you don't want to include (only relevant 
                         if you passed diags='?ALL?'). 
                          
       (10) username   - (OPTIONAL) STRING contaiing name of user who created files. 
                         This gets listed in header of resulting planedat input files. 
                
       (11) overwrite  - (OPTIONAL) BOOLEAN of whether to overwrite any existing files 
                          at 'savedir' with this name or not. If FALSE, & any files 
                          under 'savedir' would be overwritten, a new sub-directory 
                          under 'savedir' named "NEW_YYYYMMDD_HHMMSS" is created to 
                          hold all the new output files. If TRUE, only files
                          with conflicting names under 'savedir' are overwritten.
                          Default is set to FALSE (not to overwrite files). 
                                 
    OUTPUT:
    ------
       (1) One file for each day listed in 'datetimes' named "planeflight.dat.YYYYMMDD"
           written to "savedir"that can be passed to GEOS-Chem as input files for 
           the planeflight diagnostic. 
           
    """
    # -------------------------------------------------------------------------
    #                   Perform User Input Checks: 
    #--------------------------------------------------------------------------
    
    # Check that the output "savedir" directory exists: 
    if not os.path.isdir(savedir):
        raise OSError("The following directory passed for savedir='"+savedir+
                      "' does not exist.").with_traceback(sys.exc_info()[2])
    
    # Check that datetimes, lat_arr, lon_arr, and vert_arr don't contain NaNs/NaT/Inf: 
    if (pd.isna(datetimes)).any(): 
        raise ValueError('Input "datetimes" contains np.NaT or np.nan.').with_traceback(sys.exc_info()[2]) 
    if ((np.isnan(np.array(lat_arr)).any()) or (np.isinf(np.array(lat_arr)).any())): 
        raise ValueError('Input "lat_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])
    if ((np.isnan(np.array(lon_arr)).any()) or (np.isinf(np.array(lon_arr)).any())): 
        raise ValueError('Input "lon_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])
    if ((np.isnan(np.array(lat_arr)).any()) or (np.isinf(np.array(lat_arr)).any())): 
        raise ValueError('Input "vert_arr" contains np.nan or np.inf.').with_traceback(sys.exc_info()[2])
        
    # Check that all lats are in range -90 to 90 or raise error: 
    if any(abs(lat_arr)> 90):
        raise ValueError('Some |Latitudes| are > 90 (e.g. out of range)'+ 
                         'in input "lat_arr".').with_traceback(sys.exc_info()[2])
        
    # Check that all lons are in range -180 to 180 or raise error:  
    if any(abs(lon_arr)> 180):
        raise ValueError('Some |Longitudes| are > 180 (e.g. out of range)'+
                 'in input "lon_arr".').with_traceback(sys.exc_info()[2])
                                                        
    # Check that pressures/alts in vert_arr are not < 0 or raise error:  
    if any(vert_arr< 0):                                                
        raise ValueError('Some values in "vert_arr" are < 0 (e.g. out of range).'
                         ).with_traceback(sys.exc_info()[2])     
    
    # Check that no pressures are > 1100 hPa or ask if they want to proceed:                                                                         
    if ((vert_is_pres is True) and (any(vert_arr >1100))):
        print('WARNING: Some |Pressures| in "vert_arr" are > 1100 hPa.'+
              'Did you input them in a unit other than hPa by mistake?') 
        while True: 
            user_input = input('Enter "1" to proceed if this was intended.\n'+
                               'Enter "2" to exit now.')
            if user_input == 2: 
                exit(0)  
            elif user_input!=1: 
                print("Invalid input. Please enter '1' to proceed or '2' to exit.")
                
    # Check that NOT all alts <15 m, ask them if they want to proceed or not. 
    elif ((vert_is_pres is False) and (all(vert_arr< 15))):
         print('WARNING: All altitude inputs in "vert_arr" are < 15 m.'+
               'Did you input them in kilometers by mistake?') 
         user_input = input('Enter "1" to proceed if this was intended.\n'+ 
                            'Enter "2" to exit now.')
         if user_input == 2: 
             exit(0)  
         elif user_input!=1: 
             print("Invalid input. Please enter '1' to proceed or '2' to exit.")
                             
    # Check that 'lat_arr', 'lon_arr', 'vert_arr', and 'datetimes' are all the same length: 
    input_arr_lens=[len(var) for var in [datetimes, lat_arr, lon_arr, vert_arr]]
    if any([len_i != input_arr_lens[0] for len_i in input_arr_lens]): 
        raise ValueError('Inputs for "datetimes","lat_arr","long_arr" and "vert_arr" '+ 
                         'must all be the same length.'+
                         '\n\tlen(datetimes)='+str(input_arr_lens[0])+
                         '\n\tlen(lat_arr)='+str(input_arr_lens[1])+ 
                         '\n\tlen(long_arr)='+str(input_arr_lens[2])+ 
                         '\n\tlen(vert_arr)='+str(input_arr_lens[3])).with_traceback(sys.exc_info()[2])
        
    # Check that input geoschem_config.yml file exists & extract list of all 
    # advected species and the simulation type from it if it does (otherwise error) 
    adv_species, simtype=_parse_gc_config(gc_config)
    
    # Check that input for "tracers" is either our wild card string ('?ALL?') or 
    # is a list-like type of tracers that doesn't contain any tracers not tracked in model: 
    tracer_err, tracer_err_msg, tracers=_check_str_arr_inputs(tracers,'tracers',adv_species, 
                                                     'the list of advected species in geoschem_config.yml',
                                                     wildcard_str=True)
    
    # If any checks on "tracers" input were failed, throw descriptive error
    if tracer_err is True: 
        raise ValueError('Invalid input for "tracers".'+tracer_err_msg+
                         ' Input for "tracers" must either be\n:'+ 
                         '\t(1) A STRING equal to "?ALL?" to sample all'+
                         'advected species listed in geoschem_config.yml'+
                         '(e.g. tracers="?ALL/").\n''\t(2) A LIST/ARRAY of '+
                         'STRINGS containing tracers you wish to sample '+
                         '(e.g. tracers=["NO","NO2"]').with_traceback(sys.exc_info()[2]) 
            
    if len(diags)>0: # If the has user asked to include specific optional diagnostics... 
   
        # Get dictionary of ALL optional diagnostics compatiable with user's simtype. 
        diag_dict=get_compatible_input_diags(simtype, display=False)
        
        # Check that input for requested "diags" is either our wild card string ('?ALL?') or 
        # is a list/array of diagnostics including ONLY those compatiable with user's simulation type.   
        diag_err, diag_err_msg, diags =_check_str_arr_inputs(diags,'diags',list(diag_dict.keys()), 
                                                           'the list of optional diagnostics compatible'+
                                                           'with simtype="'+simtype+'" outputted from '+
                                                           'get_diag_info(simtype,return_all=False)',
                                                           wildcard_str=True) 
                                                            
        # If any checks on "diag" were failed, throw descriptive error
        if diag_err is True: 
            raise ValueError('Invalid input for "diags".'+diag_err_msg+
                              ' Input for "diags" must either be\n:\t(1) A STRING '+ 
                              'equal to "?ALL?" to include all optional diagnostics '+ 
                              'compatible with your simtype="'+simtype+'"(e.g. dias="?ALL/").\n'+
                              '\t(2) A LIST/ARRAY of STRINGS containing diags you wish '+ 
                              'to sample (e.g. diags=["GMAO_RELH","RO2"].\n To see a list of '+
                              'all optional diagnostics compatible with your simulation ' + 
                              'try: pln.get_diag_info("'+simtype+'",print_options=True)'\
                              ).with_traceback(sys.exc_info()[2]) 
        
        if type(diags)==str: 
            # If input is string (which will match wildcard if made it this far),
            # then grab all valid options from diag_dict to define list of diags to include. 
            diags=[*diag_dict]            
        
        if len(diags_minus) > 0: # If the has user asked to exclude specific optional diagnostics... 
        
            # Check that input for "diags_minus" is a list/array of diagnostics
            # and that it only includes diagnositics we recognize from list of diags: 
            diag_minus_err, diag_minus_err_msg, diags_minus =_check_str_arr_inputs(diags_minus,'diags_minus',
                                                               list(diag_dict.keys()), 
                                                               'the list of optional diagnostics compatible'+
                                                               'with simtype="'+simtype+'" outputted from '+
                                                               'get_diag_info(simtype,baddies=False)',
                                                                wildcard_str=False) 
            
            # If any checks on "diags_minus" were failed, throw descriptive error
            if diag_minus_err is True: 
                raise ValueError('Invalid input for "diags_minus".'+diag_minus_err_msg+
                                  ' Input for "diags_minus" must be a LIST/ARRAY '+ 
                                  'of STRINGS containing diags you do NOT wish '+ 
                                  'to sample (e.g. diags_minus=["GMAO_RELH","RO2"].\n To see a list of '+
                                  'all optional diagnostics compatible with your simulation ' + 
                                  'try: pln.get_diag_info("'+simtype+'",print_options=True)'\
                                  ).with_traceback(sys.exc_info()[2]) 
                    
            # Remove any diagnostics in list "diag_minus" if they apprear in "diags":    
            diags = [item for item in diags if item not in diags_minus] 
        
        # -------------------------------------------------------------------------
        # Assign "TYPE" column in output planedat files appropriately for vert inputs: 
        #--------------------------------------------------------------------------
        # GEOSCore/planeflight_mod.F90 decides input is pres/alt based on the string 
        # assigned as "TYPE". But, it only converts alts to pressures to pull sample from 
        # model if one of the following strings are passed as the "TYPE": 
        valid_alt_types=['Aacg','Aabne','Acar','Acma','Acrv','Adnd','Aesp','Aetl',
                     'Ahil','Ahip','Alef','Anha','Apfa','Arta','Asca','Asgp','Atgc','Athd','Awbi',
                     'Tbao','Tamt', 'Tcrb','Tlef', 'Tlew','Tmbo','Tmvy',
                     'Tmwo','Tnwr', 'Tsct', 'Tsgp','Tstr','Twbi','Twgc','Twkt']
         
        if vert_is_pres==True: 
            # If user indicates values in vert_arr are PRESSURE, assign TYPE as "Pinp"
            # which is NOT one of the strings that will trigger GEOS-Chem to assume input 
            # is altiutude and convert it to pressure: 
            type2assign='Pinp'
        else: 
            # If user indicates values in vert_arr are ALTITUDES, assign TYPE to be 
            # one of the strings that triggers GEOS-Chem to convert alts to pressures: 
            type2assign=valid_alt_types[0] 
        
        # Make list of all tracers/diagnostics that planedat files should include: 
        tracer_list= diags+ tracers
        ntracers=str(int(len(tracer_list))) # count the number of quantities. 
        
        # Designate a few variables that we'll use to make the header lines of the files
        today= str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) 
        spacer='-------------------------------------------------------------------------------'
        title = '  Now give the times and locations of the flight' 
        
        # Parse the passed dates, get list of individual days we need planeflight files for:  
        all_dates= datetimes.dt.date
        unq_dates=np.unique(datetimes.dt.date)
        
    # Make list of filenames we'll create & check to see if any existing files 
    # under savedir match that pattern/ would be overwritten: 
    new_files=['Planeflight.dat.'+unq_dates[i].strftime('%Y%m%d') for i in range(0, len(unq_dates))] 
    old_files = [f for f in os.listdir(savedir) if os.path.isfile(os.path.join(savedir, f)) and 'Planeflight.dat' in f]  
    dupe_files=[file for file in new_files if file in old_files] 
    
    # Create a new subdirectory under savedir to hold files if overwrite=False. 
    if len(dupe_files) > 0 and overwrite==False: 
        print('Planeflight.dat files found that would be overwritten under: '+ savedir)
        new_savedir=os.path.join(savedir, 'NEW_'+today.replace('-','').replace(':',''))
        print('Creating new directory to hold new planeflight.dat files at:'+new_savedir)
        os.mkdir(new_savedir)
        savedir=new_savedir
        
    # Loop over each unique day in date range we need to make a planeflight.dat file for: 
    for i in range(0, len(unq_dates)): 
        
        # Create a filename based on the date: 
        filename='Planeflight.dat.'+ unq_dates[i].strftime('%Y%m%d')
        
        if overwrite==True: # Delete any existing files under savedir with this name.
            if os.path.isfile(os.path.join(savedir,filename)):
                os.remove(os.path.join(savedir,filename))
                  
        # Get indexes in larger array of flight obs that took place on this date.
        inds = np.where(all_dates== unq_dates[i])[0] 
        
        if len(inds) > 0: 
            #==========================================================================
            # Build arrays of all the data we need in Planeflight.dat for this date.
            # Rount to # of decmial places allowed in GEOS-Chem input for these vars. 
            #==========================================================================
            
            # Define array with that will fill POINTS column (indexing var) 
            points = np.arange(1, len(inds)+1).astype(str) 
            # Define array holding data to fill OBS column full of 9999.000 
            obs = np.full(len(inds),9999.000)              
            # Define/Format array that will hold TYPE col (allowex max of 7 chars)
            typez = np.full(len(inds),'{typ: >6}'.format(typ=type2assign)) 
            # Define array holding DATE (GMT) for each obs point in DD-MM-YYYY format
            day = datetimes.dt.strftime('%d-%m-%Y')[inds]  
            # Define array holding TIME (GMT) for each obs point in HH:MM format: 
            tms = datetimes.dt.strftime('%H:%M')[inds]     
            # Define array holding LAT, LON, and ALT/PRE info formatted to max 2 decimals: 
            lats = np.around(lat_arr[inds], decimals=2)    
            lons = np.around(lon_arr[inds], decimals=2)
            verts= np.around(vert_arr[inds], decimals=2) #GC doesn't allow more than 2 decimals
            
            # Pack all this into a dataframe, because pandas writes tab delimited files nicely! 
            df= pd.DataFrame({'POINT': points, 'TYPE': typez, 'DD-MM-YYYY': day, 
                              'HH:MM':tms, 'LAT':lats, 'LON':lons, 'ALT/PRE': verts, 
                              'OBS':obs})
            
            # Append a line at the bottom of the DF that says its the end! 
            new_row = pd.DataFrame([{'POINT': 99999, 'TYPE': '{typ: >6}'.format(typ='END'), 
                                     'DD-MM-YYYY': '00-00-0000', 'HH:MM': '00:00', 
                                     'LAT': 0.00, 'LON': 0.00, 'ALT/PRE': 0.00, 
                                     'OBS': 0.000}])
            
            # Use pd.concat to append the new row to the existing DataFrame
            df = pd.concat([df, new_row], ignore_index=True)
            
            #=========================================================================
            # Format everything so they're in GEOS-Chem's expected format!
            #==========================================================================
            
            # Note: I got these vals used for each var from planeflight_mod.f90 i/o checks.
            df.POINT=df.POINT.map('{: >5}'.format) # POINTS has len 5 char strings 
            df.LAT=df.LAT.map('{:7.2f}'.format)    # LAT has max 6 total digits, w/ 2 decimals    
            df.LON=df.LON.map('{:7.2f}'.format)    # LON has max 6 total digits, w/ 2 decimals  
            df['ALT/PRE']=df['ALT/PRE'].map('{:7.2f}'.format)# ALT/PRE has max 6 total digits, w/ 2 decimals  
            df.OBS=df.OBS.map('{:10.3f}'.format)   # OBS has max 10 total digits, w/ 3 decimals  
            
            # Write the header lines that GEOS Chem expects to the file, considering format.
            header= '{strr: >6}'.format(strr='POINT ')+ \
                '{strr: >7}'.format(strr='TYPE ')+ \
                '{strr: >11}'.format(strr='DD-MM-YYYY ')+ \
                '{strr: >6}'.format(strr='HH:MM ')+\
                '{strr: >8}'.format(strr='LAT ')+\
                '{strr: >8}'.format(strr='LON ')+\
                '{strr: >7}'.format(strr='ALT/PRE') +\
                '{strr: >11}'.format(strr='OBS')
            
            #======================================================================
            #            Begin writing the planeflight.dat text file
            #=====================================================================
            # This is just a list containing our header lines in the right order... 
            textList = [filename, username, today, spacer, ntracers, spacer] + \
                        tracer_list + [spacer, title, spacer, header]
            
            # Open the output file and write headers line by line. 
            outF = open(savedir+filename, "w")
            for line in textList:
                outF.write(line)
                outF.write("\n")
                
            # Write the data to a temporary file, using ASCII encoding! Pandas default 
            # in Python 3 uses UTF-8 encoding, which GEOS-Chem can't read.
            df.to_csv(savedir+filename+'_0', header=False, index=None, sep=' ', mode='a',
                      encoding='ascii') 
            
            # Annoyingly b/c we use a space as a delimiter in our resulting csv file,
            # we get quotation marks around strings, so open the temporary file, read line 
            # by line, and take out the quotation marks, and write that to the 
            # actual output file and then delete the temporary file: 
            fin = open(savedir+filename+'_0', 'r')
            Lines = fin.readlines()
            for line in Lines:
                outF.write(line.replace('"','' ))
                
            outF.close() # Close the output file. 
            fin.close() # Close the tempororay file 
            os.remove(savedir+filename+'_0') # And delete the temp file. 
            
            print('Output saved at: '+ savedir + filename) # Tell where output is saved.
    
    return


def read_planelog(planelog_file: str, spdb_yaml:str, config_yaml:str, 
                 as_xarr:bool=False, convert2_molmol: bool = False, 
                 output_dir:str='', output_file:str='', overwrite:bool=False):
    """Function to read a single plane.log ouput file into either a pandas dataframe
    or xarray dataset. Handles data split across multiple lines in a plane.log file 
    if a large number of output diagnostics are requested and auto converts 
    the "YYYYMMDD" and "HHMM" date/time diagnostics into a pandas datetime object
    used as either the pandas dataframe index or the xarray coordinate. Automatically 
    corrects output units to mol/mol if tracer names rather than tracer numbers used
    if asked. Will always output metadata telling you about the vars in the outputs. 
    
    INPUT: 
    ------
        (1) planelog_file  -  STRING containing absolute path to a plane.log file. 
        
        (2) spdb_yaml      -  STRING containing absolute path to the species_database.yml 
                              file associated with your run. Used to add metadata.
                             
        (3) config_yaml    -  STRING containing path to geoschem_config.yml file 
                              associated with your run. Used to add metadata.
        
        (4) as_xarr         - (OPTIONAL) BOOLEAN indicating if you want to output the data as 
                              an xarray dataset instead of as a pandas dataframe. 
                              If this is set to TRUE, the data arrays within the 
                              output dataset will have attributes containing metadata
                              assigned to them. Otherwise, output in the xarray 
                              dataset is the same as that in the df. Default is set 
                              to FALSE (to return a pandas df instead). 
                         
        (5) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all 
                              outputted concentrations into units of mol/mol. 
                              Only variables with concentration units are 
                              converted & others are left in their native unit. 
                              Default is set to FALSE (not to do conversion).
                              If set to "TRUE" and "TRA_###" appears in plane.log headers, 
                              no conversion is actually performed to prevent against 
                              improper usage. 
                              
        (6) output_dir       - STRING containing absolute path to a directory where 
                              output file will be saved. If none is passed, then 
                              output saved in same dir as input planelog file.
                              
        (7) output_file     - STRING containing desired name of ouptput file. If none
                              is passed, the output file will be named: 
                                  "planelog_[startYYYYMMDD]_[endYYYYMMDD]"
                                  
        (8) overwrite       - BOOL indicating that if there is already a file named 
                              "output_file" in the output directory, "out_dir", whether 
                              that file should be overwritten or if a new unique name 
                              should be used. 
                            
        NOTE: Outputted advected species concentrations are in molec/cm3 if Tracer 
              Names are used, but are in mol/mol if Tracer Numbers are used in input files. 
              For more details on why this happens, see GitHub issue #796: 
                  https://github.com/geoschem/geos-chem/issues/796    

    OUTPUT: 
    -------
    If output_xarr==FALSE (Default): 
        
        (1) df        - PANDAS DATAFRAME indexed along a pd.datetime object index named 'time_UTC'. 
                         Saved at: output_dir+output_file+'.csv' 
                 
        (2) info_dict - DICTIONARY containing info about output vars & units.          
                        Saved as yaml at: output_dir+output_file+'.yaml'
                 
    If output_xarr == TRUE: 
        
        (1) ds - XARRAY DATASET indexed along a np.datetime64 type coordinate 
                 & dimension named 'time_UTC'. All data vars within the ds have 
                 individual attributes for "long_name" and "units" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. If this is done an additional 
                 attribute 'note' is included indicating this occurred.
                 Saved at: output_dir+output_file+'.nc'
            
    Regardless of the output type selected, the resulting df/ds contains columns/data arrays 
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are 
    dropped after constructing the time index/coordinate. Column/data_var names are 
    nearly identical to the diagnostics in plane.log, with the exception that any 
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I').        
    """
    
    # Read raw file in as dataframe (no unit conversions / renaming)  
    df = _read_planelog_to_df(planelog_file)
    
    # Build meta data dictionary with info about all vars in output dataframe. 
    #   This Assumes advected species units are mol/mol if tracer numbers used, but 
    #    molec/cm3 if tracer numbers used. 
    #
    # **NOTE**: This line MUST remain before renaming df columns with tracer numbers 
    #           to tracer names for output units to be correct in info_dict since 
    #           assigning them properly requires us to know if the column names 
    #           were originally tracer numbers or tracer names. 
    info_dict=_build_output_meta_dict(df, config_yaml,spdb_yaml)
    
    # Get dictionary to map tracer numbers to tracer names: 
    nums2names=_get_tracer_name_num_mapping(config_yaml, keys_are_nums=True)
    
    # Rename all columns in dataframe so columns with tracer numbers are converted 
    # to tracer names. Won't affect anything if cols are already named with tracer names. 
    df.rename(columns=nums2names, inplace=True)
     
    # Convert concentrations to mol/mol from molec/cm3 & update units in info_dict
    # accordingly. Also add note if conversion performed or not ... 
    if convert2_molmol is True: 
        # NOTE: Function will not actually do anything if conc_units!='molec/cm3' 
        # which is the case if tracer numbers were used in analogous input files. 
        df, info_dict = _convert_moleccm3_to_mr(df, info_dict)
    
    # =========================================================================
    # Configure output filename/directory: 
    # =========================================================================
    # Verify that output dirpath exists, otherwise assign output_dir to dir of input file. 
    if not os.path.isdir(output_dir):
        raise Warning('The output directory passed to read_planelog() could not '+ 
                      'be found:\n\t'+output_dir+'/n Output will be saved at '+
                      'input planelog file directory instead.')
        output_dir= os.path.dirname(planelog_file)
        
    if len(output_file)==0: 
        # Set to name of input file witout extension... 
        output_file = os.path.splitext(os.path.basename(planelog_file))[0]
    
    # Combine output directory and output filename for abs path to output files
    savefile= os.path.join(output_dir,output_file)
    
    # Save output / convert to xarray if asked & print to screen where outpput saved: 
    out = _save_outputs(df, info_dict,savefile, as_xarr=as_xarr, overwrite=overwrite)
   
    return out 
        

def read_and_concat_planelogs(planelog_dir: str, spdb_yaml:str, config_yaml:str, 
                              as_xarr:bool=False, convert2_molmol: bool = False, 
                              output_dir:str='', output_file:str='', overwrite:bool=False):
    """
    Function to read in all plane.log files within a directory and concatenate them
    into either a single pandas dataframe or into an xarray dataset indexed in UTC
    time. If you only want to read in a single file, use read_planelog() instead.
    
    INPUT: 
    -------
    (1) planelog_dir   -  STRING containing absolute path to a directory containing 
                          GEOSChem output plane.log files.
                          
    (2) spdb_yaml      -  STRING containing absolute path to the species_database.yml 
                          file associated with your run. Used to add metadata. 
                          
    (3) config_yaml    -  STRING containing path to geoschem_config.yml file 
                          associated with your run. Used to add metadata.
    
    (4) as_xarr        -  (OPTIONAL) BOOLEAN indicating if you want to output data as 
                          an xarray dataset instead of as a pandas dataframe. 
                          If this is set to TRUE, the data arrays within the 
                          output dataset will have attributes containing metadata
                          assigned to them. Otherwise, output in the xarray 
                          dataset is the same as that in the df. Default is set 
                          to FALSE (to return a pandas df instead). 
                     
    (5) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all 
                          outputted concentrations into units of mol/mol. 
                          Only variables with concentration units are converted 
                          & others are left in their native unit. Default is set 
                          to FALSE (not to do conversion).If set to "TRUE" and 
                          "TRA_###" appears in plane.log headers, no conversion 
                          is actually performed to prevent against improper usage. 
                          
    (6) output_dir       - STRING containing absolute path to a directory where 
                          concatenated output file will be saved. If none is 
                          passed, then output saved in planelog_dir. 
                          
    (7) output_file     - STRING containing desired name of ouptput file. If none
                          is passed, the output file will be named: 
                              "planelog_concat_[startYYYYMMDD]_[endYYYYMMDD]"
                              
    (8) overwrite       - BOOL indicating that if there is already a file named 
                          "output_file" in the output directory, "output_dir", whether 
                          that file should be overwritten or if a new unique name 
                          should be used. 
                            
    NOTE: Outputted advected species concentrations are in molec/cm3 if Tracer 
          Names are used, but are in mol/mol if Tracer Numbers are used in input files. 
          For more details on why this happens, see GitHub issue #796: 
              https://github.com/geoschem/geos-chem/issues/796  
              
              
    OUTPUT: 
    -------
    If output_xarr==FALSE (Default): 
        
        (1) df        - PANDAS DATAFRAME indexed along a pd.datetime object index named 'time_UTC'. 
                         Saved at: output_dir+output_file+'.csv' 
                 
        (2) info_dict - DICTIONARY containing info about output vars & units.          
                        Saved as yaml at: output_dir+output_file+'.yaml'
                 
    If output_xarr == TRUE: 
        
        (1) ds - XARRAY DATASET indexed along a np.datetime64 type coordinate 
                 & dimension named 'time_UTC'. All data vars within the ds have 
                 individual attributes for "long_name" and "units" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. If this is done an additional 
                 attribute 'note' is included indicating this occurred.
                 Saved at: output_dir+output_file+'.nc'
            
    Regardless of the output type selected, the resulting df/ds contains columns/data arrays 
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are 
    dropped after constructing the time index/coordinate. Column/data_var names are 
    nearly identical to the diagnostics in plane.log, with the exception that any 
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I'). 
    """
    # =========================================================================
    # Check User arguments & raise errors if files/dirs not found  
    # =========================================================================
    # Normalize the directory path (remove trailing slash if present)
    planelog_dir = planelog_dir.rstrip('/')
    
    # Verify that input dirpath exists, otherwise thrown an error: 
    if not os.path.isdir(planelog_dir):
        raise FileNotFoundError('The following input directory passed to'+ 
                                'read_and_concat_planelogs() could not '+ 
                                'be found:\n\t'+planelog_dir).with_traceback(sys.exc_info()[2])
    else: 
        # Verify that at least 1 file matching plane.log.* exists 
        # within the dirpath to read in, otherwise throw an error: 
        files = os.listdir(planelog_dir)
        planelog_files = [f for f in files if f.startswith('plane.log.')]
        if len(planelog_files) == 0:
            raise FileNotFoundError('No files matching "plane.log.*" were found'+
                                    'under the following directory to read in:\n\t'+ 
                                    planelog_dir).with_traceback(sys.exc_info()[2])
            
    # =========================================================================
    # Loop over all plane.log files in directory & concat data into single df
    # =========================================================================  
    for i,file in enumerate(planelog_files):  
        
        # Read raw file in as dataframe (no unit conversions / renaming)  
        df_i = _read_planelog_to_df(os.path.join(planelog_dir,file))
        
        # Add a new column to df indicating the file this data came from: 
        df_i['PlaneLog_File']= np.full(len(df_i),os.path.join(planelog_dir,file))  
        
        # Add a new column to df indicating the index of the file this data came from:  
        df_i['File_Index']= np.full(len(df_i), i+1)  
        
        if i==0 : 
            # W/ first file, define concatenated "df" w/ deep copy of df_i.
            df= df_i.copy() 
        else: 
            # For all subsequent loops, concatenate a deep copy of df_i with
            # existing data in df_all from other plane.log files. Deep copy is 
            # precaution to ensure we're not inadvertently referencing a view. 
            df = pd.concat([df, df_i.copy()])
            
    # =========================================================================
    # Build meta data dictionary with info about all vars in output dataframe.
    # ========================================================================= 
    # This Assumes advected species units are mol/mol if tracer numbers used, but 
    # molec/cm3 if tracer numbers used. 
    #
    # **NOTE**: This line MUST remain before renaming df columns with tracer numbers 
    #           to tracer names for output units to be correct in info_dict since 
    #           assigning them properly requires us to know if the column names 
    #           were originally tracer numbers or tracer names. 
    
    info_dict=_build_output_meta_dict(df, config_yaml,spdb_yaml)
    
    # Get dictionary to map tracer numbers to tracer names: 
    nums2names=_get_tracer_name_num_mapping(config_yaml, keys_are_nums=True)
    
    # Rename all columns in dataframe so columns with tracer numbers are converted 
    # to tracer names. Won't affect anything if cols are already named with tracer names. 
    df.rename(columns=nums2names, inplace=True)
     
    # =========================================================================
    # (OPTIONALLY) Convert units to mol/mol from molec/cm3 if tracer names were used. 
    # ========================================================================= 
    # Convert concentrations to mol/mol from molec/cm3 & update units in info_dict
    # accordingly. Also add note if conversion performed or not ... 
    if convert2_molmol is True: 
        # NOTE: Function will not actually do anything if conc_units!='molec/cm3' 
        # which is the case if tracer numbers were used in analogous input files. 
        df, info_dict = _convert_moleccm3_to_mr(df, info_dict)
    
    # =========================================================================
    # Configure output filename/directory: 
    # =========================================================================
    # Verify that output dirpath exists, otherwise assign out_dir to planelog_dir. 
    if len(output_dir)==0 or not os.path.isdir(output_dir): 
        output_dir=planelog_dir
        if not os.path.isdir(output_dir):
            raise Warning('The output directory passed to read_and_concat_planelogs() could not '+ 
                          'be found:\n\t'+planelog_dir+'/n Output will be saved at '+
                          'input planelog_dir instead.')

        
    if len(output_file)==0: 
        # Create filename from start/end date of planelogs if none passed by extracting
        # the dates from parsed filenames & find the min/max dates 
        dates = []; good_dates=True # Make empty list to save dates from file
        for file in planelog_files:
            date_str = file.split('.')[-1] # Get the date part (expecting 'plane.log.YYYYMMDD')
            if len(date_str) == 8:  # Ensure it's in YYYYMMDD format
                dates.append(date_str)
            else: 
                good_dates=False # Update bool to ID that not all dates could be parsed. 
                raise Warning('Dates extracted from plane.log.YYYYMMDD filenames'+
                              'were not 8 chars long as expected, which can result'+ 
                              'if you changed the format of the planeflight output'+ 
                              'filenames from their default format in your'+
                              'geoschem_config.yml file. So,the output filename of'+
                              'the concatenated data will not contain the date range'+
                              'of all data concatenated within it.').with_traceback(sys.exc_info()[2])
        if good_dates is True: 
            # If we were able to get len(8) dates from all file names, take min/max. 
            earliest= min(dates); latest = max(dates)
            # And use them to define the name of the output file: 
            output_file = 'planelog_concat_'+earliest+'_'+latest
        else: 
            # If we couldn't get dates from all file names, make dateless outputfile: 
             output_file = 'planelog_concat'
    
    # Combine output directory and output filename for abs path to output files
    savefile= os.path.join(output_dir,output_file)
    
    # Save output / convert to xarray if asked & print to screen where outpput saved: 
    out = _save_outputs(df, info_dict,savefile, as_xarr=as_xarr, overwrite=overwrite)
    
    return out 
    



