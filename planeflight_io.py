"""
Script to make readable Planelight.dat.YYYYMMDD input files for GEOS-Chem
and read in outputted planeflight.log files in as pandas dataframes with 
option to concatenate files. 

Created on Sun Mar 21 14:21:14 2021

@author: Dr. Jessica D. Haskins
GitHub: @jhaskinsPhD
Email: jessica.haskins@utah.edu 
"""
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

def _display_diags(diags:dict, collections:list): 
    """Internal helper function to parse "collections" within dictionary "diags" 
    and create a single string that can be printed to screen displaying info 
    about all the available diagnostic collections. Called within "get_diag_info().
    
    CHANGE LOG: 
    ------------
    10/01/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
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
        lines.append(f"To retrieve all diagnostics in collection:\n\tpln.get_diag_info(simtype,these_collections=['{coll}'])\n")
    
    # Join all lines together with line breaks & output single string:     
    info= '\n'.join(lines)
    
    return info  


def _get_tracer_name_num_mapping(logfile_pth, n_species, keys_are_nums:bool=True): 
    """Internal helper function to create dictionary associating tracer name with 
    tracer index numbers. Created by parsing the logfile & using info the planeflight
    diagnostic prints out to determine mapping is right even across all dif 
    versions of GEOS-Chem which may have new species w/ dif indexes. 
    
    INPUT: 
    ------
    
        (1) logfile_pth  - STRING containing absolute path to a GEOS_Chem log file 
                           for the simulation where planeflight diagnostic is 
                           TURNED ON. 
                           
        (2) n_species    - INTEGER containing the total number of advected species 
                           in the simulation. Get this from len(adv_species) from 
                           output of: 
                               adv_species, simtype= pln._parse_geoschem_config(config_yaml)
    
        (2) key_are_nums - BOOLEAN indicating if you want the output dictionary to 
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
                 
    CHANGE LOG: 
    ------------
    10/02/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
    
    """
    # Verify the logfile exists at the input path: 
    if not os.path.isfile(logfile_pth):
        raise FileNotFoundError('The input GEOS-Chem logfile could not be found at the following path:\n\t'+
                                logfile_pth).with_traceback(sys.exc_info()[2])   
        
    # Open log file and read in all lines to a list. 
    fin = open(logfile_pth, 'r')  
    lines = fin.readlines() 
    
    # Figure out where info on the tracer name/ number from Planeflight is within the log file 
    # (for each planeflight.dat) file it opens- in case dif species are tracked in dif files:  
    starts=[index for index, l in enumerate(lines) if  'P L A N E   F L I G H T   D I A G N O S T I C'==l.replace('\n','')]
    stops= [index for index, l in enumerate(lines) if  'Number of flight track points' in l]
    
    # Check that we actually found some planeflight diagnostic output in the log file to parse:  
    if len(starts) == 0: 
        raise ValueError("No planeflight diagnostic output could be found in the logfile:\n\t"+
                         logfile_pth+'\nPlease check that the string, "'+ 
                         '"P L A N E   F L I G H T   D I A G N O S T I C" ',
                         "does indeed appear in your log file.\n If your log " + 
                         "file DOES NOT contain this string, you are using this " + 
                         "function incorrectly (e.g. should only be used on log "+ 
                         "files from runs in which the planeflight diagnostic " + 
                         "was turned on).\nIf your log file DOES contain this " + 
                         "string, then this error message indicates a bug in the"+ 
                         "function _get_tracer_name_num_mapping(), so please raise" +
                         "an issue on the GCPy GitHub page.")
    
    # Check that we got an equal number start/stop indexes to parse:  
    if len(starts) !=len(stops):
        raise ValueError("Start/Stop indexes of line #s in log file to parse for"+ 
                         "planeflight diagnostic output were not equal. It's"+  
                         "possible the planeflight diagnostic output written to"+  
                         "the GEOS-Chem log file may have changed its structure,"+  
                         "thereby breaking this function's ability to determine which"+ 
                         "lines of the log file contain planeflight diagnostic"+  
                         "output. If you get this error, its likely a bug, so "+ 
                         "please raise an issue on the GCPy GitHub page.")
        
    # Create a dictionary to hold the tracer names/nums for each planeflight file
    d=dict({}) 
    
    # Loop over the pairs of start/stop indexes to parse vars traced w/ all 
    # planedat files: 
    for i in range(0,len(starts)-1): 
        
        # Loop over all lines in log file between this iteration of start/stop: 
        for l in lines[starts[i]:stops[i]]:
            # Split the line on the spaces: 
            ln_ls = l.split(' ')
            
            # Ignore any lines that are all '-----' and any empty item (e.g. what happens w/ multiple spaces): 
            ln_ls = [item.replace('\n','') for item in ln_ls if item != '' and not all(char == '-' for char in item.replace('\n',''))]
            
            # If the line has 3 items now (#, Species, PVAR), then we want the 
            # species name & PVAR associated with it to build a dict mapping the
            # tracer # referenced (TRA_PVAR#) in the planeflight output file 
            # with the name of the tracer name it corresponds to in this simulation.
            # We need to ignore the header line (which also only has 3 items in it...) 
            if ((len(ln_ls) == 3) and (ln_ls!=['#','Species','PVAR'])): 
                # Extract the tracer name ONLY & ensure it has no spaces: 
                tracer_name=ln_ls[1].replace(' ', '')
                
                # Only add this tracer name/ number combo to dict if its NOT already there
                if tracer_name not in list(d.keys()):
        
                    # Extract the PVAR index associated with this species: 
                    tracer_num=int(ln_ls[2].replace('#','').replace(' ',''))
                    
                    # ONLY species concentrations are referenced with TRA_###, 
                    # not additional diagnostics like RO2, or GMAO_TEMP. 
                    # So, only add tracer #s for species concentration vars
                    #(e.g. species conc are only those with PVAR <= # advected species)
                    if tracer_num <= n_species: 
                        d[tracer_name]=f'TRA_{tracer_num:03}'
                    else: 
                        d[tracer_name]=tracer_name
                
    if keys_are_nums is True: 
        # If set to true, reverse dictionary before outputting it so tracer 
        # numbers (e.g. "TRA_052") are keys and values are tracer names ('ClNO2').  
        name_map={value: key for key, value in d.items()}
    else:     
        # If set to False, output dict as is with keys being tracer names ('ClNO2') 
        # and values are tracer numbers (e.g. 'TRA_052'). 
        name_map=d.copy() 
        
    return name_map


def _parse_geoschem_config(config_yaml: str):
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
        
    CHANGE LOG: 
    ------------
    10/01/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
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

          
def _check_str_arr_inputs(item, item_nm:str, allowed_vals, allowed_nm:str,
                          wildcard_str:bool=False ): 
    """ Internal helper function to check that input for "item" is EITHER a 
    string wildcard  matching '?ALL?' (if wildcard_str) is set to TRUE, OR is
    array like, with all array elements as type(string), and with all elements
    are within "allowed_vars". Returns if an error is found & why the error happened,
.   Function only called within make_planeflight_inputs() to help with user 
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

def _get_std_outputs(): 
    
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
        
    # Paull info about std_output_collection & place in diag_dict: 
    diag_dict=dict({})
    for key in diags['std_output']['Diagnostics']: 
        diag_dict[key]=diags['std_output']['Diagnostics'][key]
            
    return diag_dict 

        
def _convert2_xr(df:pd.core.frame.DataFrame, conc_units:str, spdb_yaml: str):
    """Internal helper function to parse the species_database.yml file and 
    pulls out info about compounds outputted in plane.log files which it adds 
    as attributes of all xarray datarrays inside the input dataset w/ 
    matching names. This function only used if the `output_xarray` option is 
    set to TRUE in read_planelog() or read_and_concat_planelogs(). 
    
    INPUT: 
    ------ 
     (1) df- Pandas dataframe with plane.log output 
          
     (2) conc_units - STRING containing what units species concentrations are in. 
     
     (3) spdb_yaml - STRING containing absolute path to the species database yaml file 

    
    OUTPUT: 
    ------ 
    (1) ds - Xarray dataset with all columns in df now as data_arrays that each 
             include, at minimum, attributes for "FULLNAME" and "UNITS" using 
             info from the species database.yml and the planeflight_diagnostics.yaml. 
             
    CHANGE LOG: 
    ------------
    10/02/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
    """
    
    # Ensure 'time_UTC' is index of the df (messes up dims/coords of ds if not)  
    if ((df.index.name!='time_UTC') and ('time_UTC' in df.columns)):
        df = df.set_index('time_UTC')
            
    # Convert the concatenated df to an xarray dataset 
    ds = df.to_xarray() 
    
    # Ensure 'time_UTC' is set as a coordinate.
    ds = ds.set_coords('time_UTC')
    
    # Check that the species database file exists, otherwise throw an error: 
    if not os.path.isfile(spdb_yaml):
        raise FileNotFoundError('The following file could not be found:\n\t'+
                                spdb_yaml).with_traceback(sys.exc_info()[2])
    
    # Load the species database as a dict from the species_database.yml file: 
    with open(spdb_yaml, 'r') as f:
        spdb=yaml.load(f, Loader=yaml.FullLoader)
    
    # Extract info about any other additional diagnostics outputted in plane.log: 
    diags= get_diag_info(print_options=False, return_all=True) 

    # Extract info about standard output vars in plane.log: 
    std_out=_get_std_outputs()
    
    # Combine std_out and diags dicts: 
    diags.update(std_out)

    # Loop over all vars outputted in plane.log: 
    for sp in list(ds.data_vars): 
        
        # If we have info about the var in the species database dictionary:: 
        if sp in list(spdb.keys()):
            
            # Dump all info about this var into "tracer_info" dict
            tracer_info=dict({}) 
            for key in list(spdb[sp].keys()):
                
                value=spdb[sp][key]
                # Ensure values aren't bools, but strs (xarr won't save bool attrs!) 
                if isinstance(value, bool): value=str(value)
                tracer_info[key]=value

            # Set units based on input (for species concentrations) 
            tracer_info['UNITS']=conc_units 
            
            # Add all info from tracer info dict as attributes of this var in the output ds: 
            for key in list(tracer_info.keys()): 
                value=tracer_info[key]
                # Ensure values aren't bools, but strs (xarr won't save bool attrs!) 
                if isinstance(value, bool): value=str(value)
                ds[sp].attrs[key] = value
                
            # Delete "tracer_info" so no info about older tracers gets assigned to wrong var.
            del tracer_info 
            
        # Otherwise, if we have info about that var in the optional diagnotistic dictionary: 
        elif sp in list(diags.keys()): 
            # Add all info about the diagnostic as an attribute of this var in the output ds 
            for key in list(diags[sp].keys()): 
                value= diags[sp][key]
                # Ensure values aren't bools, but strs (xarr won't save bool attrs!) 
                if isinstance(value, bool): value=str(value)
                ds[sp].attrs[key] = value

    return ds





def _convert_moleccm3_to_mr(df:pd.core.frame.DataFrame, conc_units:str):
    """Concentrations in plane.log files are outputted in units of molec cm-3 
    if **TRACER NAMES** rather than **TRACER NUMBERS** were used in the analogous 
    input files. This function uses the Pressure/Temperature outputs to convert 
    *ONLY CONCENTRATION* outputsfrom molec cm-3 into mol mol-1. Skipping of vars 
    that are not concentrations is automated herein. 
    
    INPUT: 
    ------ 
      (1) df         -  PANDAS DATAFRAME containing plane.log output.
      
      (2) conc_units -  STRING containing units of all species concentrations within 
                        the input dataframe. 
              
    OUTPUT: 
    ------ 
      (1) df         -  Same as input, but with all species concentration columns
                        now expressed in units of mol/mol instead of molec/cm3.
                        
      (2) conc_units -  Same as input, but with new units that species concentration
                        columns are expressed in, if updated. 
                 
    Change Log: 
    ----------- 
    09/24/24:  Written by Jessica D. Haskins (Email: jessica.haskins@utah.edu | Github: @jdhask)
    """
    
    # Only do conversion if input conc_units are molec/cm3
    if conc_units !='molec/cm3': 
        print('Species concentrations are not in molec/cm3. No conversion performed') 
        return df, conc_units 
    else: 
        #==========================================================================
        # Check to make sure Pressure/ Temperature variables exist in the output 
        # (since they're required to convert to mol/mol from molec/cm3):
        #==========================================================================
        if 'PRESS' in df.columns: 
            pres_var='PRESS'
        elif 'GAMO_PRES' in df.columns:
            pres_var='GAMO_PRES'
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
        # Use n= PV/ RT to get # of mols in V=1 cm3 of air at ambient Pres/Temp:  
        #==========================================================================
        # Get all vars in correct units for n= PV/RT to come out in mols: 
        T= df.GMAO_TEMP.values    # Temp in K 
        R=8314.462                # Gas Constant in L Pa K-1 mol-1
        P=df[pres_var].values*100 # Pres in hPa, convert to Pa: 1 Pa =  100*P(hPa) 
        V=1e-3                    # Vol=1 cm3, convert to L: 1 cm3 = 1e-3 L
        nmols=(P*V)/(R*T)         # Number of mols in 1 cm3 of air at ambient P/T... 
        
        # Define unit_conversion as a single multiplication factor:
            
        #   X molec (CMPD)   {      1 mol         1 cm3 (air)  }   Y mols (CMPD)
        #  --------------- * { -------------- * -------------- } =  -------------
        #    1 cm3 (air)     { 6.022e23 molec    n mols (air)  }     mols (air)
        #
        #                                   ^^
        #                    {       moleccm3_to_molmol        }
        moleccm3_to_molmol = (1/(6.022e23*nmols))
    
        #==========================================================================
        # Make a list of all output that's NOT in molec/cm3 to ignore in conversion
        #==========================================================================
        # This list should consist of Lat/Lon/Pres/Temp/Time/Point/Type/Index #s, 
        # and **ALL** optional diagnostic quantities. Most diagnotics are not concentration
        # units, but those that are (e.g. 'RO2', 'AN', 'NOy',and 'GMAO_AVGW') are 
        # **NOT** output in molec/cm3, but actually in mol/mol, so ALL optional 
        # diagnostics should NOT be converted from molec/cm3 to mol/mol.  
        
        dont_convert = ['POINT', 'TYPE', 'YYYYMMDD', 'HHMM', 'LAT', 'LON', 'PRESS', 
                       'OBS','T-IND', 'P-I', 'I-IND', 'J-IND', 'TIME_LT',
                       'PlaneLog_File']+list(get_diag_info(print_options=False, return_all=True).keys())
        
        dont_convert=dont_convert+[col for col in df.columns for typ in [np.str_, str] if type(df[col].values[0])!= typ]
        
        # Convert all values to units of mol/mol for all outputs not in "dont_convert"
        for item in df.columns:
            if item not in dont_convert: 
                df[item]=df[item]*moleccm3_to_molmol
                
        # Update output conc_units to be 'mol/mol' now... 
        conc_units= 'mol/mol'
        
        return df, conc_units

###############################################################################
# MAIN FUNCTIONS- INTENDED TO BE USED 
###############################################################################        
def get_diag_info(simtype:str='', print_options:bool=False, return_all:bool=False, 
                  these_collections:list=[]): 
    """Function to load all optional planeflight diagnotics in an output dicitonary
    as keys with info on diagnostic's FULLNAME/UNITS stored in values. Can be used  
    that are compatiable with a user's simulation type & optionally print to screen. 
    
    INPUT: 
    ------ 
    
    (1) simtype           - (OPTIONAL-ISH) STRING containing user's simulation type. 
                            Must be one of the following (case insensitive): 
                                ['','fullchem','aerosol','carbon','Hg','POPs','tagO3',
                                 'TransportTracers','metals','CH4','CO2','tagCO']
                            Used to determine what available diagnostics are 
                            compatiable with your simulation type. Default is blank 
                            string, '', in which case, all possible diagnostics 
                            will be returned (e.g. return_all is set to TRUE, 
                            even if passed as FALSE, but warning will print). 
    
    (2) print_options     - (OPTIONAL): BOOL indicating if diagnotics should be 
                            printed to screen. if set to FALSE, only dictionary 
                            returned (default behavior). If set to TRUE, a dictionary 
                            with diagnostics will be returned, but info will be
                            ALSO printed to screen & nicely formatted. 
    
    (3) return_all        - (OPTIONAL): BOOL indicating if ALL diagnostics, not 
                            *just* those compatible with user's input "simtype" 
                            are printed/returned in output dict. 
    
    (4) these_collections - (OPTIONAL): LIST of STRINGS containing diagnostic 
                            collection names that you would like to retrieve. 
                            Useful if you only want to get specific collections 
                            of diags, not just ALL diags compatiable with simtype: 
                            Valid options to pass include: 
                              ['aer_uptake', 'aodb', 'aodc', 'aq_aer', 'broken', 
                               'chem_fams', 'defaults', 'gmao_ice', 'gmao_met', 
                               'hg', 'htep', 'tomas']
                            If print_options is set to TRUE, value to pass for 
                            "these_collections" is listed for each. 
                                
    OUTPUT: 
    ------
    
    (1) diag_dict         - DICTIONARY with keys corresponding to diagnostic names
                            used in planeflight.dat input files / plane.log output 
                            iles with VALUES being a dictionary with "FULLNAME" and 
                            "UNITS" keys describing what the diagnostic is / its units. 
                            This info is used to print info to screen / assign attrs 
                            in output dataset if _convert2_xr() is called within 
                            read_planelog() or read_and_concat_planelogs(). 
                 
    CHANGE LOG: 
    ----------
        10/02/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)  
    """ 
        
    # Note: Info in planeflight_diagnostics.yml compiled from within planeflight_mod.F90 
    #       and from the Planefligth Diagnostic ReadTheDocs page: 
    # https://geos-chem.readthedocs.io/en/stable/gcclassic-user-guide/planeflight.html
    
    
    # Check that 'simtype' is one of the allowed values (e.g. any option assigned for 
    # 'sim_name' in createRunDir.sh or blank (returns all not broken!))
    allowed=['fullchem','aerosol','carbon','Hg','POPs','tagO3','TransportTracers',
             'metals','CH4','CO2','tagCO','']
    if simtype.lower() not in [name.lower() for name in allowed]: 
        raise ValueError("Invalid input for 'simtype'='"+simtype+"'"+ 
                          'Valid options include (case insensitive): \n\t'+ 
                          ','.join(allowed)).with_traceback(sys.exc_info()[2])
        
    # User can only skip passing 'simtype' in the case that "return_all is True". 
    # Check if they have used it incorrectly and set return_all to true in this case. 
    if len(simtype)==0 and return_all==False: 
        print('If no "simtype" is passed to pln.get_diag_info(), then ' + 
              '"return_all" must be set to True. Since you did not pass'+ 
              'an arg for "simtype", "return_all" will be set to TRUE '+ 
              'and all diagnostics will be returned.')
        
        return_all=True
    
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
        
        # Always remove the std_output collection. This is just included in planeflight.yml
        # so we can assign attrs for those vars when _convert_to_xarr() is called. 
        # But, since this function is intended to retrieve INPUT diagnotisic options, 
        # always remove it (since you can't reques this in input files). 
        # We use this from planeflight.yml in other places... 
        diags.pop('std_output') 
        
    # Define dict, initially holding strings of all available diagnostic collections : 
    ok_diags=list(diags.keys()) 

    # Remove invalid collections from "ok_diags" based on user's simulation type
    # but don't remove any options if "return_all" is passed: 
    if ((simtype.lower()!='hg') and  (return_all is False)): ok_diags.remove('hg') 
    if ((simtype.lower()!='tomas') and (return_all is False)): ok_diags.remove('tomas')
    if ((simtype.lower()=='carbon') and (return_all is False)): 
        for item in ['aer_uptake','htep','chem_fams']: ok_diags.remove(item) 
    
    # The "broken" collection includes diags that do not work with any simtype b/c
    # they need source code updates in GEOS-Chem... Only return if "all" is requested. 
    if return_all is False: ok_diags.remove('broken') 
    
    # If user asked for specific output diags by listing collections: 
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
                              '\nDo pln.get_diag_info(simtype,print_options=True) to display '+
                              'only the diagnostics compatiable with your simulation type.'
                              ).with_traceback(sys.exc_info()[2]) 
            
        # Define collections so it only includes compatiable collections: 
        collections=[item for item in these_collections if item in ok_diags]
            
    else: # If specific collections weren't requested, return all compatiable options: 
        collections=ok_diags 
        
    # Parse the now curated collections list and build output diags dict: 
    diag_dict=dict({})
    for c in collections: 
        for key in diags[c]['Diagnostics']: 
            diag_dict[key]=diags[c]['Diagnostics'][key]
    
    if print_options is True: 
        # Decide to parse all diags or just those compatible/ those requested. 
        if return_all is True: 
            parse_this= list(diags.keys())
        else: 
            parse_this= collections
        # Format info about diags in parse_this to print to screen: 
        info=_display_diags(diags, parse_this)
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
           
    CHANGE LOG: 
    ------------
    10/02/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
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
    adv_species, simtype=_parse_geoschem_config(gc_config)
    
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
   
        # Get dictionary of all optional diagnostics compatiable with user's simtype. 
        diag_dict= get_diag_info(simtype, print_options= False, return_all=False)  
        
        # Check that input for "diags" is either our wild card string ('?ALL?') or 
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

def read_planelog(filename: str, output_xarr:bool=False, spdb_yaml:str='', 
                  convert2_molmol:bool=False, simtype:str=''):
    """Function to read a single plane.log ouput file into either a pandas dataframe
    or xarray dataset. Handles data split across multiple lines in a plane.log file 
    if a large number of output diagnostics are requested and auto converts 
    the "YYYYMMDD" and "HHMM" date/time diagnostics into a pandas datetime object
    used as either the pandas dataframe index or the xarray coordinate. 
    
    INPUT: 
    ------
        (1) filename       - STRING containing absolute path to a plane.log file. 
        
        (2) output_xarr    - (OPTIONAL) BOOLEAN indicating if you want to output the data as 
                              an xarray dataset instead of as a pandas dataframe. 
                              If this is set to TRUE, the data arrays within the 
                              output dataset will have "LONG_NAME" and "UNITS" attributes
                              assigned to them. But otherwise, output in the xarray 
                              dataset is the same as that in the df. Default is set 
                              to FALSE (to return a pandas df instead). 
                          
        (3) spdb_yaml       - (OPTIONAL) STRING containing absolute path to the 
                              species_database.yml file associated with your run. 
                              **Only relevant if output_xarr is set to TRUE to assign 
                              "FULLNAME" and other attributes to each data_array 
                              within the output dataset.** 
        
        (4) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all 
                              outputted concentrations into units of mol/mol. 
                              Only variables with concentration units are 
                              converted & others are left in their native unit. 
                              Default is set to FALSE (not to do conversion).
                              If set to "TRUE" and "TRA_###" appears in plane.log headers, 
                              no conversion is actually performed to prevent against 
                              improper usage. 
                            
              
        NOTE: Users should only set convert2_molmol to TRUE if **Tracer Names**, (e.g. ISOP) 
              rather than **Tracer Numbers ** (e.g. TRA_001) were used to list 
              the desired diagnostic outputs in their plane.dat input files. This 
              is because outputted concentrations are in molec/cm3 if Tracer 
              Names are used, but are in mol/mol if Tracer Numbers are used. 
              For more details on why this happens, see GitHub issue #796: 
                  https://github.com/geoschem/geos-chem/issues/796   
                   
    OUTPUT: 
    -------
    If output_xarr==FALSE (Default): 
         (1) df -  A PANDAS DATAFRAME indexed along a pd.datetime object 
                 index named 'time_UTC' constructed from the "YYYYMMDD" and 
                 "HHMM" outputs in plane.log. No info about UNITS or LONG_NAME
                 are included in the output df, but can be viewed /grabbed by calling: 
                     
                 diag_info= get_diag_info(sim_type='fullchem', print_options=True)

    If output_xarr == TRUE: 
        (1) ds - An XARRAY DATASET indexed along a np.datetime64 type coordinate 
                 & dimension named 'time_UTC' constructed from the "YYYYMMDD" and 
                 "HHMM" outputs in plane.log. All data vars within the ds have 
                 individual attributes for "FULLNAME" and "UNITS" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. 
                 
    Regardless of the output type selected, the resulting df/ds contains columns/data arrays 
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are 
    dropped after constructing the time index/coordinate. Column/data_var names are 
    nearly identical to the diagnostics in plane.log, with the exception that any 
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I').
    Outputs remain in their native units, unless convert2_molmol is set to TRUE 
    in which case concentrations are assumed to be in molec/cm3 and are converted 
    to mol/mol at ambient P/T.             
                 
    USAGE: 
    ------ 
         Can be used to read in data from a single plane.log file as follows: 
             
             # Define file paths: 
             filename= '/my/path/to/plane.log.20160101'
             spd_yaml='/my/path/to/species_database.yml'
             
             # For pandas df output: 
             df= read_planelog(filename, convert2_molmol=True)
             
             # For xarray ds output: 
             ds= read_planelog(filename, output_xarr=True, spdb_yaml=spdb_yaml,
                               convert2_molmol=True)
             
         This function is also called within "read_and_concat_planelogs()" to 
         concatenate output from all plane.log files in a directory into a 
         single df/ds as well. 
         
    CHANGE LOG: 
    ------------
    05/20/21- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
    
    09/25/24- JDH removed depency on astropy ascii library & changed method of 
             dealing with data spread across multiple lines to be more efficnet/
             better contained. Also added optional convert2_molmol functionality 
             within this function, instead of only in read_and_concat_planelogs(). 
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
    
    #==========================================================================
    #  (OPTIONALLY) CONVERT CONC UNITS/ CONVERT TO XARRAY DS & DO OUTPUT 
    #==========================================================================
    # Figure out if tracer #s are used (species concs outputted in mol/mol) or if 
    # tracer names are used (species conc outputted in molec/cm3): 
    if any('TRA_' in col for col in df.columns): 
        conc_units='mol/mol'
    else: 
        conc_units='molec/cm3' 
       
    # Convert concentrations to mol/mol from molec/cm3 & update conc_units accordingly. 
    if convert2_molmol is True: 
       
        # Note: Function will not actually do anything if conc_units!='molec/cm3' 
        df, conc_units = _convert_moleccm3_to_mr(df, conc_units)
    
    # Convert df to xarray ds if asked. 
    if output_xarr is True: 
        # Function handles conversion to ds, attaches "FULLNAME" and "UNITS" 
        # attributes to individual data vars plus any other info from sp database
        # and from planeflight_diagnostics.yml 
        ds=_convert2_xr(df, conc_units, spdb_yaml)
        
        return ds
        
    else: 
        # Otherwise return output as pandas dataframe
        return df 

def read_and_concat_planelogs(dirpath: str, output_xarr:bool=False, spdb_yaml:str='',
                   convert2_molmol: bool = False ):
    """
    Function to read in all plane.log files within a directory and concatenate them
    into either a single pandas dataframe or into an xarray dataset indexed in UTC
    time. If you only want to read in a single file, use read_planelog() instead.
    
    INPUT: 
    -------
        (1) dirpath        -  STRING containing filepath to a directory containing 
                              GEOSChem output plane.log files
     
        (2) output_xarr    - (OPTIONAL) BOOLEAN indicating if you want to output the data as 
                              an xarray dataset instead of as a pandas dataframe (default). 
                              If this is set to TRUE, the data arrays within the 
                              output dataset will have "LONG_NAME" and "UNITS" attributes
                              assigned to them. But otherwise, output in the xarray 
                              dataset is the same as that in the df. Default is set 
                              to FALSE (to return a pandas df instead). 
                          
        (3) spdb_yaml       - (OPTIONAL) STRING containing absolute path to the 
                              species_database.yml file associated with your run. 
                              **Only relevant if output_xarr is set to TRUE to assign 
                              "FULLNAME" and other attributes to each data_array 
                              within the output dataset.** 
        
        (4) convert2_molmol - (OPTIONAL) BOOLEAN indicating if you wish to convert all 
                              outputted concentrations into units of mol/mol. 
                              Only variables with concentration units are 
                              converted & others are left in their native unit. 
                              Default is set to FALSE (not to do conversion).
                              If set to "TRUE" and "TRA_###" appears in plane.log headers, 
                              no conversion is actually performed to prevent against 
                              improper usage. 
              
        NOTE: Users should only set convert2_molmol to TRUE if **Tracer Names**, (e.g. ISOP) 
              rather than **Tracer Numbers ** (e.g. TRA_001) were used to list 
              the desired diagnostic outputs in their plane.dat input files. This 
              is because outputted concentrations are in molec/cm3 if Tracer 
              Names are used, but are in mol/mol if Tracer Numbers are used. 
              For more details on why this happens, see GitHub issue #796: 
                  https://github.com/geoschem/geos-chem/issues/796
                  
                  
    OUTPUT: 
    -------
    If output_xarr==FALSE (Default): 
         (1) df -  A PANDAS DATAFRAME indexed along a pd.datetime object 
                 index named 'time_UTC' constructed from the "YYYYMMDD" and 
                 "HHMM" outputs in plane.log. No info about UNITS or FULLNAME
                 are included in the output df, but can be viewed /grabbed by calling: 
                     
                 diag_info= get_diag_info(sim_type='fullchem', print_options=True)

    If output_xarr == TRUE: 
        (1) ds - An XARRAY DATASET indexed along a np.datetime64 type coordinate 
                 & dimension named 'time_UTC' constructed from the "YYYYMMDD" and 
                 "HHMM" outputs in plane.log. All data vars within the ds have 
                 individual attributes for "FULLNAME" and "UNITS" assigned to them
                 and will reflect unit conversion changes to mol/mol from molec/cm3
                 if convert2_molmol is set to TRUE. 
                 
    Regardless of the output type selected, the resulting df/ds contains columns/data arrays 
    corresponding to all output diagnostics except "YYYYMMDD" and "HHMM" which are 
    dropped after constructing the time index/coordinate. Column/data_var names are 
    nearly identical to the diagnostics in plane.log, with the exception that any 
    '-' characters are replaced with '_' (e.g. 'P-I' in plane.log becomes 'P_I').
    Outputs remain in their native units, unless convert2_molmol is set to TRUE 
    and 'TRA_XXX' appars in the output headers, in which case concentrations 
    are assumed to be in molec/cm3 and are converted to mol/mol at ambient P/T.  
    
    USAGE: 
    ------ 
         Can be used to read in & concatenate data from multiple plane.log files
         underneath a directory as follows: 
             
             # Define file paths: 
             dirpath= '/my/path/to/planelog_dir/'
             spd_yaml='/my/path/to/species_database.yml'
             
             # For pandas df output: 
             df= read_and_concat_planelogs(dirpath, convert2_molmol=True)
             
             # For xarray ds output: 
             ds= read_planelog(filename, output_xarr=True, spdb_yaml=spdb_yaml,
                               convert2_molmol=True)
             
    CHANGE LOG: 
    ----------
    10/02/24- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)        
    """
    # =========================================================================
    # Check User arguments & raise errors if files/dirs not found  
    # =========================================================================
    
    # Normalize the directory path (remove trailing slash if present)
    dirpath = dirpath.rstrip('/')
    
    # Verify that input dirpath exists, otherwise thrown an error: 
    if not os.path.isdir(dirpath):
        raise FileNotFoundError('The following directory passed to'+ 
                                'read_and_concat_planelogs() could not '+ 
                                'be found:\n\t'+dirpath).with_traceback(sys.exc_info()[2])
    else: 
        # Verify that at least 1 file matching plane.log.* exists 
        # within the dirpath to read in, otherwise throw an error: 
        files = os.listdir(dirpath)
        planelog_files = [f for f in files if f.startswith('plane.log.')]
        if len(planelog_files) == 0:
            raise FileNotFoundError('No files matching "plane.log.*" were found'+
                                    'under the following directory to read in:\n\t'+ 
                                    dirpath).with_traceback(sys.exc_info()[2])
            
    # =========================================================================
    # Loop over all plane.log files in directory & concat data into single df
    # =========================================================================  
    for i,file in enumerate(planelog_files):  
        # Read in data from each file as pandas df & don't do unit conversion yet. 
        # This is done after concatenating all dfs, if requested. 
        df_i= read_planelog(os.path.join(dirpath,file),output_xarr=False,
                            convert2_molmol=False)
        
        # Add a new column to df indicating the file this data came from: 
        df_i['PlaneLog_File']= np.full(len(df_i),os.path.join(dirpath,file))  
        
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
    # Extract the dates from parsed filenames & find the min/max dates to create 
    # the full path to an output file name where all concatenated data is stored. 
    # =========================================================================
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
        # And use them to define the full path to the output file. 
        outfile = os.path.join(dirpath,'planelog_concat_'+earliest+'_'+latest)
    else: 
        # If we couldn't get dates from all file names, make dateless outputfile: 
         outfile = os.path.join(dirpath,'planelog_concat') 
         
    #==========================================================================
    #  (OPTIONALLY) CONVERT CONC UNITS/ CONVERT TO XARRAY DS &  SAVE/ OUTPUT 
    #==========================================================================
    # Figure out if tracer #s are used (species concs outputted in mol/mol) or if 
    # tracer names are used (species conc outputted in molec/cm3): 
    if any('TRA_' in col for col in df.columns): 
        conc_units='mol/mol'
    else: 
        conc_units='molec/cm3' 
        
    # Convert concentrations to mol/mol from molec/cm3 & update conc_units accordingly. 
    if convert2_molmol is True: 
        # Note: Function will not actually do anything if conc_units!='molec/cm3' 
        #df, conc_units = _convert_moleccm3_to_mr(df, conc_units)
        return df
    
    # Convert df to xarray ds if asked. 
    if output_xarr is True: 
        if spdb_yaml=='': 
            raise ValueError('If output_xarr is TRUE, you must pass an argument for'+ 
                             '"spdb_yaml" (str) containing the path to your ' +
                             "simulations's species_database.yml file so that "+
                             'attributes for each data array in the output data set '+
                             'can be filled with info about each species outputted.'+ 
                             'If this functionality is not desired, then simply set'+
                             '"output_xarr=False" to get output as a pandas dataframe '+
                             'with no attributes attached.')
            
        # Function handles conversion to ds, attaches "FULLNAME" and "UNITS" 
        # attributes to individual data vars plus any other info from sp database
        # and from planeflight_diagnostics.yml. Does check for spdb_yaml existance within.
        ds=_convert2_xr(df, conc_units, spdb_yaml)
        
        # Save output xarray dataset as netcdf file
        ds.to_netcdf(outfile+'.nc') 
        
        # Print  output file path to screen & return ds
        print('Concatenated plane.log data saved as xarray dataset at: \n\t'+outfile+'.nc')
        
        return  ds
    else: 
        # Save output pandas dataframe as a csv file: 
        df.to_csv(outfile+'.csv')
        
        # Print output file path to screen & return df
        print('Concatenated plane.log data saved as pandas dataframe at: \n\t'+outfile+'.csv')
        
        return df 
    



