#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:04:44 2024

@author: u6044586
"""

import sys
import os 
import numpy as np 
import yaml 

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
    05/20/21- Written by Prof. Jessica D. Haskins (Github: @jhaskinsPhD)
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

        (1) name_map  -    DICTIONARY with key/value pairs of tracer names (ClNO2) & 
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


def create_num_inputs(in_dir, out_dir, logfile_pth, config_yaml): 
    # Parse the config file & get # of advectoed species: 
    adv_species, simtype= _parse_geoschem_config(config_yaml)
                                                 
    # Parse log file and get dict of tracer name/number pairs: 
    name_map= _get_tracer_name_num_mapping(logfile_pth, len(adv_species), keys_are_nums=False)

    files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
    
    for f in files: 
        # Open file and read in all lines to a list. 
        fin = open(os.path.join(in_dir,f), 'r')  
        lines = fin.readlines() 
        
        new_lns=[]
        for l in lines: 
            got_it=False
            for sp in list(name_map.keys()): 
                if l.replace(' ', '').replace('\n','')==sp and got_it==False: 
                    new_lns.append(name_map[sp]+'\n')
                    got_it=True
            if got_it==False: 
                new_lns.append(l)
    
        # Create new file with 
        outF = open(os.path.join(out_dir,f), 'w')
        for line in new_lns:
            outF.write(line)
            
        del fin, lines, outF, new_lns
        
    return 

def create_num_outputs(in_dir, out_dir, logfile_pth, config_yaml): 

    def update_headers(header_line,name_map ): 
        # We are also checking for leading and trailing spaces
        headers=header_line.split() 
        out_chars = list(header_line)
        for header in headers:
            if header in name_map:
                key = header
                value = name_map[key]
                key_index = header_line.find(key)
        
                while key_index != -1:
                    # Ensure that the key is matched as a whole word surrounded by spaces
                    is_full_match = (
                        (key_index == 0 or header_line[key_index - 1] == ' ') and
                        (key_index + len(key) == len(header_line) or header_line[key_index + len(key)] == ' ')
                    )
                    
                    if is_full_match:
                        # Replace the found key with the corresponding value in out_chars
                        for i in range(len(key)):
                            out_chars[key_index + i] = ' '  # Clear original key
                            
                        for i in range(len(value)):
                            if i < len(value):
                                out_chars[key_index + i] = value[i]  # Insert the new value
                        
                    # Look for the next occurrence
                    key_index = header_line.find(key, key_index + 1)
        
        # Join the list back into a string
        out_str = ''.join(out_chars)
        
        return out_str
    
    # Parse the config file & get # of advected species: 
    adv_species, simtype= _parse_geoschem_config(config_yaml)
                                
    # Parse log file and get dict of tracer name/number pairs: 
    name_map= _get_tracer_name_num_mapping(logfile_pth, len(adv_species), keys_are_nums=False)

    files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]

    for f in files: 
        with open(os.path.join(in_dir,f), 'r') as infile, open(os.path.join(out_dir,f), 'w') as outfile:
    
            # Open file and read in all lines to a list. 
            header_line1 = infile.readline() 
            header_line2 = infile.readline() 
            data_lines=infile.readlines()
        
            new_header1=update_headers(header_line1,name_map)
            new_header2=update_headers(header_line2,name_map)
        
            # Create new file with 
            outfile.write(new_header1)
            outfile.write(new_header2)
            for line in data_lines: 
                outfile.write(line)
    return 

    
config_yaml='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base/geoschem_config_jan_apr.yml'
spdb_path='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base/species_database.yml'
logfile_pth= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base/geoschem_2017_jan_apr.log'


i_names='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/test_plane_data/inputs_names/'
i_nums='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/test_plane_data/inputs_nums/'
o_names='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/test_plane_data/outputs_names/'
o_nums='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/gcpy_campaigns/test_plane_data/outputs_nums/'



## Create input fils with tracer #s from tracer name inputs: 
#create_num_inputs(i_names, i_nums, logfile_pth, config_yaml)   

## Create output files with tracer #s from tracer name outputs: 
#create_num_outputs(o_names, o_nums, logfile_pth, config_yaml):

with open(spdb_path, 'r') as f:
    spdb=yaml.load(f, Loader=yaml.FullLoader)
























