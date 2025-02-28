#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 17:41:45 2025

@author: u6044586
"""

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
