#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:00:33 2025

@author: u6044586
"""
import os 
import yaml 
import warnings 
import sys 







    

  

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



























