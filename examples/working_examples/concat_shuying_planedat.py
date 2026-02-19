#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 11:10:03 2024

@author: u6044586
"""
import os 
import sys
import numpy as np 
import xarray as xr 
import matplotlib.pyplot as plt 
import fnmatch
from scipy.odr import ODR, Model, Data, RealData

# Add Jessica's library to make/read in planeflight input files: 
sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/planeflight_io/')
import planeflight_io as pln

###############################################################################
#   Functions to concatenate planelog output data: 
###############################################################################
def get_paths(camp, nested): 
    # Set paths to output rundir files & real campaign data:  
    if camp=='UWFPS': 
        if nested is False: 
            # UWFPS 2017 Global 2x2.5: 
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base'
            mod_pth= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_rs_ssa'
        else: 
            # UWFPS 2017 Nested Grid
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/UWFPS_base'
            mod_pth= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/UWFPS_rs'
            
        dat_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/UWFPS_2017_1min_Merged.nc'
        dvars={'ClNO2':'ClNO2_pptv', 'HCl':'HCl_pptv', 'Cl2':'Cl2_pptv','HOCl':'HOCl_pptv','BrCl':'BrCl_pptv',
               'N2O5':'N2O5_pptv'}
    elif camp=='SNACK': 
        if nested is False: 
            # Kalamazoo 2018 Global 2x2.5: 
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_Kalamazoo_base'
            mod_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_full_merra_47_2x25_2018_Kalamazoo'
        else: 
            # Kalamazoo 2018 Nested Grid
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/Kalamazoo_base'
            mod_pth= '/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/Kalamazoo_rs'
        # Path to the actual campaign data: 
        dat_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/Kalamazoo_2018_30min_Merged.nc'
        dvars={'ClNO2':'ClNO2_pptv', 'N2O5':'N2O5_pptv'}
        
    elif camp=='WINTER': 
        if nested is False: 
            # WINTER 2015 Global 2x2.5: 
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_WINTER_base'
            mod_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/gc_fullchem_47layer_merra2_WINTER/'
        else: 
            # WINTER 2015 Nested Grid: 
            base_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/winte_base/'   
            mod_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/GEOS_CHEM/GC_RunDirs/nested_grid/winter_rs'
            
        dat_pth='/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao/WINTER_2015_1min_Merged.nc'
        dvars={'ClNO2':'CIMS1_UW_ClNO2', 'HCl':'CIMS1_UW_ClH', 'Cl2':'CIMS1_UW_Cl2','HOCl':'CIMS1_UW_ClHO','N2O5':'CIMS1_UW_N2O5'}
        
    return base_pth, mod_pth, dat_pth, dvars


def find_concat_files(rundir_pths):
    # Turn input into list if str passed: 
    if type(rundir_pths)==str: 
        rundir_pths=[rundir_pths] 
        
    # Create output list to hold results
    files=[] 
    
    # Loop over rundirpaths and fined the planelog concatenated file: 
    for pth in rundir_pths:
        # Set path to planlog outputs in this rundir_path: : 
        planelog_dir=os.path.join(pth,'Plane_Outputs')
        
        # Find the concatenated planelog output filename: 
        ofile =[os.path.join(planelog_dir, f) for f in os.listdir(planelog_dir) if '.nc' in f and 'planelog' not in f]
        
        files=files+ofile
        
    # Return only 1 file as str if only 1 rundir_pth passed rather than list: 
    if len(files)==1: files=files[0]
        
    return files 

def merge_n_check(camp:str, do_base:bool, remerge:bool, nested:bool, comp:bool=False, to_comp:str='clno2'): 
    # Get paths to base run, modified run, real data and name of ClNO2 in data: 
    base_pth, mod_pth, dat_pth, dvars = get_paths(camp, nested)
    
    # Decide if doing concat on base model output or with road salt model output: 
    if do_base==True: 
        pth=base_pth
    else: 
        pth=mod_pth
            
    # Set paths to the output planelog dir, the config and species database files in each rundir: 
    planelog_dir=pth+'/Plane_Outputs/'
    config_yaml=pth+'/geoschem_config.yml'
    spdb_yaml=pth+'/species_database.yml'
    
    if remerge==True: 
        # Concatenate all planelog output files into single xarray ds with metadata attached:  
        model=pln.read_and_concat_planelogs(planelog_dir, spdb_yaml, config_yaml, as_xarr=True, 
                                         convert2_molmol=True, overwrite=True, 
                                         output_dir='/uufs/chpc.utah.edu/common/home/haskins-group1/users/jhask/')
    else: 
        # Find the concatenated planelog output filename: 
        ofile =find_concat_files(pth)
     
        # Load concatenated planedat output: 
        model=xr.open_dataset(ofile)
    
    if comp == True: 
        # Load real data & compare:  
        data=xr.open_dataset(dat_pth)
            
        # Compare model/output data: 
        mx=np.nanmax([np.nanmax(model[to_comp]*1e12), np.nanmax(data[dvars[to_comp]])])
        plt.scatter(data[dvars[to_comp]],model[to_comp]*1e12,color='r', s=5)
        # Plot 1:1 line: 
        r=np.arange(0,mx)
        plt.plot(r,r,color='k') 
        # Set axis lims & labels: 
        plt.xlim([0,mx]); plt.ylim([0,mx])
        plt.xlabel(r'Obs. ['+to_comp+'] (pptv)')
        plt.ylabel(r'Model ['+to_comp+'] (pptv)')
        plt.title(to_comp)
        plt.legend() 
        plt.grid()
        plt.show()
        
        return data, model
    else: 
        return model
###############################################################################
#   Functions to compare base/modified runs & check against data fast: 
###############################################################################
def get_stats(x,y,mx): 
    # Define the linear model
    def linear_model(B, x):
        return B[0] * x + B[1]
    
    # Set up the model
    model = Model(linear_model)
    
    # Prepare the data
    data = RealData(x, y)
    
    # Set up ODR with data and model
    odr = ODR(data, model, beta0=[1., 0.])
    
    # Run the regression
    output = odr.run()
    
    # Extract the parameters
    slope, intercept = output.beta
    
    # Calculate the predicted y values
    y_pred = linear_model([slope, intercept], x)
    
    # Calculate R^2 value
    ss_total = np.sum((y - np.mean(y))**2)
    ss_residual = np.sum((y - y_pred)**2)
    r_sq = 1 - (ss_residual / ss_total)
    
    #Get y to plot: 
    x_line =np.arange(0,mx)
    y_line = linear_model([slope, intercept], x_line)
    
    # Stats: 
    stats=f'm={slope:0.2f}; b={intercept:0.2f}; r$^2$={r_sq:0.2f}'
    
    return x_line, y_line, stats

def sum_cly(to_sum, ds):
    sum_arr = xr.zeros_like(ds[to_sum[0]])
    dims=ds[to_sum[0]].dims[0]
    coords=ds[to_sum[0]].coords
    # Sum all the specified variables
    for var in to_sum:
        sum_arr += ds[var]
    ds['Cly']=xr.DataArray(data=sum_arr, dims=dims, coords=coords)
    
    return ds

def compare_base_modified(camp,list_comp): 
    
    # Get paths to base run, modified run, real data and name of ClNO2 in data: 
    base_pth, mod_pth, dat_pth, dvars = get_paths(camp)
    
    # Get full paths to base/modified model planelog concatenated files: 
    ofiles =find_concat_files([base_pth,mod_pth])
    
    # Load in base/modified planelog concatenated data & Real data 
    base=xr.open_dataset(ofiles[0]) # without road salt
    mod=xr.open_dataset(ofiles[1])  # with road salt
    data=xr.open_dataset(dat_pth)   # data 
    
    # Create a figure and a grid of subplots
    fig, ax = plt.subplots(nrows=3, ncols=2,figsize=(8, 8))
    
    # Calc total Cly in base/model/data:
    cly=[c for c in list_comp if c not in ['N2O5','HNO3']] 
    print(cly)
    base=sum_cly(cly, base)
    mod=sum_cly(cly, mod)
    data=sum_cly([dvars[v] for v in cly], data)
    
    # Add Cly to list to compare & data vars dict: 
    list_comp.append('Cly')
    dvars['Cly']='Cly'
    
    ct=0 # initialize counter to index list_comp
    for row in range(0,3):
        for col in range(0,2): 
            ax_i=ax[row,col]
            to_comp=list_comp[ct]
            ct=ct+1
    
            # Get var name of thing to compare: 
            mvar=to_comp
            if to_comp in list(dvars.keys()):
                dvar=dvars[to_comp]
            else: 
                print('Var {to_comp} not in data variables. Includes: {data.data_vars}')
                sys.exit()
                
            # Get max of model & data to set axis lims: 
            mx=np.nanmax([np.nanmax(mod[mvar]*1e12), np.nanmax(base[mvar]*1e12), np.nanmax(data[dvar])])
            
            # Get stats of comparison: 
            for i,run in enumerate([base, mod]):
                if i==0: 
                    lbl=r'W/o RS: '; color='r'
                else:
                    lbl=r'W/ RS: '; color='b'
                if camp=='WINTER': 
                    inds=np.where(data['RAF_PSXC']>900)
                    data_var=data[dvar][inds[0]]
                    run_var=run[mvar][inds[0]]
                else: 
                    data_var=data[dvar]
                    run_var=run[mvar]
                      
                ok=np.where(np.isnan(data_var)==0)
                x= data_var[ok]
                y=run_var[ok]*1e12
            
                    
                x_line, y_line, stats= get_stats(x,y,mx)
                
                # Compare model/output data with/without road salt:  
                ax_i.plot(x_line,y_line,color=color)
                ax_i.scatter(x,y,color=color, s=5, label=lbl+stats)
        
        
            # Plot 1:1 line: 
            r=np.arange(0,mx)
            ax_i.plot(r,r,color='k') 
            
            # Set axis lims & labels: 
            ax_i.set_xlim([0,mx]); ax_i.set_ylim([0,mx])
            ax_i.set_title(to_comp)
            ax_i.set_xlabel('Obs. ['+mvar+'] (pptv)')
            ax_i.set_ylabel('Model ['+mvar+'] (pptv)')
            ax_i.legend() 
            ax_i.grid()
    plt.tight_layout()  
    plt.show()
    
    return base, mod, data


###############################################################################
# Decide which campaigns / runs to concatenate / compare results
###############################################################################
camp='UWFPS'
do_base=True
remerge=True
nested=True

# #Reconcat the output plane files: 
model=merge_n_check(camp, do_base, remerge, nested, comp=False) 
    
# #Compare the base and modified results: 
#compare_base_modified(camp,['ClNO2','HCl','Cl2','N2O5','HOCl'])









