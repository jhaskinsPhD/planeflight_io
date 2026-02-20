#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 15:32:12 2024

@author: Dr. Jessica D. Haskins
GitHub: @jhaskinsPhD
Email: jessica.haskins@utah.edu
"""
import yaml
import os 

# The following collections have been tested & work in GEOS-Chem: 
defaults= dict({"Collection":"Default Diagnostics", 
           "Notes": "Diagnostics in this collection are included by default in all files generated using this function.",
           "Diagnostics": {
                "TIME_LT":  {"LONG_NAME":"Local time",   "UNITS":"hours"},
                "GMAO_TEMP":{"LONG_NAME":"Temperature", "UNITS":"K"}, 
                "GMAO_RELH":{"LONG_NAME":"Relative humidity",  "UNITS":"%"}, 
                "GMAO_PRES":{"LONG_NAME":"Pressure at center of grid box", "UNITS":"hPa"},
                "GMAO_IIEV":{"LONG_NAME":"GEOS-Chem grid box Longitude index", "UNITS":"unitless "}, 
                "GMAO_JJEV":{"LONG_NAME":"GEOS-Chem grid box Latitude index",  "UNITS":"unitless"}, 
                "GMAO_LLEV":{"LONG_NAME":"GEOS-Chem grid box index Level index","UNITS":"unitless"}, 
           }})

fams = dict({"Collection":"Chemical Family Diagnostics", 
           "Notes": "",
           "Diagnostics": {
               "RO2": {"LONG_NAME":"Concentration of RO2 family", "UNITS":"v/v (dry)"}, 
               "AN": {"LONG_NAME":"Concentration of AN family", "UNITS":"v/v (dry)"}, 
               "NOy": {"LONG_NAME":"Concentration of NOy family", "UNITS":"v/v (dry)"}
           }})

gmao = dict({"Collection":"GMAO Grid Box Ice Fraction Diagnostics", 
           "Notes": "",
           "Diagnostics": {
               "GMAO_TEMP":{"LONG_NAME":"Temperature", "UNITS":"K"}, 
               "GMAO_ABSH":{"LONG_NAME":"Absolute humidity", "UNITS":"unitless"}, 
               "GMAO_SURF":{"LONG_NAME":"Aerosol surface area", "UNITS":"cm2/cm3"}, 
               "GMAO_PSFC":{"LONG_NAME":"Surface pressure",  "UNITS":"hPa"}, 
               "GMAO_UWND":{"LONG_NAME":"Zonal winds", "UNITS":"m/s"}, 
               "GMAO_VWND":{"LONG_NAME":"Meridional winds",  "UNITS":"m/s"}, 
               "GMAO_PSLV":{"LONG_NAME":"Sea level pressure",   "UNITS":"hPa"},
               "GMAO_THTA":{"LONG_NAME":"Potential temperature",  "UNITS":"K"}
           }})

ice= dict({"Collection":"GMAO Grid Box Ice Fraction Diagnostics", 
           "Notes": "",
           "Diagnostics": {
               "GMAO_ICE00" : {"LONG_NAME":"Fraction of each grid box that has 0% to +10% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE10" : {"LONG_NAME":"Fraction of each grid box that has 10% to +20% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE20" : {"LONG_NAME":"Fraction of each grid box that has 20% to +30% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE30" : {"LONG_NAME":"Fraction of each grid box that has 30% to +40% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE40" : {"LONG_NAME":"Fraction of each grid box that has 40% to +50% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE50" : {"LONG_NAME":"Fraction of each grid box that has 50% to +60% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE60" : {"LONG_NAME":"Fraction of each grid box that has 60% to +70% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE70" : {"LONG_NAME":"Fraction of each grid box that has 70% to +80% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE80" : {"LONG_NAME":"Fraction of each grid box that has 80% to +90% of sea ice coverage", "UNITS":"unitless"},
               "GMAO_ICE90" : {"LONG_NAME":"Fraction of each grid box that has 90% to +100% of sea ice coverage", "UNITS":"unitless"}
           }})

aods = dict({"Collection":"Column Aerosol Optical Depth Diagonstics", 
           "Notes": "",
           "Diagnostics": {
               "AODC_SULF":{"LONG_NAME":"Column aerosol optical depth for sulfate",   "UNITS":"unitless"},
               "AODC_BLKC":{"LONG_NAME":"Column aerosol optical depth for black carbon",  "UNITS":"unitless"},
               "AODC_ORGC":{"LONG_NAME":"Column aerosol optical depth for organic carbon",  "UNITS":"unitless"},
               "AODC_SALA":{"LONG_NAME":"Column aerosol optical depth for accumulation mode sea salt",  "UNITS":"unitless"},
               "AODC_SALC":{"LONG_NAME":"Column aerosol optical depth for coarse mode sea salt",   "UNITS":"unitless"},
               "AODC_DUST":{"LONG_NAME":"Column aerosol optical depth for dust",  "UNITS":"unitless"},
               "AODC_TOT": {"LONG_NAME":"Column aerosol optical depth for all aerosol types",  "UNITS":"unitless"},
               "AODB_SULF":{"LONG_NAME":"Column aerosol optical depth for sulfate below aircraft",  "UNITS":"unitless"},
               "AODB_BLKC":{"LONG_NAME":"Column aerosol optical depth for black carbon below aircraft",  "UNITS":"unitless"},
               "AODB_ORGC":{"LONG_NAME":"Column aerosol optical depth for organic carbon below aircraft",  "UNITS":"unitless"},
               "AODB_SALA":{"LONG_NAME":"Column aerosol optical depth for accumulation mode sea salt below aircraft",  "UNITS":"unitless"},
               "AODB_SALC":{"LONG_NAME":"Column aerosol optical depth for coarse mode sea salt below aircraft",  "UNITS":"unitless"},
               "AODB_DUST":{"LONG_NAME":"Column aerosol optical depth for dust below aircraft",  "UNITS":"unitless"},
               "AODB_TOT": {"LONG_NAME":"Column aerosol optical depth for all aerosols below aircraft",  "UNITS":"unitless"}
           }})

isor=dict({"Collection":"ISORROPIA/HTEP Diagnostics", 
           "Notes": "",
           "Diagnostics": {
                "ISOR_HPLUS":{"LONG_NAME":"ISORROPIA [H+]",   "UNITS":"M"},
                "ISOR_PH":   {"LONG_NAME":"ISORROPIA pH (Note: In a non-ideal system, pH can be negative!)",  "UNITS":"unitless"},
                "ISOR_AH2O": {"LONG_NAME":"ISORROPIA [aerosol water]",   "UNITS":"ug/m3 (air)"},
                "ISOR_HSO4": {"LONG_NAME":"ISORROPIA [bifulfate]",  "UNITS":"M"}
           }})

aqaer=dict({"Collection":"Aqueous Aerosol Diagnostics", 
           "Notes": "", 
           "Diagnostics": {
               "AQAER_RAD"   : {"LONG_NAME":"Aqueous aerosol radius" , "UNITS":"cm"},
               "AQAER_SURF"  : {"LONG_NAME":"Aqueous aerosol surface area", "UNITS":"cm2/cm3"}
          }})

gamm=dict({"Collection":"Aerosol Uptake Diagnostics", 
           "Notes": "Untested diagnostics. Will not work w/ v<12.8.0 prior to Bates et al.,2019 Isoprene Chemistry!!", 
           "Diagnostics": {
               "GAMM_EPOX"  : {"LONG_NAME":"Uptake coefficient for EPOX", "UNITS":"unitless"},
               "GAMM_IMAE"  : {"LONG_NAME":"Uptake coefficient for IMAE", "UNITS":"unitless"},
               "GAMM_ISOPN" : {"LONG_NAME":"Uptake coefficient for ISOPN", "UNITS":"unitless"},
               "GAMM_DHDN"  : {"LONG_NAME":"Uptake coefficient for DHDN", "UNITS":"unitless"},
               "GAMM_GLYX"  : {"LONG_NAME":"Uptake coefficient for GLYX", "UNITS":"unitless"},
           }}) 

hg_sim=dict({"Collection":"Mercury Diagnostics", 
           "Notes": "Untested diagnostics. Will ONLY work with Hg/PoPs simulations (not in fullchem)!", 
           "Diagnostics": {
               "HG2_FRACG":{"LONG_NAME":"Fraction of Hg(II) in the gas phase",   "UNITS":"unitless"},
               "HG2_FRACP":{"LONG_NAME":"Fraction of Hg(II) in the particle phase",  "UNITS":"unitless"},
            }})

tomas= dict({"Collection":"TOMAS Diagnostics", 
           "Notes": "Untested diagnostics. If you want these, you'll need to add your own TMS_nnn's manually!", 
           "Diagnostics": 
               {'TMS_nnn': {"LONG_NAME": "Nucleation rates (TOMAS)", "UNITS":"No units listed on readthedocs..."}
            }})
            
broken=dict({"Collection":"Production Rates / Reaction Rates", 
            "Notes": "THESE DIAGNOSTICS DO NOT CURRENTLY WORK. GEOS-CHEM SOURCE CODE NEEDS UPDATING", 
            "Diagnostics": {
               'PROD_xxxx': {"LONG_NAME": 'Production Rate of Reaction # xxxx', "UNITS":"molec/cm3/s"},
               'REA_nnn': {"LONG_NAME": 'Rate of Reaction # xxxx', "UNITS":"molec/cm3/s"}
           }})


# Combine all dictionaries into a single dictionary
all_diags = {
    "defaults": defaults,
    "fams": fams,
    "gmao": gmao,
    "ice": ice,
    "aods": aods,
    "isor": isor,
    "aqaer": aqaer,
    "gamm": gamm,
    "hg": hg_sim,
    "tomas": tomas,
    "broken": broken
}



# Function to save all dictionaries to a single YAML file
def save_dictionaries_to_yaml(filename, data):
    with open(filename, 'w') as f:
        yaml.dump(data, f)

# Save to a YAML file
this_dir= os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
savepath=os.path.join(this_dir, 'planeflight_diagnostics.yml')
save_dictionaries_to_yaml(savepath, all_diags)
