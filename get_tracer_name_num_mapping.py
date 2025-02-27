# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 16:26:47 2022

@author: jhask
"""
def get_tracer_name_num_mapping(gckpp_Parameters):
    """Function to read the gckpp_Parameters.F90 file & create dict of tracer #s & index #s. """ 
    file = open(gckpp_Parameters, 'r')
    lines = file.readlines();
    
    idd=dict({}); idd_rev=dict({}); 
    for l,ln in enumerate(lines): 
         if 'INTEGER, PARAMETER ::' in ln: 
             ln=ln.replace('INTEGER, PARAMETER ::','')
             lnls=ln.split('=')
             species=lnls[0].replace('ind_','').replace(' ','').replace('\n','')
             index=int(lnls[1].replace(' ' ,'').replace('\n',''))
             if species not in list(idd.keys()): idd[species]=index 
             if index not in list(idd_rev.keys()):idd_rev[index]=species 
    
    file.close() 
    
    return idd, idd_rev
      
idd, idd_rev = get_tracer_name_num_mapping('C:/Users/jhask/OneDrive/Desktop/gckpp_Parameters.F90')


these=['CO','O3','H2O','H2O2',
'SO2','DMS','MSA','SO4','SO4s','SO4H1','SO4H2','SO4D1','SO4D2','SO4D3','SO4D4',
'HNO3','NO','NO2','N2O5', 'HNO2', 'NIT','NITs','NITD1','NITD1','NITD3','NITD4',
'NH3','NH4',
'HCl','HOCl','SALACL','SALCCL','CH3Cl','HBr','HOBr','BrSALA','BrSALC','BrCl','CH3Br','AERI',
'DST1','DST2','DST3','DST4','BCPI','BCPO','pFe',
'CH2O', 'IEPOXA','IEPOXB','IEPOXD','MBO','GLYX','MGLY',
'SOAGX','SOAIE','SOAME','SOSMG','LVOCOA','ASOAG1','ASOG2','ASOG3','ASOA1','ASOA2','ASOA3','ASOAN','SOAP','SOAS']

n=list(); nope=list()
for item in these: 
    if item in list(idd.keys()): 
        n.append(idd[item])
    else:
        nope.append(item)

n=sorted(n)
n=[str(item) for item in n]
print(' '.join(n))
