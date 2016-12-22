# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:47:02 2016

@author: jennyf
"""

### LET'S DO THIS!
import os 

from map_gc_rono2 import *
from read_airborne_data import *
from read_gc_data import *


def set_files(month):
    # SEAC4RS data in this month?
    SEAC4RS = False
    
    filename=[]
    for m in month:
### KLUDGE
###        filename.append('RG+unlumpPRNO3/gc.2013'+str(m).zfill(2)+'01.nc')
        filename.append('RG/gc.2013'+str(m).zfill(2)+'01.nc')
        if m == 8 or m == 9:
            SEAC4RS=True
            
    return filename, SEAC4RS


def read_files_hippo():
    
    # Read HIPPO data - only do this once!
    print 'Reading HIPPO observations'
    hippo_obs = read_hippo('data/HIPPO_discrete_continuous_merge_20121129.tbl')
    return hippo_obs
    
    
def read_files_seac4rs(SEAC4RS=False):
    
    # READ SEAC4RS DATA HERE - only do this once!
    if SEAC4RS:        
    
        first=True
        DataDir = 'data/SEAC4RS/'
    
        # Loop over files to read in data
        for file in os.listdir(DataDir):
            
           print "Reading "+DataDir+file
       
           if first:
               seac4rs_obs = read_ict(DataDir+file)                 # read data
               first=False                                          # reset logical
           else:
               tmp_data = read_ict(DataDir+file)                    # read data
               seac4rs_obs=numpy.ma.concatenate((seac4rs_obs,tmp_data))    # concatenate, maintaining missing values
    
    return seac4rs_obs


#### PROFILES OF DIFFERENT SPECIES FOR DIFFERENT SPECIES & DIFFERENT MONTHS
def profiles(vars,hippo_obs,seac4rs_obs,altrange=[0,12],region='NPac'): 

    for v in vars:
        
        if v == 'ETNO3_C2H6' or v == 'PRNO3_C3H8':
            maxdata=0.010
        elif v == 'C1-C3_RONO2':
            maxdata=50
        elif v == 'MENO3' or v == 'ETNO3' or v == 'PRNO3':
            maxdata=20
        else:
            maxdata=None
    
        # Get HIPPO data
        lonrange=[-180,180]
    
        # N Pacific
        if region == 'NPac':
            latrange=[20,50]
        # Eq Pacific
        elif region == 'EqPac':
            latrange=[-20,20]
        # S Pacific
        elif region == 'SPac':    
            latrange=[-50,-10]
        else:
            latrange=[-90,90]
        
        hdata = extract_hippo(hippo_obs,gcname_to_hname(v),
                              lonrange=lonrange,latrange=latrange,
                              altrange=altrange,dayrange=dayrange)
                              
        GC = profile_gc(v,filename,savefig=False,title='HIPPO '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        airdata=hdata,airname='HIPPO_',maxdata=maxdata)
                              
        # Get SEAC4RS data
        lonrange=[-130,-75]
        latrange=[19,55]
        sdata = extract_seac4rs(seac4rs_obs,gcname_to_seacname(v),
                              lonrange=lonrange,latrange=latrange,
                              altrange=altrange,dayrange=dayrange)
        
        GC = profile_gc(v,filename,savefig=False,title='SEAC4RS '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        airdata=sdata,airname='SEAC4RS_',maxdata=maxdata)
        
    return GC

##### MAPS OF DIFFERENT SPECIES @ DIFFERENT ALTITUDES & DIFFERENT MONTHS
def maps(vars, hippo_obs, seac4rs_obs, altrange=[0,12]):

    for v in vars:
    
        # Get HIPPO data
        hdata = extract_hippo(hippo_obs,gcname_to_hname(v),
                              altrange=altrange,dayrange=dayrange)
                              
        # Get SEAC4RS data
        sdata = extract_seac4rs(seac4rs_obs,gcname_to_seacname(v),
                                altrange=altrange,dayrange=dayrange)
                                
        # Combine into single dataset
        airdata = airdata_combine(hdata,sdata)
        
        # plot it up    
        #GC=map_gc(v,filename,airdata=hdata,savefig=False,
        #          airname='HIPPO_',altrange=altrange)
        #GC=map_gc(v,filename,airdata=sdata,savefig=False,
        #          airname='SEAC4RS_',altrange=altrange)
    
        GC=map_gc(v,filename,airdata=airdata,savefig=False,
                  airname='SEAC4RS+HIPPO_',altrange=altrange,
                  mindata=0,maxdata=20)
    
    return



# Set-up: Dates
doy0=[1,32,60,91,121,152,182,213,244,274,305,335]
doy1=[31,59,90,120,151,181,212,243,273,304,334,365]
month=range(1,13)
#month=[8,9]
dayrange=[doy0[month[0]-1],doy1[month[-1]-1]]

# Set-up: Files
filename, SEAC4RS = set_files(month)
try:
    hippo_obs
except NameError:
    hippo_obs = read_files_hippo()    
try:
    seac4rs_obs
except NameError:
    seac4rs_obs = read_files_seac4rs(SEAC4RS)


# Set-up: Parameters

#vars=['MeNO3','EtNO3','PrNO3']
vars=['C1-C3_RONO2']
#vars=['MeNO3','EtNO3','iPrNO3','nPrNO3']
#vars=['ETNO3','C2H6','ETNO3_C2H6']
#vars=['PRNO3','C3H8','PRNO3_C3H8']
#vars=['C2H6','C3H8']


##### MAPS OF DIFFERENT SPECIES @ DIFFERENT ALTITUDES & DIFFERENT MONTHS
altrange=[8,12]
#maps(vars,hippo_obs,seac4rs_obs,altrange=altrange)


#### PROFILES OF DIFFERENT SPECIES FOR DIFFERENT SPECIES & DIFFERENT MONTHS
region='NPac'
profiles(vars,hippo_obs,seac4rs_obs,region=region)

