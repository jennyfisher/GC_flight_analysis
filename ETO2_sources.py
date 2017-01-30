#!/apps/python3/3.5.2/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 17:47:02 2016

@author: jennyf
"""
## For debugging as needed
#import pdb

### LET'S DO THIS!
import os 
import sys
import argparse

from read_airborne_data import *
from read_gc_data import *
from etno3_profiles import *

# Set up arguments to the script
parser = argparse.ArgumentParser()

# positional arguments

# optional arguments
parser.add_argument("--months", type = int, nargs='+', default = [7,9], help="start/end of month(s) to plot -- defaults to 7-9")
parser.add_argument("--datadir", dest = "datadir", default = "/short/m19/jaf574/data/", help="location of data directories on system; defaults to NCI setup")
parser.add_argument("--FRAPPE", help="Overplot FRAPPE data?", action="store_true")
parser.add_argument("--SEAC4RS", help="Overplot SEAC4RS data?", action="store_true")
parser.add_argument("--HIPPO", help="Overplot HIPPO data?", action="store_true")
parser.add_argument("--summary", help="Use summary version of legend (small contributions not shown?", action="store_true")

# store in args
args = parser.parse_args()

def read_files_hippo():

    # Read HIPPO data - only do this once!
    print ('Reading HIPPO observations')
    hippo_obs = read_hippo(args.datadir+'HIPPO/HIPPO_discrete_continuous_merge_20121129.tbl')

    return hippo_obs
    
    
def read_files_seac4rs():
    
    first=True
    DataDir = args.datadir+'SEAC4RS/'

    # Loop over files to read in data
    for file in os.listdir(DataDir):
        
       print ("Reading "+DataDir+file)
   
       if first:
           seac4rs_obs = read_ict(DataDir+file)                 # read data
           first=False                                          # reset logical
       else:
           tmp_data = read_ict(DataDir+file)                    # read data
           seac4rs_obs=numpy.ma.concatenate((seac4rs_obs,tmp_data))    # concatenate, maintaining missing values
    
    return seac4rs_obs

def read_files_frappe():
    
    first=True
    DataDir = args.datadir+'FRAPPE/mrg/'

    # Loop over files to read in data
    for file in os.listdir(DataDir):
        
       print ("Reading "+DataDir+file)
   
       if first:
           frappe_obs = read_ict(DataDir+file)                 # read data
           first=False                                          # reset logical
       else:
           tmp_data = read_ict(DataDir+file)                    # read data
           frappe_obs=numpy.ma.concatenate((frappe_obs,tmp_data))    # concatenate, maintaining missing values
    
    return frappe_obs

#### PROFILES OF DIFFERENT SPECIES FOR DIFFERENT SPECIES & DIFFERENT MONTHS
def profiles(month,hippo_obs=None,seac4rs_obs=None,frappe_obs=None,altrange=[0,12],region='NPac',dalt=0.5): 

    # For ETNO3
    v = 'ETNO3'
    maxdata=30
    
    # Get HIPPO data
    if args.HIPPO:
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
                          
        GC = etno3_profiles(month,SummaryLegend = args.summary,
                        savefig=False,title='HIPPO '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        airdata=hdata,airname='HIPPO_',maxdata=maxdata,dalt=dalt)
                          
    # Get SEAC4RS data
    if args.SEAC4RS:
        lonrange=[-130,-75]
        latrange=[19,55]
        sdata = extract_seac4rs(seac4rs_obs,gcname_to_seacname(v),
                              lonrange=lonrange,latrange=latrange,
                              altrange=altrange,dayrange=dayrange)
        
        GC = etno3_profiles(month,SummaryLegend = args.summary,
                        savefig=False,title='SEAC4RS '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        airdata=sdata,airname='SEAC4RS_',maxdata=maxdata,dalt=dalt)
    
    # Get FRAPPE data
    if args.FRAPPE:
        lonrange=[-110,-100]
        latrange=[37,42]
        fdata = extract_seac4rs(frappe_obs,gcname_to_seacname(v),
                              lonrange=lonrange,latrange=latrange,
                              altrange=altrange,dayrange=dayrange)
        
        GC = etno3_profiles(month,SummaryLegend = args.summary,
                        savefig=False,title='FRAPPE '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        airdata=fdata,airname='FRAPPE_',maxdata=maxdata,dalt=dalt)

    # Model only at this point??
    if GConly:
        lonrange=[-180,180]
        latrange=[-90,90]
        GC = etno3_profiles(month,SummaryLegend = args.summary,
                        savefig=False,title='Global Model '+v,
                        lonrange=lonrange,latrange=latrange,altrange=altrange,
                        maxdata=maxdata)

            
    return GC


# Set-up: Dates
doy0=[1,32,60,91,121,152,182,213,244,274,305,335]
doy1=[31,59,90,120,151,181,212,243,273,304,334,365]
if len(args.months) == 1:
    month = args.months 
else:
    month = numpy.arange(args.months[0],args.months[-1]+1)
dayrange=[doy0[month[0]-1],doy1[month[-1]-1]]

# Set up empty arrays as default
hippo_obs = None
seac4rs_obs = None
frappe_obs = None

# Set-up: Files
if args.FRAPPE:
   frappe_obs = read_files_frappe()
if args.SEAC4RS:
   seac4rs_obs = read_files_seac4rs()
if args.HIPPO:
   hippo_obs = read_files_hippo()

# Model only?
if args.FRAPPE or args.SEAC4RS or args.HIPPO:
    GConly = False
else:
    GConly = True



#### PROFILES OF DIFFERENT SPECIES FOR DIFFERENT SPECIES & DIFFERENT MONTHS
region='NPac'
profiles(month,hippo_obs=hippo_obs,seac4rs_obs=seac4rs_obs,frappe_obs=frappe_obs,region=region)
