# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 10:07:07 2016

@author: jennyf
"""

import numpy
import netCDF4

from PseudoNetCDF import PNC

def read_gc_nc(varname,filename,ppt=True):
    
    """Read a GEOS-Chem netcdf file created using BPCH2COARDS and extract
       lon, lat, altitude (from boxheight), and data field.
       
       Default is to convert to ppt from ppb but this can be overwritten
       with ppb=False"""

    f = netCDF4.Dataset(filename,'r')

    # Essential variables first
    glon = f.variables['lon']
    glat = f.variables['lat']
    gbxht = f.variables['BXHGHT_S__BXHEIGHT']
    
    # use boxheight variable to calculate altitude
    galt = numpy.zeros(gbxht.shape)
    nlev=galt.shape[1]
    for ll in range(nlev):
        # alt is cumulative sum, convert to m
        galt[:,ll,:,:] = numpy.sum(gbxht[:,0:ll,:,:],axis=1)*1e-3

    # conv is factor to convert from ppb if needed
    # defaults to ppt but can be set otherwise, as in special cases
    if ppt:
        conv = 1e3
        unit = 'ppt'
    else:
        conv = 1e0

    # Get data - deal with "special" fields first
    if (varname.upper() == "C1-C3_RONO2"):
        gdata =  ( numpy.array(f.variables['IJ_AVG_S__MENO3']) +
                   numpy.array(f.variables['IJ_AVG_S__ETNO3']) +
                   numpy.array(f.variables['IJ_AVG_S__IPRNO3']) +
                   numpy.array(f.variables['IJ_AVG_S__NPRNO3']) )
    elif (varname.upper() == "PRNO3"):
        gdata =  ( numpy.array(f.variables['IJ_AVG_S__IPRNO3']) +
                   numpy.array(f.variables['IJ_AVG_S__NPRNO3']) )
    elif (varname.upper() == "NOX"):
        gdata =  ( numpy.array(f.variables['IJ_AVG_S__NO']) +
                   numpy.array(f.variables['IJ_AVG_S__NO2']) )
    elif (varname.upper() == "ETNO3_C2H6"):
        gdata =  ( numpy.array(f.variables['IJ_AVG_S__ETNO3']) /
                   numpy.array(f.variables['IJ_AVG_S__C2H6']) )
        conv = 1e0
        unit = 'ppt / ppt'
    elif (varname.upper() == "PRNO3_C3H8"):
        gdata =  ( (numpy.array(f.variables['IJ_AVG_S__IPRNO3']) +
                    numpy.array(f.variables['IJ_AVG_S__NPRNO3'])) /
                   numpy.array(f.variables['IJ_AVG_S__C3H8']) )
        conv = 1e0
        unit = 'ppt / ppt'
    else:
        gdata = numpy.array(f.variables['IJ_AVG_S__'+varname.upper()])
    
    # apply unit conversion of required
    gdata = gdata * conv    
    
    return {"data":gdata, "lon":glon, "lat":glat, "alt":galt, "bh":gbxht,
            "unit":unit}

def read_gc_bpch(varname,filename,ppt=True):
    
    """Read a GEOS-Chem BPCH file created directly from model and extract
       lon, lat, altitude (from boxheight), and data field.
       
       Default is to convert to ppt from ppb but this can be overwritten
       with ppb=False"""

    args = PNC("--format=bpch,nogroup=('IJ-AVG-$','BXHGHT-$',)",filename)
    f = args.ifiles[0]

    # Essential variables first
    glon = f.variables['longitude']
    glat = f.variables['latitude']
    gbxht = f.variables['BXHEIGHT']
    
    # use boxheight variable to calculate altitude
    galt = numpy.zeros(gbxht.shape)
    nlev=galt.shape[1]
    for ll in range(nlev):
        # alt is cumulative sum, convert to m
        galt[:,ll,:,:] = numpy.sum(gbxht[:,0:ll,:,:],axis=1)*1e-3

    # conv is factor to convert from ppb if needed
    # defaults to ppt but can be set otherwise, as in special cases
    if ppt:
        conv = 1e3
        unit = 'ppt'
    else:
        conv = 1e0

    # Get data - deal with "special" fields first
    if (varname.upper() == "C1-C3_RONO2"):
        gdata =  ( numpy.array(f.variables['MENO3']) +
                   numpy.array(f.variables['ETNO3']) +
                   numpy.array(f.variables['IPRNO3']) +
                   numpy.array(f.variables['NPRNO3']) )
    elif (varname.upper() == "PRNO3"):
        gdata =  ( numpy.array(f.variables['IPRNO3']) +
                   numpy.array(f.variables['NPRNO3']) )
    elif (varname.upper() == "NOX"):
        gdata =  ( numpy.array(f.variables['NO']) +
                   numpy.array(f.variables['NO2']) )
    elif (varname.upper() == "ETNO3_C2H6"):
        gdata =  ( numpy.array(f.variables['ETNO3']) /
                   numpy.array(f.variables['C2H6']) )
        conv = 1e0
        unit = 'ppt / ppt'
    elif (varname.upper() == "PRNO3_C3H8"):
        gdata =  ( (numpy.array(f.variables['IPRNO3']) +
                    numpy.array(f.variables['NPRNO3'])) /
                   numpy.array(f.variables['C3H8']) )
        conv = 1e0
        unit = 'ppt / ppt'
    else:
        gdata = numpy.array(f.variables[varname.upper()])
    
    # apply unit conversion of required
    gdata = gdata * conv    
    
    return {"data":gdata, "lon":glon, "lat":glat, "alt":galt, "bh":gbxht,
            "unit":unit}

def extract_gc_2d_lat_lon(gcdatablock,lev=None,altrange=None):
    
    """Extracts a 2-D lat-lon slice of GEOS-Chem output for either:
       - a given level (lev=X]) or
       - a given altitude range (altrange=[X,Y])
      
      Default behaviour is to extract surface level (lev=0)"""
    
    gcdata=gcdatablock["data"]
    gcalt =gcdatablock["alt"]

    if hasattr(lev, "__len__") and (altrange is not None):
        print ("lev takes precedence over altrange; using lev!")
        altrange=None
    elif (lev is None) and (altrange is None):
        print ("no altitude specified; using surface!")
        lev = 0

    if (altrange is not None):

        # extract altitude range
        nlev=gcalt.shape[1]
        first = True
        for ll in range(nlev):            
            if numpy.mean(gcalt[:,ll,:,:]) <= altrange[0]:
                lev1=ll
            if (numpy.mean(gcalt[:,ll,:,:]) >= altrange[-1]) & (first):
                lev2=ll
                first = False

        lev = [lev1,lev2]

    # Extract appropriate levels
    gcdata = gcdata[0,lev,:,:]
    
    # average data over levels if required
    if hasattr(lev, "__len__"):
        gcdata = gcdata.mean(axis=0)

        
    return gcdata        

def extract_gc_1d_alt(gcdatablock,lonrange=None,latrange=None,
                      binned=False,dalt=1.0):
    
    """Extracts a 1-D altitude alice of GEOS-Chem output for a given:
       - Latitude range (latrange=[X,Y]) and
       - Longitude range (lonrange=[X,Y])
       
       Altitude array is also averaged over the same region.
       
       Default behaviour is to average over entire globe"""
       
    # convert to numpy arrays
    gcdata=gcdatablock["data"]
    gcalt =gcdatablock["alt"]
    gclon =gcdatablock["lon"]
    gclat =gcdatablock["lat"]
       
    if (lonrange is None):
        print ("No longitudes specified: using global!")
        lonrange = [-180,180]
    if (latrange is None):
        print ("No latitudes specified: using global!")
        latrange = [-90,90]
    if len(lonrange) != 2:
        print ("Specify longitudes as [lon1,lon2]")
        return
    if len(latrange) != 2:
        print ("Specify latitudes as [lat1,lat2]")
        return
        
    # convert to numpy arrays
    lonrange=numpy.array(lonrange)
    latrange=numpy.array(latrange)
       
    # find longitude & latitude bounds
    ind1 = [(gclon >= lonrange.min()) & (gclon <= lonrange.max())]
    ind2 = [(gclat >= latrange.min()) & (gclat <= latrange.max())]
    
    # subselect data & alt
    gcdata = gcdata[:,:,:,ind1[0]]
    gcdata = gcdata[:,:,ind2[0],:]
    gcalt  = gcalt [:,:,:,ind1[0]]
    gcalt  = gcalt[:,:,ind2[0],:]
    
    # at this point, average altitude array so it is 1-d
    gcalt  = numpy.mean(gcalt, axis=(0,2,3))    

    # Calculate median, IQR over 1-km altitude bins
    if (binned):
        altbins = numpy.arange(gcalt.max()+dalt,step=dalt)
        nalt = len(altbins)-1
        gcout = numpy.zeros((3,nalt))
        for ll in range(nalt):
            altl = altbins[ll]
            gcdtmp = gcdata[:,(gcalt > altl) & (gcalt < altl+dalt),:,:]
            if gcdtmp.shape[1] > 0:
                gcout[:,ll] = numpy.percentile(gcdtmp,(25,50,75))
            else:
                gcout[:,ll] = numpy.nan, numpy.nan, numpy.nan
        gcalt = altbins[:-1]+dalt/2.
    else:    
        # take mean over lat, lon, & time dimensions
        gcout = numpy.percentile(gcdata,(25,50,75),axis=(0,2,3))
            
    return gcout, gcalt
           
