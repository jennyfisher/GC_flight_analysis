# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 18:10:29 2016

@author: jennyf
"""

import numpy

def read_ict(fname):
    
    # Read an ICARTT format ascii data file and return data array
    # fname:    filename
    # varnames: return variable names?

    # Read data file
    file = open(fname,'r')

    # Read first header line to get header count
    header = file.readline()
    n_hdr  = long((header.split(', '))[0])

    # Read 12th line to get missing values
    for i in range(11):
        missing_values = file.readline()

    # Close file
    file.close()
    
    # Now that we know number of header rows, we can read the data
    data = numpy.genfromtxt(fname,
                            skip_header=n_hdr-1,        # skip header lines
                            delimiter=', ',             # file delimiter
                            names=True,                 # get column names
#                           missing_values='-9999.0',   # find missing values
                            missing_values=missing_values,   # find missing values
                            filling_values=numpy.nan,   # set missing to nan
                            usemask=True)               # needed for missing_values


    return(data)

#------------------------------------------------------------------------------
## SEAC4RS


def gcname_to_seacname(argument):
    
    """ Converts from GEOS-Chem species name to SEAC4RS data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MeNO3": "MeONO2_WAS",
        "MeNO3".lower(): "MeONO2_WAS",
        "MeNO3".upper(): "MeONO2_WAS",
        "EtNO3": "EtONO2_WAS",
        "EtNO3".lower(): "EtONO2_WAS",
        "EtNO3".upper(): "EtONO2_WAS",
        "iPrNO3": "iPrONO2_WAS",
        "iPrNO3".lower(): "iPrONO2_WAS",
        "iPrNO3".upper(): "iPrONO2_WAS",
        "nPrNO3": "nPrONO2_WAS",
        "nPrNO3".lower(): "nPrONO2_WAS",
        "nPrNO3".upper(): "nPrONO2_WAS",
        "C2H6": "Ethane_WAS",
        "c2h6": "Ethane_WAS",
        "C3H8": "Propane_WAS",
        "c3h8": "Propane_WAS"
    }
    return switcher.get(argument, origname)

def extract_seac4rs(obs,varname,
                  altrange=None,latrange=None,lonrange=None,dayrange=None):
    
    """ Extracts SEAC4RS data & lon/lat/alt fields for a given species.
    Requires:
       obs = full SEAC4RS data array
       varname = name of variable in SEAC4RS data file
       
   Optional:
       altrange = range of altitude to limit data extracted (default=all)"""
    
    # Get essential fields
    slon = obs["LONGITUDE"]
    slat = obs["LATITUDE"]
    salt = obs["ALTP"]
    sday = obs["JDAY"]
    
    # Fix SEAC4RS longitude ([-180,180] from [0,360])
    slon[slon > 180] = slon[slon > 180] - 360.
    
    
    # Get data - deal with "special" fields first
    if (varname.upper() == "C1-C3_RONO2"):
        sdata = ( obs["MeONO2_WAS"]    +
                  obs["EtONO2_WAS"]   +
                  obs["iPrONO2_WAS"]+
                  obs["nPrONO2_WAS"] )
        sunit = 'ppt'
    elif (varname.upper() == "PRNO3"):
        sdata = ( obs["iPrONO2_WAS"]+
                  obs["nPrONO2_WAS"] )
        sunit = 'ppt'
    elif (varname.upper() == "ETNO3_C2H6"):
        sdata = ( obs["EtONO2_WAS"] /
                  obs["Ethane_WAS"] )
        sunit = 'ppt / ppt'
    elif (varname.upper() == "PRNO3_C3H8"):
        sdata = ( (obs["iPrONO2_WAS"]+
                   obs["nPrONO2_WAS"]) /
                  obs["Propane_WAS"] )
        sunit = 'ppt / ppt'
    else:
        sdata = obs[varname]
        sunit = 'ppt'
    
    # Restrict altitude range, if required
    if altrange is not None:
        ind = numpy.where( (salt >= altrange[0]) &
                           (salt <= altrange[1]) )
        slon = slon[ind]
        slat = slat[ind]
        salt = salt[ind]
        sday = sday[ind]
        sdata =sdata[ind]
    # Restrict latitude range, if required
    if latrange is not None:
        ind = numpy.where( (slat >= latrange[0]) &
                           (slat <= latrange[1]) )
        slon = slon[ind]
        slat = slat[ind]
        salt = salt[ind]
        sday = sday[ind]
        sdata =sdata[ind]
    # Restrict longitude range, if required
    if lonrange is not None:
        ind = numpy.where( (slon >= lonrange[0]) &
                           (slon <= lonrange[1]) )
        slon = slon[ind]
        slat = slat[ind]
        salt = salt[ind]
        sday = sday[ind]
        sdata =sdata[ind]
    # Restrict date range, if required
    if dayrange is not None:
        ind = numpy.where( (sday >= dayrange[0]) &
                           (sday <= dayrange[1]) )
        slon = slon[ind]
        slat = slat[ind]
        salt = salt[ind]
        sday = sday[ind]
        sdata =sdata[ind]
        
    # Screen out ULOD, LLOD values
    if sdata.min() < -111111.:
        ind = numpy.where(sdata > -111111.)
        slon = slon[ind]
        slat = slat[ind]
        salt = salt[ind]
        sday = sday[ind]
        sdata =sdata[ind]

    return {"data":sdata, "lon":slon, "lat":slat, "alt":salt,"jday":sday,
            "unit":sunit}

    
#------------------------------------------------------------------------------
## HIPPO
def read_hippo(fname):
    
    """Reads the HIPPO data files (.tbl) and returns array of all data"""
    
    # Open data file
    #file = open(fname,'r')
    
    # Read file, including variable names
    data = numpy.genfromtxt(fname,
                            delimiter=' ',
                            names=True)
                            
    return(data)
    
def gcname_to_hname(argument):
    
    """ Converts from GEOS-Chem species name to HIPPO data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MeNO3": "MeONO2_AW",
        "MeNO3".lower(): "MeONO2_AW",
        "MeNO3".upper(): "MeONO2_AW",
        "EtNO3": "EthONO2_AW",
        "EtNO3".lower(): "EthONO2_AW",
        "EtNO3".upper(): "EthONO2_AW",
        "iPrNO3": "i_PropONO2_AW",
        "iPrNO3".lower(): "i_PropONO2_AW",
        "iPrNO3".upper(): "i_PropONO2_AW",
        "nPrNO3": "n_PropONO2_AW",
        "nPrNO3".lower(): "n_PropONO2_AW",
        "nPrNO3".upper(): "n_PropONO2_AW",
        "C2H6": "Ethane_AW",
        "c2h6": "Ethane_AW",
        "C3H8": "Propane",
        "c3h8": "Propane"
    }
    return switcher.get(argument, origname)

def extract_hippo(obs,varname,
                  altrange=None,latrange=None,lonrange=None,dayrange=None):
    
    """ Extracts HIPPO data & lon/lat/alt fields for a given species.
    Requires:
       obs = full HIPPO data array
       varname = name of variable in HIPPO data file
       
   Optional:
       altrange = range of altitude to limit data extracted (default=all)"""
    
    # Get essential fields
    hlon = obs["GGLON"]
    hlat = obs["GGLAT"]
    halt = obs["GGALT"]*1e-3 # convert from m to km
    hday = obs["jd"]
    year = obs["Year"]
    
    # fix hday so it is doy (year independent)
    hday = hday - (year-2009)*365.
         
    # Get data - deal with "special" fields first
    if (varname.upper() == "C1-C3_RONO2"):
        hdata = ( obs["MeONO2_AW"]    +
                  obs["EthONO2_AW"]   +
                  obs["i_PropONO2_AW"]+
                  obs["n_PropONO2_AW"] )
        hunit = 'ppt'
    elif (varname.upper() == "PRNO3"):
        hdata = ( obs["i_PropONO2_AW"]+
                  obs["n_PropONO2_AW"] )
        hunit = 'ppt'
    elif (varname.upper() == "ETNO3_C2H6"):
        hdata = ( obs["EthONO2_AW"] /
                  obs["Ethane_AW"] )
        hunit = 'ppt'
    elif (varname.upper() == "PRNO3_C3H8"):
        hdata = ( (obs["i_PropONO2_AW"]+
                   obs["n_PropONO2_AW"]) /
                  obs["Propane"] )
        hunit = 'ppt / ppt'
    else:
        hdata = obs[varname]
        hunit = 'ppt'
    
    # Restrict altitude range, if required
    if altrange is not None:
        ind = numpy.where( (halt >= altrange[0]) &
                           (halt <= altrange[1]) )
        hlon = hlon[ind]
        hlat = hlat[ind]
        halt = halt[ind]
        hday = hday[ind]
        hdata = hdata[ind]
    # Restrict latitude range, if required
    if latrange is not None:
        ind = numpy.where( (hlat >= latrange[0]) &
                           (hlat <= latrange[1]) )
        hlon = hlon[ind]
        hlat = hlat[ind]
        halt = halt[ind]
        hday = hday[ind]
        hdata = hdata[ind]
    # Restrict longitude range, if required
    if lonrange is not None:

# To do: straddle date line        
#        # straddle date line?
#        DL = False
#        if (lonrange[0] > 0) and (lonrange[1] < 0):
#            hlon[hlon < 0] = hlon[hlon < 0]+360
#            lonrange[11] = lonrange[11]+360
#            DL = True
#        
        ind = numpy.where( (hlon >= lonrange[0]) &
                           (hlon <= lonrange[1]) )
#                           
#        if DL:
#            hlon[hlon > 180] = hlon[hlon > 180]-360
                           
        hlon = hlon[ind]
        hlat = hlat[ind]
        halt = halt[ind]
        hday = hday[ind]
        hdata = hdata[ind]
    # Restrict date range, if required
    if dayrange is not None:
        ind = numpy.where( (hday >= dayrange[0]) &
                           (hday <= dayrange[1]) )
        hlon = hlon[ind]
        hlat = hlat[ind]
        halt = halt[ind]
        hday = hday[ind]
        hdata = hdata[ind]
        
    return {"data":hdata, "lon":hlon, "lat":hlat, "alt":halt,"jday":hday,
            "unit":hunit}
    
#------------------------------------------------------------------------------
## combine multiple
def airdata_combine(data1,data2):
    
    """Combine two airborne datasets after extracting"""

    lon = numpy.ma.concatenate((data1["lon"],data2["lon"]))
    lat = numpy.ma.concatenate((data1["lat"],data2["lat"]))
    alt = numpy.ma.concatenate((data1["alt"],data2["alt"]))
    jday = numpy.ma.concatenate((data1["jday"],data2["jday"]))
    data = numpy.ma.concatenate((data1["data"],data2["data"]))
    
    return {"data":data, "lon":lon, "lat":lat, "alt":alt,"jday":jday}