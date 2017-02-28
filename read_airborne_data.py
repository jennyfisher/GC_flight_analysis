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
    n_hdr  = int((header.split(', '))[0])

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
## FRAPPE
def gcname_to_frappname(argument):
    
    """ Converts from GEOS-Chem species name to FRAPPE data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "MeONO2_WAS",
        "ETNO3": "EtONO2_WAS",
        "IPRNO3": "iPrONO2_WAS",
        "NPRNO3": "nPrONO2_WAS",
        "ANS": "ANs_LIF",
        "C2H6": "Ethane_WAS",
        "C3H8": "Propane_WAS",
        "ALD2": "Acetaldehyde_MixingRatio",
        "ACET": "AcetonePropanal_MixingRatio",
        "NO2": "NO2_MixingRatio",
        "PAN": "PAN",
        "O3": "O3_MixingRatio",
    }
    return switcher.get(argument.upper(), origname)

## SEAC4RS
def gcname_to_seacname(argument):
    
    """ Converts from GEOS-Chem species name to SEAC4RS data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "MeONO2_WAS",
        "ETNO3": "EtONO2_WAS",
        "IPRNO3": "iPrONO2_WAS",
        "NPRNO3": "nPrONO2_WAS",
        "ANS": "ANs_TDLIF",
        "C2H6": "Ethane_WAS",
        "C3H8": "Propane_WAS",
        "ALD2": "Acetaldehyde",
        "ACET": "Acetone_Propanal",
        "NO2": "NO2_TDLIF",
        "PAN": "PAN_GTCIMS",
        "O3": "O3_ESRL",
    }
    return switcher.get(argument.upper(), origname)

## DC3
def gcname_to_dc3name(argument):
    
    """ Converts from GEOS-Chem species name to DC-3 data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "MeONO2_WAS",
        "ETNO3": "EtONO2_WAS",
        "IPRNO3": "iPrONO2_WAS",
        "NPRNO3": "nPrONO2_WAS",
        "ANS": "ANs_TDLIF",
        "C2H6": "Ethane_WAS",
        "C3H8": "Propane_WAS",
        "ALD2": "Acetaldehyde_PTRMS",
        "ACET": "Acetone_Propanal_PTRMS",
        "NO2": "NO2_TDLIF",
        "PAN": "PAN_GTCIMS",
        "O3": "O3_ESRL",
    }
    return switcher.get(argument.upper(), origname)

## CalNex
def gcname_to_calname(argument):
    
    """ Converts from GEOS-Chem species name to CalNex data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "MeONO2",
        "ETNO3": "EtONO2",
        "IPRNO3": "i_PrONO2",
        "NPRNO3": "n_PrONO2",
        "C2H6": "Ethane",
        "C3H8": "Propane",
        "ALD2": "Acetaldehyde",
        "ACET": "Acetone",
        "NO2": "NO2_ppbv",
        "PAN": "PAN_ppbv",
        "O3": "O3_ppbv",
    }
    return switcher.get(argument.upper(), origname)

## ARCTAS
def gcname_to_arcname(argument):
    
    """ Converts from GEOS-Chem species name to ARCTAS data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "MeONO2",
        "ETNO3": "EtONO2",
        "IPRNO3": "iPrONO2",
        "NPRNO3": "nPrONO2",
        "ANS": "ANs_UCB",
        "C2H6": "Ethane",
        "C3H8": "Propane",
        "ALD2": "Acetaldehyde_TOGA",
        "ACET": "Acetone_TOGA",
        "NO2": "NO2_NCAR",
        "PAN": "GT_PAN",
        "O3": "O3",
    }
    return switcher.get(argument.upper(), origname)

## TexAQS
def gcname_to_texname(argument):
    
    """ Converts from GEOS-Chem species name to TexAQS data field name"""
    
    # Default return original name if not found (may be special case)
    origname=argument
    switcher = {
        "MENO3": "Methyl_nitrate_Atlas",
        "ETNO3": "Ethyl_nitrate_Atlas",
        "IPRNO3": "Isopropyl_nitrate_Atlas",
        "NPRNO3": "n_Propylnitrate_Atlas",
        "NO2": "NO2_ppbv_Ryerson",
        "PAN": "PAN_pptv_Flocke",
    }
    return switcher.get(argument.upper(), origname)


def extract_seac4rs(obs,varname,
                  altrange=None,latrange=None,lonrange=None,dayrange=None,altunit="km",
		  lonname = "LONGITUDE",latname="LATITUDE",altname="ALTP",dayname="JDAY"):
    
    """ Extracts SEAC4RS data & lon/lat/alt fields for a given species.
    Requires:
       obs = full SEAC4RS data array
       varname = name of variable in SEAC4RS data file
       
   Optional:
       altrange = range of altitude to limit data extracted (default=all)"""
    
    # Get essential fields
    #slon = obs["LONGITUDE"]
    #slat = obs["LATITUDE"]
    #salt = obs["ALTP"]
    #sday = obs["JDAY"]
    slon = obs[lonname]
    slat = obs[latname]
    salt = obs[altname]
    sday = obs[dayname]

    # Fix altitudes - convert from m to km
    if altunit=="m": salt = salt*1e-3
    
    # Fix SEAC4RS longitude ([-180,180] from [0,360])
    slon[slon > 180] = slon[slon > 180] - 360.
    
    sunit = 'ppt'

    # Get data - deal with "special" fields first
    if (varname.upper() == "C1-C3_RONO2"):
        sdata = ( obs["MeONO2_WAS"]    +
                  obs["EtONO2_WAS"]   +
                  obs["iPrONO2_WAS"]+
                  obs["nPrONO2_WAS"] )
    elif (varname.upper() == "NOX"):
        sdata = ( obs["NO2_ESRL"] +
                  obs["NO_ESRL"]  )*1e3
    elif (varname.upper() == "PRNO3"):
        sdata = ( obs["iPrONO2_WAS"]+
                  obs["nPrONO2_WAS"] )
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
        "MENO3": "MeONO2_AW",
        "ETNO3": "EthONO2_AW",
        "IPRNO3": "i_PropONO2_AW",
        "NPRNO3": "n_PropONO2_AW",
        "C2H6": "Ethane_AW",
        "C3H8": "Propane",
    }
    return switcher.get(argument.upper(), origname)

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
