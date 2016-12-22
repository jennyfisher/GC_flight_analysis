# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 13:08:18 2016

@author: jennyf
"""

# Standard packages
import numpy
import matplotlib.pyplot as pyplot
import matplotlib
from mpl_toolkits.basemap import Basemap, addcyclic

# My packages
from read_gc_data import *

matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['font.sans-serif'] = 'Arial' #- for illustrator...
matplotlib.rcParams['font.sans-serif'] = 'Bitstream Vera Sans'

def map_gc_diff(varname,filename1,filenam2, percent=False, savefig=False,
           lev=None,altrange=None,cbar=None,
           figname=None):
    
    """
    Equivalent to CTM_PLOTDIFF basic functionality
    
    If more than one filename is specified, use the mean over all files.
    
    To save the figure, use savefig=True.
    To specify figure name set a value for figname.
    
    To specify colorbar range, use mindata=X and maxdata=Y keywords."""
    
    # Read data & extract variables
    
    # Average over multiple files
    first = True
    for f in filename1:
        print "Old Data from file: "+f
        GC = read_gc_nc(varname,f)    
        tdata = extract_gc_2d_lat_lon(GC,lev=lev,altrange=altrange)
        if first:
            data = tdata
            first = False
        else:
            data = data + tdata
    data1 = data / len(filename1)
    
    first = True
    for f in filename2:
        print "New Data from file: "+f
        GC = read_gc_nc(varname,f)    
        tdata = extract_gc_2d_lat_lon(GC,lev=lev,altrange=altrange)
        if first:
            data = tdata
            first = False
        else:
            data = data + tdata
    data2 = data / len(filename2)

    
    lon = GC["lon"][:]
    lat = GC["lat"][:]
    
    diff = data2 - data1
    pdiff = numpy.zeros(diff.shape)
    pdiff[abs(data1) > 0] = 100. * diff[abs(data1) > 0] / data1[abs(data1) > 0]
    
    if percent:
        data = pdiff
    else:
        data = diff
    
    # Get min and max data values from GC if needed
    if cbar is None:
        maxdata = numpy.max(abs(data))
        mindata = -1.0 * maxdata
    elif not hasattr(cbar, "__len__"): 
        maxdata = cbar
        mindata = -1.0 * maxdata
    elif len(cbar)==2:
        maxdata = cbar.max()
        mindata = cbar.min()
    else:
        print "Specify cbar as 1 or 2 values"
        return

    
    # repeat last data column to avoid white space    
    data,lon = addcyclic(data, lon)    
    
    # transform lon/lat to coordinate grid
    lon,lat = numpy.meshgrid(lon,lat)
    
    # Set up the map
    map=Basemap(projection='robin',lon_0=-130,lat_0=0)
    f=pyplot.figure()
    
    # plot the data
    col=map.pcolormesh(lon,lat,data,latlon=True,
                       cmap=pyplot.get_cmap('coolwarm'),
                       vmin=mindata,vmax=maxdata)
    
    # improve the map
    map.drawcoastlines()
    map.drawcountries()

    # add a color bar
    cb = map.colorbar(col, "bottom")
    if percent:
        cb.set_label(varname.upper()+', %')
    else:
        cb.set_label(varname.upper()+', ppt')
                
    if ( savefig ):
        # filename
        if figname is None:
            figname=varname
            if altrange is not None:
                figname=figname+'_'+str(altrange[0])+'-'+str(altrange[1])
            else:
                figname=figname+'_L'+str(lev)
        pyplot.savefig('img/'+figname+'.pdf')
        pyplot.close()
    else:
        pyplot.show()
    
    return f


## Test some plots
month=range(1,13)
filename1=[]
filename2=[]
for m in month:
    filename1.append('RF/gc.2013'+str(m).zfill(2)+'01.nc')
    filename2.append('RG/gc.2013'+str(m).zfill(2)+'01.nc')
    
v='PAN'    
altrange=[0,2]
percent=False


# Call it
GC=map_gc_diff(v,filename1,filename2,savefig=False,
          altrange=altrange,percent=percent,cbar=10)

