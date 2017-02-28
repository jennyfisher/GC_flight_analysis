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
from scipy.stats import binned_statistic

# for testing
import time
#

# My packages
from read_airborne_data import *
from read_gc_data import *

matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['font.sans-serif'] = 'Arial' #- for illustrator...
matplotlib.rcParams['font.sans-serif'] = 'Bitstream Vera Sans'

def boxes(airlon,airlat,pairdata):
    
    # for now follow zitely, 20° lon x 10°lat bins
    lonedge = numpy.arange(-180,181,5)
    latedge = numpy.arange(-90,91,4)
    
    newdata = numpy.zeros([len(lonedge),len(latedge)])
    count   = numpy.zeros([len(lonedge),len(latedge)])
    newlon  = numpy.zeros([len(lonedge),len(latedge)])
    newlat  = numpy.zeros([len(lonedge),len(latedge)])
    
    for nn in range(len(pairdata)):
        ii = (numpy.where(airlon[nn] <= lonedge))[0][0]
        jj = (numpy.where(airlat[nn] <= latedge))[0][0]
        newdata[ii,jj] = newdata[ii,jj] + pairdata[nn]
        count[ii,jj] = count[ii,jj] + 1
        
        newlon[ii,jj] = lonedge[ii]+2.5
        newlat[ii,jj] = latedge[jj]+2
        
    newdata = newdata / count

    newdata = newdata.flatten()
    newlon = newlon.flatten()
    newlat = newlat.flatten()
        
    #newdata = (numpy.histogram2d(airlon, airlat, [lonedge, latedge], weights=pairdata) /
    #            numpy.histogram2d(airlon, airlat, [lonedge, latedge], weights=pairdata) )
    
    return newlon, newlat, newdata

    
    

def map_gc(varname,filename,savefig=False,airdata=None,airname='',
           usebox=True,lev=None,altrange=None,
           mindata=None,maxdata=None,outline=True,
           figname=None,title=''):
    
    """
    This function makes a map of a given GEOS-Chem variable at either a
    given level (lev=X) or over a given range of altitudes (altrange=[X,Y]).
    
    For now, map is global. Aircraft data can be overplotted using the
    keyword argument airdata, where airdata should include variables "lon",
    "lat", and "data", with bad data marked as numpy.nan).
    
    If more than one filename is specified, use the mean over all files.
    
    To save the figure, use savefig=True.
    To specify figure name set a value for figname.
    
    To specify colorbar range, use mindata=X and maxdata=Y keywords."""
    
    # Read data & extract variables
    first = True
    
    # Average over multiple files
    for f in filename:
        GC = read_gc_nc(varname,f)    
        tdata = extract_gc_2d_lat_lon(GC,lev=lev,altrange=altrange)
        if first:
            data = tdata
            first = False
        else:
            data = data + tdata
    data = data / len(filename)
    lon = GC["lon"][:]
    lat = GC["lat"][:]
    
    
    # Get min and max data values from GC if needed
    if mindata is None:
        mindata = numpy.min(data)
    if maxdata is None:
        maxdata = numpy.max(data)
    
    # repeat last data column to avoid white space    
    data,lon = addcyclic(data, lon)    
    
    # transform lon/lat to coordinate grid
    lon,lat = numpy.meshgrid(lon,lat)
    
    # Set up the map
    map=Basemap(projection='robin',lon_0=-130,lat_0=0)
    f=pyplot.figure()
    
    # plot the data
    col=map.pcolormesh(lon,lat,data,latlon=True,
                       cmap=pyplot.get_cmap('YlOrRd'),
                       vmin=mindata,vmax=maxdata)
    
    # improve the map
    map.drawcoastlines()
    map.drawcountries()

    # add a color bar
    cb = map.colorbar(col, "bottom")
    cb.set_label(varname.upper()+', ppt')

    # get the aircraft data to overplot
    if airdata is not None:
        airlon=airdata["lon"][:]
        airlat=airdata["lat"][:]
        pairdata =airdata["data"][:]
        
        # get rid of NaN values for plotting
        airlon = airlon[~numpy.isnan(pairdata)]
        airlat = airlat[~numpy.isnan(pairdata)]
        pairdata = pairdata[~numpy.isnan(pairdata)]
        
        # get rid of masked data for plotting, etc.
        if numpy.ma.is_masked(pairdata):
            airlon   = airlon[~pairdata.mask]
            airlat   = airlat[~pairdata.mask]
            pairdata = pairdata[~pairdata.mask]
        
        if usebox:
            airlon, airlat, pairdata = boxes(airlon,airlat,pairdata)
            airlon = airlon[~numpy.isnan(pairdata)]
            airlat = airlat[~numpy.isnan(pairdata)]
            pairdata = pairdata[~numpy.isnan(pairdata)]
        
        # convert to map coordinates
        airlon,airlat = map(airlon,airlat)
        # add data locations to map
        if outline:
            map.scatter(airlon,airlat,c="white")

        # add data points to map    
        map.scatter(airlon,airlat,c=pairdata,
                    cmap=pyplot.get_cmap('YlOrRd'),
                    vmin=mindata,vmax=maxdata,
                    edgecolors='none')
                
    if ( savefig ):
        # filename
        if figname is None:
            figname=airname+varname
            if altrange is not None:
                figname=figname+'_'+str(altrange[0])+'-'+str(altrange[1])
            else:
                figname=figname+'_L'+str(lev)
        pyplot.savefig('img/'+figname+'.pdf')
        pyplot.close()
    else:
        pyplot.show()
    
    return f

def profile_gc(varname,filename,savefig=False,airdata=None,airname='',
               latrange=None,lonrange=None,altrange=None,dalt=1.0,
               mindata=None,maxdata=None,outline=True,
               figname=None,title='',filename2=[],ftype="netcdf"):
    
    """
    This function makes an altitude profile of a given GEOS-Chem variable
    averaged globally or over a region (specified with latrange, lonrange.
    
    If more than one filename is specified, use the mean over all files.    

    To show a second set of model files, use filename2=[file1,file2,...].
    
    To specify altitude bounds for profile, use altrange=[X,Y].
    
    To save the figure, use savefig=True.
    To specify figure name set a value for figname.
    
    To specify data axis range, use mindata=X and maxdata=Y keywords."""
    
    # Average over multiple files
    first = True
    for f in filename:
        if ftype == "netcdf":
            GC = read_gc_nc(varname,f)
        elif ftype == "bpch":
            GC = read_gc_bpch(varname,f)
        else:
            sys.exit("file type not allowed!")

        unit = GC["unit"]
        tdata, talt = extract_gc_1d_alt(GC,lonrange=lonrange,latrange=latrange,
                                  binned=True,dalt=dalt)
        if first:
            data = tdata
            alt = talt
            first = False
        else:
            # sometimes something weird happens with altitude dimension...
            if data.shape[-1] < tdata.shape[-1]:
                data = data + tdata[:,:data.shape[-1]]
                alt = alt + talt[:alt.shape[-1]]
            elif data.shape[-1] > tdata.shape[-1]:
                data[:,:tdata.shape[-1]] = data[:,:tdata.shape[-1]] + tdata
                alt[:talt.shape[-1]] = alt[:talt.shape[-1]] + talt
            else:
                data = data + tdata
                alt = alt + talt
    data = data / len(filename)    
    alt = alt / len(filename)    
        
    # KLUDGE! fix issue where one alt bin undefined
    if numpy.isnan(data[1,8]):
        data[:,8] = numpy.mean([data[:,7],data[:,9]],axis=0)
    
    # Second model profile?
    if len(filename2) > 0:
        first = True
        for f in filename2:
            if ftype == "netcdf":
                GC2 = read_gc_nc(varname,f)
            elif ftype == "bpch":
                GC2 = read_gc_bpch(varname,f)
            else:
                sys.exit("file type not allowed!")
            unit = GC2["unit"]
            tdata, talt = extract_gc_1d_alt(GC2,lonrange=lonrange,latrange=latrange,
                                      binned=True,dalt=dalt)
            if first:
                data2 = tdata
                alt2 = talt
                first = False
            else:
                # sometimes something weird happens with altitude dimension...
                if data2.shape[-1] < tdata.shape[-1]:
                    data2 = data2 + tdata[:,:data2.shape[-1]]
                    alt2 = alt2 + talt[:alt2.shape[-1]]
                elif data2.shape[-1] > tdata.shape[-1]:
                    data2[:,:tdata.shape[-1]] = data2[:,:tdata.shape[-1]] + tdata
                    alt2[:talt.shape[-1]] = alt2[:talt.shape[-1]] + talt
                else:
                    data2 = data2 + tdata
                    alt2 = alt2 + talt
        data2 = data2 / len(filename2)    
        alt2 = alt2 / len(filename2)    

        # KLUDGE! fix issue where one alt bin undefined
        if numpy.isnan(data2[1,8]):
            data2[:,8] = numpy.mean([data2[:,7],data2[:,9]],axis=0)
    
    # Get min and max data values from GC if needed
    if mindata is None:
        mindata = numpy.nanmin(data[0,:])*.9
    if maxdata is None:
        maxdata = numpy.nanmax(data[2,:])*1.1
    
    # Get altitude range if needed
    if altrange is None:
        altrange=[alt.min(),alt.max()]

    f=pyplot.figure()
    
    # plot the data
    pyplot.plot(data[1,:],alt,'-r',lw=2)
    pyplot.fill_betweenx(alt,data[0,:],data[2,:],color='red',alpha=0.5)
    pyplot.xlim([mindata,maxdata])
    pyplot.ylim(altrange)
    pyplot.minorticks_on()
    pyplot.ylabel('Altitude (km)')
    pyplot.xlabel(varname+' ('+unit+')')
    pyplot.title(title)
    
    # plot a second model profile?
    if len(filename2) > 0:
        pyplot.plot(data2[1,:],alt2,'-b',lw=2)
        pyplot.fill_betweenx(alt2,data2[0,:],data2[2,:],color='blue',alpha=0.5)
    
    # get the aircraft data to overplot
    if airdata is not None:
        airlon=airdata["lon"][:]
        airlat=airdata["lat"][:]
        airalt=airdata["alt"][:]
        pairdata =airdata["data"][:]
        
        # get rid of NaN values for plotting
        airlon = airlon[~numpy.isnan(pairdata)]
        airlat = airlat[~numpy.isnan(pairdata)]
        airalt = airalt[~numpy.isnan(pairdata)]
        pairdata = pairdata[~numpy.isnan(pairdata)]

        # TODO: figure out how to treat masked values properly! but for now
        airlon = airlon[pairdata > -7777]
        airlat = airlat[pairdata > -7777]
        airalt = airalt[pairdata > -7777]
        pairdata = pairdata[pairdata > -7777]

        # bin data by altitude
        pairmed = binned_statistic(airalt,pairdata,'median',
                                   bins=numpy.arange(altrange[1]+dalt,step=dalt))
        pair25 = binned_statistic(airalt,pairdata,
                                  statistic=lambda y: numpy.percentile(y, 25),
                                  bins=numpy.arange(altrange[1]+dalt,step=dalt))
        pair75 = binned_statistic(airalt,pairdata,
                                  statistic=lambda y: numpy.percentile(y, 75),
                                  bins=numpy.arange(altrange[1]+dalt,step=dalt))
        pairalt = pairmed.bin_edges[:-1]+dalt/2.


        # plot...
        pyplot.plot(pairmed.statistic,pairalt,color='black')
        pyplot.fill_betweenx(pairalt,pair25.statistic,pair75.statistic,
                             color='gray',alpha=0.5)
                             
    if ( savefig ):
        # filename
        if figname is None:
            figname=airname+varname
            if altrange is not None:
                figname=figname+'_'+str(altrange[0])+'-'+str(altrange[1])
            else:
                figname=figname+'_L'+str(lev)
        pyplot.savefig('img/'+figname+'.pdf')
        pyplot.close()
    else:
        pyplot.show()
    
    return f
