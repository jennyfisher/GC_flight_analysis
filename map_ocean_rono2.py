# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 14:31:43 2016

@author: jennyf
"""

import numpy
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
import netCDF4
from mpl_toolkits.basemap import Basemap, addcyclic


def map_ocean_rono2():
    
    '''Quick scripts to map current ocean values of RONO2'''
    

    f = netCDF4.Dataset('model/Chl_2003.geos.4x5.nc','r')
    
    lon = f.variables['lon']
    lat = f.variables['lat']
    chl = f.variables['CHLO_A_S__CHLO_A']
    
    meno3 = numpy.zeros(chl.shape)
    etno3 = numpy.zeros(chl.shape)
    
    for ll in range(len(lat)):
        
        if abs(lat[ll]) < 10.0:
            tmpchl=chl[:,ll,:]
            tmpmeno3=meno3[:,ll,:]
            tmpetno3=etno3[:,ll,:]
            
            tmpmeno3[numpy.where(tmpchl > 0.1)] = 300.0
            tmpmeno3[numpy.where(tmpchl < 0.1)] = 25.0
            tmpmeno3[numpy.where(tmpchl <= 0)] = 0.0
    
            tmpetno3[numpy.where(tmpchl > 0.1)] = 60.0
            tmpetno3[numpy.where(tmpchl < 0.1)] = 5.0
            tmpetno3[numpy.where(tmpchl <= 0)] = 0.0
    
            meno3[:,ll,:]=tmpmeno3
            etno3[:,ll,:]=tmpetno3        
            
        elif abs(lat[ll]) < 20.0:
            tmpchl=chl[:,ll,:]
            tmpmeno3=meno3[:,ll,:]
            tmpetno3=etno3[:,ll,:]
            
            tmpmeno3[:,:] = 5.0
            tmpmeno3[numpy.where(tmpchl <= 0)] = 0.0
            meno3[:,ll,:]=tmpmeno3
    
            tmpetno3[:,:] = 2.0
            tmpetno3[numpy.where(tmpchl <= 0)] = 0.0
            etno3[:,ll,:]=tmpetno3
    
                    
    # cut to one month only
    meno3 = meno3[0,:,:]
    etno3 = etno3[0,:,:]
    
    # mask zeros
    meno3 = numpy.ma.array(meno3,mask=meno3 <= 0)
    etno3 = numpy.ma.array(etno3,mask=etno3 <= 0)
    
    meno3,lon0 = addcyclic(meno3, lon)
    etno3,lon = addcyclic(etno3, lon)
            
    # transform lon/lat to coordinate grid
    lon,lat = numpy.meshgrid(lon,lat)
        
    # Set up the map
    map=Basemap(projection='robin',lon_0=-130,lat_0=0)
    
    f=pyplot.figure()
        
    # plot the data
    col=map.pcolormesh(lon,lat,meno3,latlon=True,vmin=1,vmax=300,
                       cmap=pyplot.get_cmap('PuRd'),norm=colors.LogNorm())
        
    # improve the map
    map.drawcoastlines()
    map.fillcontinents()
    
    # add a color bar
    cb = map.colorbar(col, "bottom")
    cb.set_label('MeNO3, ppt')
    cb.set_ticks([5, 25, 300])
    cb.set_ticklabels([5, 25, 300])
    pyplot.show()
    
    f=pyplot.figure()
        
    # plot the data
    col=map.pcolormesh(lon,lat,etno3,latlon=True,vmin=1,vmax=60,
                       cmap=pyplot.get_cmap('PuRd'),norm=colors.LogNorm())
        
    # improve the map
    map.drawcoastlines()
    map.fillcontinents()
    
    # add a color bar
    cb = map.colorbar(col, "bottom")
    cb.set_label('EtNO3, ppt')
    cb.set_ticks([2, 5, 60])
    cb.set_ticklabels([2, 5, 60])
    pyplot.show()
    
    
def map_ocean_rono2_flux(type='net',
                         mindatam=None,maxdatam=None,
                         mindatae=None,maxdatae=None):
    
    f = netCDF4.Dataset('model/RG/gc.20130101.nc','r')
    
    lon = f.variables['lon']
    lat = f.variables['lat']
    meno3_op = f.variables['ACETSRCE__MeNO3op']
    meno3_ol = f.variables['ACETSRCE__MeNO3ol']

    etno3_op = f.variables['ACETSRCE__EtNO3op']
    etno3_ol = f.variables['ACETSRCE__EtNO3ol']
    
    if type == 'P':
        meno3=meno3_op
        etno3=etno3_op
        title='Flux out of ocean'
    elif type == 'L':
        meno3=meno3_ol
        etno3=etno3_ol
        title='Flux into ocean'
    else:
        meno3=meno3_op-meno3_ol
        etno3=etno3_op-etno3_ol
        title='Net Flux'

    # cut to one month only
    meno3 = meno3[0,:,:]
    etno3 = etno3[0,:,:]
    
    # mask zeros
    meno3 = numpy.ma.array(meno3,mask=meno3 <= 0)
    etno3 = numpy.ma.array(etno3,mask=etno3 <= 0)
 
    meno3,lon0 = addcyclic(meno3, lon)
    etno3,lon = addcyclic(etno3, lon)
    
    # Get min and max data values from GC if needed
    if mindatam is None:
        mindatam = numpy.min(meno3)
    if maxdatam is None:
        maxdatam = numpy.max(meno3)
    if mindatae is None:
        mindatae = numpy.min(etno3)
    if maxdatae is None:
        maxdatae = numpy.max(etno3)

            
    # transform lon/lat to coordinate grid
    lon,lat = numpy.meshgrid(lon,lat)
        
    # Set up the map
    map=Basemap(projection='robin',lon_0=-130,lat_0=0)
    
    f=pyplot.figure()
        
    # plot the data
    col=map.pcolormesh(lon,lat,meno3,latlon=True,
                       vmin=mindatam,vmax=maxdatam,
                       cmap=pyplot.get_cmap('PuRd'))
        
    # improve the map
    map.drawcoastlines()
    map.fillcontinents()
    
    # add a color bar
    cb = map.colorbar(col, "bottom")
    cb.set_label('MeNO3 Flux, molec/cm2/s')
    pyplot.title(title)
  
    f=pyplot.figure()
        
    # plot the data
    col=map.pcolormesh(lon,lat,etno3,latlon=True,
                       vmin=mindatae,vmax=maxdatae,
                       cmap=pyplot.get_cmap('PuRd'))
        
    # improve the map
    map.drawcoastlines()
    map.fillcontinents()
    
    # add a color bar
    cb = map.colorbar(col, "bottom")
    cb.set_label('EtNO3 Flux, molec/cm2/s')
    pyplot.title(title)
        
    
    
map_ocean_rono2_flux(type='P',mindatam=0,maxdatam=2e8)
map_ocean_rono2_flux(type='L',mindatam=0,maxdatam=2e7)
