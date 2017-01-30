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
from matplotlib import cm

# for testing
import time
#

# My packages
from read_airborne_data import *
from read_gc_data import *

matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['font.sans-serif'] = 'Arial' #- for illustrator...
matplotlib.rcParams['font.sans-serif'] = 'Bitstream Vera Sans'


def etno3_profiles(month,savefig=False,airdata=None,airname='',
               latrange=None,lonrange=None,altrange=None,dalt=1.0,
               mindata=None,maxdata=None,outline=True,
               figname=None,title='',SummaryLegend=False):
    
    """
    This function makes an altitude profile of ETNO3 from model runs with
    a range of ETO2 source reactions switched off,
    averaged globally or over a region (specified with latrange, lonrange.
    
    To specify altitude bounds for profile, use altrange=[X,Y].
    
    To save the figure, use savefig=True.
    To specify figure name set a value for figname.
    
    To specify data axis range, use mindata=X and maxdata=Y keywords."""
    
    # Average over multiple files
    first = True

    # Set up directories with the various runs
    rundir = "/short/m19/jaf574/SEAC4RS/runs/run.C1-C3_RONO2.RG+"
    moddir = ["JPLrates/",     # all reactions
              "testETO2src2/", # no RCO3 reactions
              "testETO2src3/", # no RCHO photolysis
              "testETO2src4/", # no R4O2 reactions
              "testETO2src5/", # no R4N2 photolysis
              "testETO2src6/", # no MEK bphotolysis
              "testETO2src7/", # no LIM+O3 or Br + C2H6
              "testETO2src8/", # no C2H6+NO3
              "testETO2src9/", # no ETP+OH
              "testETO2src10/"] # no C2H6+OH

    labels = ['ALL','RCO3+X','RCHO+hv','R4O2+X','R4N2+hv','MEK+hv',
              'LIM+OH & C2H6+Br','C2H6+NO3','ETP+OH','C2H6+OH','ocean']

    # Set up colors for stacked contribution plot
    colors=cm.Set1(numpy.linspace(0, 1, len(moddir)+1))

    # index for the legend if using summary version only
    LegIdx = [-1,-2,-3,-6,-8,-9,-10]

    # loop over model run directories
    for index,d in enumerate(moddir):

        for m in month:
            filename = rundir+d+'bpch/ctm.bpch.2013'+str(m).zfill(2)+'01'
            GC = read_gc_bpch('ETNO3',filename)

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
        data = data / len(month)    
        alt = alt / len(month)    
            
        # KLUDGE! fix issue where one alt bin undefined
        if numpy.isnan(data[1,8]):
            data[:,8] = numpy.mean([data[:,7],data[:,9]],axis=0)

        # Get the correct contribution for this source
        if index == 0:

            # Get min and max data values from GC if needed
            if mindata is None:
                mindata = numpy.nanmin(data[0,:])*.9
            if maxdata is None:
                maxdata = numpy.nanmax(data[2,:])*1.1
            
            # Get altitude range if needed
            if altrange is None:
                altrange=[alt.min(),alt.max()]
        
            # Set up figure/plot
            f=pyplot.figure()
            ax = pyplot.subplot(111)
    
            # plot the data
            pyplot.plot(data[1,:],alt,'-r',lw=2)
            pyplot.plot(data[2,:],alt,':r',lw=1)
            pyplot.xlim([mindata,maxdata])
            pyplot.ylim(altrange)
            pyplot.minorticks_on()
            pyplot.ylabel('Altitude (km)')
            pyplot.xlabel('ETNO3 (pptv)')
            pyplot.title(title)

        # for all other runs plot stacked contribution
        if index != 0:
           pyplot.fill_betweenx(alt,data[1,:],olddata[1,:],color=colors[index],label=labels[index])

        # Save prior data value
        olddata = data

        # Reset first variable for each directory
        first = True
    
    # Remainder is C2H6+OH (eventually)
    pyplot.fill_betweenx(alt,numpy.zeros(len(alt)),olddata[1,:],color=colors[-1],label=labels[-1])

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if SummaryLegend :
        pyplot.legend([handles[i] for i in LegIdx],[labels[i] for i in LegIdx],loc='best')
    else:
        pyplot.legend(handles[::-1], labels[::-1], loc='best')
    
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
        airlon = airlon[pairdata >= -99999]
        airlat = airlat[pairdata >= -99999]
        airalt = airalt[pairdata >= -99999]
        pairdata = pairdata[pairdata >= -99999]

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
        pyplot.plot(pairmed.statistic,pairalt,color='black',lw=3)
        pyplot.plot(pair25.statistic,pairalt,color='black',ls=':')
        pyplot.plot(pair75.statistic,pairalt,color='black',ls=':')
                             
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
