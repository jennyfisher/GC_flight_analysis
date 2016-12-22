# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:24:59 2016

@author: jennyf
"""

import pandas
import numpy
import matplotlib.pyplot as pyplot

def FastJ_addX(species, QY=1.0,
               source='Higgins',
               T=298,
               plot=True,
               xlim=[289,345]):

    # Units of cross section depend on source
    if (source=='Higgins'):
        scale=1e-20
    else:
        scale=1e0

    # Define wavelength bins used in GEOS-Chem Fast J (v9 only!)
    lambda_bins = numpy.array([289.00, 298.25, 307.45, 312.45,
                               320.30, 345.00, 412.45, 850.00])
    
    # Read solar wavelength data provided with Fast J
    solar_file = 'data/FastJ/solar-p05nm-UCI.dat'
    
    solar = pandas.read_table(solar_file,
                              sep=' ',
                              skiprows=2,
                              header=None,
                              names=['W','F'],
                              skipinitialspace=True)
    
    lambda_solar = solar['W'].values
    flux_solar   = solar['F'].values
    
    # Also read cross section data from Higgins et al. 2014
    xsec_file = 'data/'+source+'_XSEC/'+species+'.csv'
    
    xsec = pandas.read_table(xsec_file,
                             sep=',',
                             header=0,
                             skip_blank_lines=True)
                             
    lambda_xsec = xsec['Wavelength_nm'].values
    xsec_meas   = xsec['xsec_1e-20_cm2_molec_-1'].values*scale
    
    # Some species have T-dependence to deal with...
    if (T != 298 ):
        
        Tdep_file='data/MPI_XSEC/'+species+'_Tdep.csv'
        
        Tdep = pandas.read_table(Tdep_file,
                                 sep=',',
                                 header=0,
                                 skip_blank_lines=True)
                                 
        lambda_Tdep = Tdep['lambda'].values
        Bval        = Tdep['B*1e3'].values*1e-3
        
        Bval = numpy.interp(lambda_xsec,lambda_Tdep,Bval)
        
        # Calculate cross section at given T
        xsec_meas = xsec_meas * numpy.exp(Bval*(T-298.))

    
    # Restrict solar data to range of xsec data
    flux_solar=flux_solar[(lambda_solar >= lambda_xsec.min()) & 
                          (lambda_solar <= lambda_xsec.max())]
    lambda_solar=lambda_solar[(lambda_solar >= lambda_xsec.min()) & 
                              (lambda_solar <= lambda_xsec.max())]
    
    # Interpolate measured cross sections to same resolution as solar flux
    xsec_interp = numpy.interp(lambda_solar,lambda_xsec,xsec_meas)
    
    # Pre-compute the flux * xsec * quantum yield at each wavelength
    prod = xsec_interp * flux_solar * QY
    
    # assign bin # (1-7) to each solar bin
    # note that the bin will be 0 for values < min(bins)
    solar_bin_index = numpy.digitize(lambda_solar,lambda_bins)
    
    # Set up array for final weighted cross sections
    xsec_wm = numpy.zeros(len(lambda_bins)-1)
    
    # loop over bins
    for w in range(len(lambda_bins)-1):
        
        # find relevant points (first bin is 1 not 0)
        num   = prod[(solar_bin_index == w+1)]
        denom = flux_solar[(solar_bin_index == w+1)]
    
        if len(num) > 0:
            xsec_wm[w] = num.sum()/denom.sum()
    
    
    # Plot?
    if (plot):
        weff = numpy.array([294, 303, 310, 316, 333, 380, 574])
        
        pyplot.plot(lambda_solar,xsec_interp,'o',color='orange')
        pyplot.plot(lambda_xsec,xsec_meas,'o',color='red')
        pyplot.plot(weff,xsec_wm,'o',color='blue')
        
        for w in lambda_bins:
            pyplot.plot([w,w],[0,xsec_meas.max()],'--',color='gray')
        
        pyplot.xlim(xlim)
        pyplot.ylim([0,2e-20])
        pyplot.title(species+' Cross Sections from '+source)
        pyplot.show()
        
            
    return(xsec_wm)

def compare_MPI_Higgins():
    xsec_MeNO3 = FastJ_addX('MeNO3',QY=1.0,source='MPI')
    xsec_MeNO3_H = FastJ_addX('MeNO3',QY=1.0,source='Higgins')
    
    xsec_EtNO3 = FastJ_addX('EtNO3',QY=1.0,source='MPI')
    xsec_EtNO3_H = FastJ_addX('EtNO3',QY=1.0,source='Higgins')
    
    xsec_nPrNO3 = FastJ_addX('nPrNO3',QY=1.0,source='MPI')
    xsec_nPrNO3_H = FastJ_addX('nPrNO3',QY=1.0,source='Higgins')
    
    xsec_iPrNO3 = FastJ_addX('iPrNO3',QY=1.0,source='MPI')
    xsec_iPrNO3_H = FastJ_addX('iPrNO3',QY=1.0,source='Higgins')
    
    
    # Print values & differences
    numpy.set_printoptions(precision=3,suppress=False)
    print 'MeNO3 MPI: ', xsec_MeNO3[0:5]
    print 'MeNO3 Hig: ',xsec_MeNO3_H[0:5]
    print 'Difference: ',100*(xsec_MeNO3_H - xsec_MeNO3)/xsec_MeNO3
    print ''
    print 'EtNO3 MPI: ', xsec_EtNO3[0:5]
    print 'EtNO3 Hig: ',xsec_EtNO3_H[0:5]
    print 'Difference: ',100*(xsec_EtNO3_H - xsec_EtNO3)/xsec_EtNO3
    print ''
    print 'n-PrNO3 MPI: ', xsec_nPrNO3[0:5]
    print 'n-PrNO3 Hig: ',xsec_nPrNO3_H[0:5]
    print 'Difference: ',100*(xsec_nPrNO3_H - xsec_nPrNO3)/xsec_nPrNO3
    print ''
    print 'i-PrNO3 MPI: ', xsec_iPrNO3[0:5]
    print 'i-PrNO3 Hig: ',xsec_iPrNO3_H[0:5]
    print 'Difference: ',100*(xsec_iPrNO3_H - xsec_iPrNO3)/xsec_iPrNO3
    
    # mean for PrNO3 across i-, n-
    xsec_PrNO3 = numpy.mean([[xsec_nPrNO3],[xsec_iPrNO3]],axis=0)
    xsec_PrNO3_H = numpy.mean([[xsec_nPrNO3_H],[xsec_iPrNO3_H]],axis=0)
    print ''
    print 'PrNO3 MPI: ', xsec_PrNO3[0:5]
    print 'PrNO3 Hig: ',xsec_PrNO3_H[0:5]
    print 'Difference: ',100*(xsec_PrNO3_H - xsec_PrNO3)/xsec_PrNO3


def test_MPN():    
    xsec_MPN = FastJ_addX('MPN',QY=1.0,source='MPI')
    print 'MPN',xsec_MPN

def calc_298_240():    

    species_list=['MeNO3','EtNO3','iPrNO3','nPrNO3']

    for species in species_list:
    
        if (species != 'nPrNO3'):
            xsec_240 = FastJ_addX(species,QY=1.0,source='MPI',T=240,plot=False)
            print species+' @ 240K: ',xsec_240[0:6]

        xsec_298 = FastJ_addX(species,QY=1.0,source='MPI',T=298,plot=False)
        print species+' @ 298K: ',xsec_298[0:6]
        
        if (species =='iPrNO3'):
            i_save=xsec_298
        elif (species == 'nPrNO3'):
            n_save=xsec_298
        
    # calculate mean for PrNO3 @ 298K
#    print 'PrNO3 @ 298K: ', \
#          numpy.ravel(numpy.mean([[i_save],[n_save]],axis=0))[0:6]
        


# What should I do right now?        
#compare_MPI_Higgins() # evaluate differences between rec vals & Higgins
#test_MPN()            # test using MPN as test species
calc_298_240()        # calculate values at 298K & 240K

