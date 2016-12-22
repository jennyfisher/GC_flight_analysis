# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:24:59 2016

@author: jennyf
"""

import pandas
import numpy

def X_NO2(WW):
    
    WNM = numpy.array([240.964, 243.902, 246.914, 250.000, 253.165, 256.410, 259.740,
                    263.158, 266.667, 270.270, 273.973, 277.778, 281.690, 285.714,
                    289.855, 294.118, 298.507, 303.030, 307.692, 312.5, 317.5, 322.5,
                    327.5, 332.5, 337.5, 342.5, 347.5, 352.5, 357.5, 362.5, 367.5,
                    372.5, 377.5, 382.5, 387.5, 392.5, 397.5, 402.5, 407.5, 412.5,
                    417.5, 422.5, 427.5])
                    
    X294K = numpy.array([5.77,  2.79,  1.62,  0.998, 1.05,  1.28, 1.58, 2.05, 2.64, 3.24,
                         4.07, 5.21, 6.23, 7.59, 9.51, 11.5, 13.2, 16.1, 18.8, 21.6, 25.3,
                         28.7, 31.7, 35.8, 40.2, 41.8, 46.2, 49.7, 50.9, 54.9, 56.1, 59.0,
                         59.3, 60.1, 63.0, 59.7, 64.4, 58.2, 62.4, 59.1, 59.9, 57.0])
                         
    loc = (WNM[0:len(WNM)-1]<WW) & (WNM[1:len(WNM)]>WW)
    
    XNO2 = X294K[loc]*1e-20
    
    # avoid zero value outside allowable range!
    if XNO2.size == 0:
        XNO2=0.
    
    return(XNO2)


fastj_dir = '/Users/jennyf/Downloads/UCI_fastJ_addX_73C/'
wavel_file= fastj_dir + 'wavel-bins.dat'

# First read the pratmos wavelength bins
wbin=pandas.read_table(wavel_file,
                    sep=' ',
                    skiprows=1,
                    header=None,
                    nrows=78,
                    usecols=[1],
                    skipinitialspace=True)
wbin=wbin.values.ravel()


# Next read the FAST-JX bins
widths = [3]
for n in range(16):
    widths.append(5)
srb=pandas.read_fwf(wavel_file,
                      sep=' ',
                      skiprows=81,
                      header=None,
                      nrows=8,
                      usecols=range(1,16),
                      widths=widths,
                      na_values='.')
srb=numpy.nan_to_num(srb.values)

# Next read the bin assignments
ijx=pandas.read_table(wavel_file,
                    sep=' ',
                    skiprows=90,
                    header=None,
                    nrows=78-15,
                    usecols=[1],
                    skipinitialspace=True)

ijx = ijx.values.ravel()
ijxn = numpy.zeros(15)
ijx=numpy.hstack((ijxn,ijx)).astype(int)

# New file! read solar
solar_file = fastj_dir + 'solar-p05nm-UCI.dat'

solar = pandas.read_table(solar_file,
                          sep=' ',
                          skiprows=2,
                          header=None,
                          names=['W','F'],
                          skipinitialspace=True)

W = solar['W'].values
F = solar['F'].values

# assign bin # (1-77) to each solar bin
# note that the bin will be 0 for values < min(wbin)
ibinj = numpy.digitize(W,wbin)

Fbin = numpy.zeros(len(wbin)-1)
Abin = numpy.zeros(len(wbin)-1)

# loop over all solar values to get weighted output
for i in range(len(F)):
    if ibinj[i] > 0:
        XNO2 = X_NO2(W[i])
        if XNO2 > 0:
            ind = ibinj[i]
            Fbin[ind-1] = Fbin[ind-1] + F[i]
            Abin[ind-1] = Abin[ind-1] + F[i]*XNO2
        
# now loop over bins to get weighted mean
for i in range(len(Fbin)):
    if Fbin[i] > 0:
        Abin[i] = Abin[i]/Fbin[i]
        
FFbin = numpy.zeros(18)
AAbin = numpy.zeros(18)

# loop over bin assignments
for i in range(15,len(ijx)):
    ind = ijx[i]
    print i, ind, Fbin[i]
    FFbin[ind-1] = FFbin[ind-1] + Fbin[i]
    AAbin[ind-1] = AAbin[ind-1] + Fbin[i]*Abin[i]
    
# first 15 w/ SRB??
for i in range(srb.shape[0]):
    for j in range(srb.shape[1]):
        FFbin[j] = FFbin[j]+Fbin[i]*srb[i,j]
        AAbin[j] = AAbin[j]+Fbin[i]*Abin[i]*srb[i,j]
    
# now loop again to get weighted mean
for i in range(len(FFbin)):
    if FFbin[i] > 0:
        AAbin[i] = AAbin[i]/FFbin[i]
        
print AAbin