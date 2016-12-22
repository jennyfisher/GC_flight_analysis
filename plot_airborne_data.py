# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 16:49:15 2016

@author: jennyf
"""

import os
import numpy
import matplotlib.pyplot as pyplot
from read_ict import read_ict

# Name of directory with data
DataDir = 'data/DISCOVERAQ/'
DataDir = 'data/TexasAQS/'

# Set up logical for reading
first=True

# Loop over files to read in data
for file in os.listdir(DataDir):
   
   if first:
       data = read_ict(DataDir+file)                 # read data
       first=False                                   # reset logical
   else:
       tmp_data = read_ict(DataDir+file)             # read data
       data=numpy.ma.concatenate((data,tmp_data))    # concatenate, maintaining missing values

# variable names
VarNames = data.dtype.names

# look for a specific variable - for now NOy
for name in VarNames:
    if 'noy' in name.lower():
        NOyname = name
    if 'alt' in name.lower():
        ALTname = name

       
## plot as a function of pressure
#pyplot.plot(data['ANs_Cohen'],data['FMS_ALT_PRES_nav'],'o',color='red')
#pyplot.title('ANs')
#pyplot.show()

pyplot.plot(data[NOyname],data[ALTname],'o',color='red')
pyplot.xlabel('NOy')
pyplot.ylabel('altitude')
pyplot.show()
