#!/usr/bin/python

# Copyright (c) and all credits to Russ Poldrack
# https://github.com/poldrack/fmriqa

# Adapted 2015-08-26, Jordi Huguet, Dept. Radiology AMC Amsterdam

####################################
__author__      = 'Jordi Huguet'  ##
__dateCreated__ = '20150826'      ##
__version__     = '0.1'           ##
__versionDate__ = '20150826'      ##
####################################

# MAD.py
# Script for computing the Median Absolute Deviation (depicts statistical dispersion) along given axis of an array:
#    
# 	median(abs(a - median(a))) / c (i.e. compute the median of the deviation/distance of a to its meadian )
# 	c = 0.6745 is the constant to convert from MAD to std; it is used by default (??)
#
# TO DO:
# - ...
# - 

import numpy

def MAD(a, c=0.6745, axis=0):
   
    good = (a==a)
    a = numpy.asarray(a, numpy.float64)
    if a.ndim == 1:
        d = numpy.median(a[good])
        m = numpy.median(numpy.fabs(a[good] - d) / c)
    else:
        d = numpy.median(a[good], axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(a[good],0,axis)
        else:
            aswp = a[good]
        m = numpy.median(numpy.fabs(aswp - d) / c, axis=0)

    return m
	

