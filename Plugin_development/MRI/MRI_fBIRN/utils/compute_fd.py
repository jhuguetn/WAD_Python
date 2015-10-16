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

# compute_fd.py
# Script for implementing the frame displacement measure described by Power et al.(2011):
# 	- http://www.ncbi.nlm.nih.gov/pubmed/22019881
# 	- http://www.sciencedirect.com/science/article/pii/S1053811911011815
#
# TO DO:
# - ...
# - 

import numpy
import os
import sys

# load the datafile
def compute_fd(motpars):

    # compute absolute displacement
    dmotpars=numpy.zeros(motpars.shape)
    
    dmotpars[1:,:]=numpy.abs(motpars[1:,:] - motpars[:-1,:])
    
    # convert rotation to displacement on a 50 mm sphere
    # mcflirt returns rotation in radians
    # from Jonathan Power:
    # The conversion is simple - you just want the length of an arc that a rotational
    # displacement causes at some radius. Circumference is pi*diameter, and we used a 5
    # 0 mm radius. Multiply that circumference by (degrees/360) or (radians/2*pi) to get the 
    # length of the arc produced by a rotation.
    
    
    headradius=50
    disp=dmotpars.copy()
    disp[:,0:3]=numpy.pi*headradius*2*(disp[:,0:3]/(2*numpy.pi))
    
    FD=numpy.sum(disp,1)
    
    return FD