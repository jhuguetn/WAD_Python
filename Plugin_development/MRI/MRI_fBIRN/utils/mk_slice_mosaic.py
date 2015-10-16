#!/usr/bin/python

# Copyright (c) and all credits to Russ Poldrack
# https://github.com/poldrack/fmriqa

# Adapted 2015-08-26, Jordi Huguet, Dept. Radiology AMC Amsterdam

####################################
__author__      = 'Jordi Huguet'  ##
__dateCreated__ = '20150826'      ##
__version__     = '0.1'           ##
__versionDate__ = '20150828'      ##
####################################

# mk_slice_mosaic.py
# Script for, given an image, make an image mosaic and save as a file
#
# TO DO:
# - ...
# - 

import sys
#import nibabel as nib
import numpy
import matplotlib
# by default matplotlib ships configured to work with a graphical user interface which may require an X11 connection (errors running on background!)
# 	+info:: http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
matplotlib.use('Agg') 
import matplotlib.pyplot as plot

def mk_slice_mosaic(imgdata,outfile,title,contourdata=[],ncols=6,colorbar=True):
    if imgdata.shape[0]==imgdata.shape[1]==imgdata.shape[2]:
        min_dim=2
    else:
        min_dim=numpy.where(numpy.min(imgdata.shape[0:3])==imgdata.shape[0:3])[0]

        slice_dims=numpy.where(numpy.min(imgdata.shape[0:3])!=imgdata.shape[0:3])[0]        
        
    #if not len(imgdata.shape)==3:
    #    sys.stdout.write(__doc__)
    if 1:
        nrows=int(numpy.ceil(numpy.float(imgdata.shape[min_dim])/ncols))
        mosaic=numpy.zeros((nrows*imgdata.shape[slice_dims[0]],ncols*imgdata.shape[slice_dims[1]]))
        if not contourdata==[]:
            contourmosaic=numpy.zeros((nrows*imgdata.shape[slice_dims[0]],ncols*imgdata.shape[slice_dims[1]]))
        ctr=0
        
        for row in range(nrows):
            rowstart=row*imgdata.shape[slice_dims[0]]
            rowend=(row+1)*imgdata.shape[slice_dims[0]]
            for col in range(ncols):
                if ctr<imgdata.shape[min_dim]:
                    colstart=col*imgdata.shape[slice_dims[1]]
                    colend=(col+1)*imgdata.shape[slice_dims[1]]
                    
                    if min_dim==0:
                        imgslice=imgdata[ctr,:,::-1].T
                        if not contourdata==[]:
                            contourslice=contourdata[ctr,:,::-1].T             
                    elif min_dim==1:
                        imgslice=imgdata[:,ctr,::-1].T
                        if not contourdata==[]:
                            contourslice=contourdata[:,ctr,::-1].T       
                    elif min_dim==2:
                        imgslice=imgdata[:,::-1,ctr].T
                        if not contourdata==[]:
                            contourslice=contourdata[:,::-1,ctr].T          
                    try:
                        mosaic[rowstart:rowend,colstart:colend]=imgslice
                    except:
                        mosaic[rowstart:rowend,colstart:colend]=imgslice.T
                        
                    if not contourdata==[]:
                        try:
                            contourmosaic[rowstart:rowend,colstart:colend]=contourslice
                        except:
                            #contourmosaic[rowstart:rowend,colstart:colend]=contourslice.T <--- crashes: could not broadcast input array from shape (96,96,106) into shape (96,96)
                            contourmosaic[rowstart:rowend,colstart:colend]=contourslice.T[:,:,0]
                    ctr+=1
    elif dir=='saggital':
        nrows=int(numpy.ceil(numpy.float(imgdata.shape[1])/ncols))
        mosaic=numpy.zeros((nrows*imgdata.shape[2],ncols*imgdata.shape[1]))
        if not contourdata==[]:
            contourmosaic=numpy.zeros((nrows*imgdata.shape[0],ncols*imgdata.shape[2]))
        ctr=0
        for row in range(nrows):
            rowstart=row*imgdata.shape[2]
            rowend=(row+1)*imgdata.shape[2]
            for col in range(ncols):
                if ctr<imgdata.shape[1]:
                    colstart=col*imgdata.shape[0]
                    colend=(col+1)*imgdata.shape[0]
                    mosaic[rowstart:rowend,colstart:colend]=numpy.flipud(imgdata[ctr,:,:].squeeze().T)
                    if not contourdata==[]:
                        contourmosaic[rowstart:rowend,colstart:colend]=contourdata[ctr,::-1,:].T
                    ctr+=1

    plot.figure(figsize=(8,8))
    plot.imshow(mosaic,cmap=plot.cm.gray)
    if not title=='':
        plot.title(title)
    if colorbar:
        plot.colorbar()
    if not contourdata==[]:
        plot.contour(contourmosaic,colors='red')
    plot.savefig(outfile,bbox_inches='tight')
    plot.close()
