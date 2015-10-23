#!/usr/bin/python

# Copyright (c) and all credits to Russ Poldrack
# https://github.com/poldrack/fmriqa

# Created 2015-08-26, Jordi Huguet, Dept. Radiology AMC Amsterdam

####################################
__author__	  = 'Jordi Huguet'    ##
__dateCreated__ = '20150826'	  ##
__version__	 = '0.1.5'		      ##
__versionDate__ = '20151016'	  ##
####################################

# fBIRN_lib.py
# Script for performing fMRI quality control (fBIRN QA tools)
#
# TO DO:
# - ...

#__version__ = 20151016

# IMPORT FUNCTIONS
#import os
import sys
import ctypes #provides C compatible data types, allows calling functions from DLLs or shared libraries
flags = sys.getdlopenflags()  #?? (ctypes stuff)
sys.setdlopenflags(flags|ctypes.RTLD_GLOBAL) #?? (ctypes stuff)
#import argparse
#import traceback
#import numpy
import nibabel
#from compute_fd import *
from statsmodels.tsa.tsatools import detrend
#import statsmodels.api
#import matplotlib
# by default matplotlib ships configured to work with a graphical user interface which may require an X11 connection (errors running on background!)
# 	+info:: http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
#matplotlib.use('Agg') # should be specified right after matplotlib is imported
#import matplotlib.pyplot as plt
import sklearn.cross_validation
from Plugin_development.MRI.MRI_fBIRN.utils.mk_slice_mosaic import *
#from matplotlib.mlab import psd
#from plot_timeseries import plot_timeseries
from Plugin_development.MRI.MRI_fBIRN.utils import MAD
import skimage.filters # import scikit-image filters module for on-the-fly mask creation
import skimage.morphology # import scikit-image morphology module

sys.setdlopenflags(flags)  #?? (ctypes stuff)

# Assumptions
# . threshold for spike detection
AJKZ_thresh=25
# . sampling frequency (samples per time unit) at computing the power spectral density
TR = 2.5

def error_and_exit(msg):
	print msg
	sys.stdout.write(__doc__)
	sys.exit(2)

# SIGNAL FLUCTUATION NOISE RATIO
def fBIRN_SFNR(imgdata,maskfile=None,verbose=False,plot_data=True):

	#(1) Preliminars

	#flag for storing sfnr as nii image
	save_sfnr=False

	#get number slices
	nslices=imgdata.shape[2]

	#get number dynamics/voxels/volumes
	ndynamics=imgdata.shape[3]

	#(2) Generate mask image

	#compute the mask image
	if maskfile==None :
		# prior to compute the mask, compute the mean of the whole dynamics to get a single 3D blob
		# compute otsu threshold, and get mask and non-mask data
		voxmean=numpy.mean(imgdata,3)
		maskThreshold = skimage.filters.threshold_otsu(voxmean)
		maskdata=numpy.where(voxmean>maskThreshold, 1, 0)
		# remove black wholes
		maskdata=skimage.morphology.closing(maskdata)
		nonmaskdata=numpy.where(voxmean>maskThreshold, 0, 1)
		maskvox=numpy.where(maskdata>0)
		nonmaskvox=numpy.where(maskdata==0)
		#arr = numpy.memmap( os.path.join(qadir,'mm.raw'), dtype='float32', mode='w+', shape=(96, 96, 35))
		#arr[:,:,:] = maskdata
		#maskimg=nibabel.Nifti1Image(arr,img.get_affine())
		#maskimg.to_filename(os.path.join(qadir,'maskimg.nii.gz'))
		if verbose:
			print '[info] nmaskvox: %i' %len(maskvox[0])
			print '[info] nnonmaskvox: %i' %len(nonmaskvox[0])
	else :
		maskimg=nibabel.load(maskfile)
		maskdata=maskimg.get_data()
		nonmaskdata=numpy.where(maskdata==1, 0, 1)
		maskvox=numpy.where(maskdata>0)
		nonmaskvox=numpy.where(maskdata==0)
		if verbose:
			print '[info] nmaskvox: %i' %len(maskvox[0])
			print '[info] nnonmaskvox: %i' %len(nonmaskvox[0])


	# ToDo: EXCLUDED by now
	# load motion parameters and compute FD and identify bad vols for
	# potential scrubbing (ala Power et al.)

	#motpars=numpy.loadtxt(motfile)
	#fd=compute_fd(motpars)
	#numpy.savetxt(os.path.join(qadir,'fd.txt'),fd)


	#(3) Compute std, mean and coeff of variation values

	# compute the voxels mean (from all dynamics, create a mean volume of the pixels of the array of 3D elements)
	voxmean=numpy.mean(imgdata,3)

	# compute the voxels standard deviation (std) from all dynamics, calculated across all equally positioned pixels in the array of 3D elements
	# std :: variation or dispersion of a set (the spread of a distribution measurement)
	voxstd=numpy.std(imgdata,3)

	# compute the voxel coefficient of variation
	# cv ::  standardized measure of dispersion of a probability distribution (extent of variability in relation to the mean of the population)
	voxcv=voxstd/numpy.abs(voxmean)

	# at computing the coef. of variation --> RuntimeWarning: invalid value encountered in divide (NaN)
	# replace in the 3D object the NaN values for zeros
	voxcv[numpy.isnan(voxcv)]=0
	# replace in the 3D object the larger values for 1s
	voxcv[voxcv>1]=1


	# (4) compute timepoint statistics

	# create arrays for storing per-dynamic computed statistical values for mean,median,mad,cv and SNR (after applying mask)
	maskmedian=numpy.zeros(ndynamics)
	maskmean=numpy.zeros(ndynamics)
	maskmad=numpy.zeros(ndynamics)
	maskcv=numpy.zeros(ndynamics)
	imgsnr=numpy.zeros(ndynamics)

	# loop over dynamics to compute statistical values and store them in the arrays for plotting
	for t in range(ndynamics):
		tmp=imgdata[:,:,:,t]

		tmp_brain=tmp[maskvox]
		tmp_nonbrain=tmp[nonmaskvox]

		maskmad[t]=MAD.MAD(tmp_brain)
		maskmedian[t]=numpy.median(tmp_brain)
		maskmean[t]=numpy.mean(tmp_brain)
		maskcv[t]=maskmad[t]/maskmedian[t]
		# define SNR as mean(masked(Image))/std(noise)
		imgsnr[t]=maskmean[t]/numpy.std(tmp_nonbrain)


	# spike detection - (Greve et al./fBIRN)
	# 1. Remove mean and temporal trend from each voxel.
	# 2. Compute temporal Z-score for each voxel.
	# 3. Average the absolute Z-score (AAZ) within a each slice and time point separately.
	# This gives a matrix with number of rows equal to the number of slices (nSlices)
	# and number of columns equal to the number of time points (nFrames).
	# 4. Compute new Z-scores using a jackknife across the slices (JKZ). For a given time point,
	# remove one of the slices, compute the average and standard deviation of the AAZ across
	# the remaining slices. Use these two numbers to compute a Z for the slice left out
	# (this is the JKZ). The final Spike Measure is the absolute value of the JKZ (AJKZ).
	# Repeat for all slices. This gives a new nSlices-by-nFrames matrix (see Figure 8).
	# This procedure tends to remove components that are common across slices and so rejects motion.
	if verbose:
		print '[info] computing spike stats...'

	detrended_zscore=numpy.zeros(imgdata.shape)
	detrended_data=numpy.zeros(imgdata.shape)

	for i in range(len(maskvox[0])):
		tmp=imgdata[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]
		tmp_detrended=detrend(tmp)
		detrended_data[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]=tmp_detrended
		detrended_zscore[maskvox[0][i],maskvox[1][i],maskvox[2][i],:]=(tmp_detrended - numpy.mean(tmp_detrended))/numpy.std(tmp_detrended)

	loo=sklearn.cross_validation.LeaveOneOut(nslices)
	AAZ=numpy.zeros((nslices,ndynamics))

	for s in range(nslices):
		for t in range(ndynamics):
			AAZ[s,t]=numpy.mean(numpy.abs(detrended_zscore[:,:,s,t]))

	JKZ=numpy.zeros((nslices,ndynamics))

	if verbose:
		print '[info] computing outliers...'

	for train,test in loo:
		for tp in range(ndynamics):
			train_mean=numpy.mean(AAZ[train,tp])
			train_std=numpy.std(AAZ[train,tp])
			JKZ[test,tp]=(AAZ[test,tp] - train_mean)/train_std

	AJKZ=numpy.abs(JKZ)
	spikes=[]
	if numpy.max(AJKZ)>AJKZ_thresh:
		if verbose : print '[warning] Possible spike: Max AJKZ = %f'%numpy.max(AJKZ)
		spikes=numpy.where(numpy.max(AJKZ,0)>AJKZ_thresh)[0]
	#if len(spikes)>0:
	#	numpy.savetxt(os.path.join(qadir,'spikes.txt'),spikes)

	voxmean_detrended=numpy.mean(detrended_data,3)
	voxstd_detrended=numpy.std(detrended_data,3)

	# compute per-voxel SFNR!
	voxsfnr=voxmean/voxstd
	# compute mean SFNR!
	meansfnr=numpy.mean(voxsfnr[maskvox])

	# create plots

	datavars={'imgsnr':imgsnr,'meansfnr':meansfnr,'spikes':spikes}#,'badvols':badvols_expanded_index}
	if plot_data:
		print '[info] Plotting data to image files!'
		#trend=plot_timeseries(maskmean,'Mean signal (unfiltered)',os.path.join(qadir,'maskmean.png'),
		#				plottrend=True,ylabel='Mean MR signal')
		#datavars['trend']=trend
		#plot_timeseries(maskmad,'Median absolute deviation (robust SD)',
		#				os.path.join(qadir,'mad.png'),ylabel='MAD')

		# [DEBUG] store SNR for comparing with precomputed mask (BET)
		#plot_timeseries(imgsnr,'SNR',
		#				os.path.join(qadir,'snr.png'),plottrend=True,ylabel='SNR')

		# ToDo: EXCLUDED
		#plot_timeseries(DVARS,'DVARS (root mean squared signal derivative over brain mask)',
		#				os.path.join(qadir,'DVARS.png'),plotline=0.5,ylabel='DVARS')

		# ToDo: EXCLUDED
		#plot_timeseries(fd,'Framewise displacement',os.path.join(qadir,'fd.png'),
		#				badvols_expanded_index,'Timepoints to scrub (%d total)'%len(badvols),
		#				plotline=0.5,ylims=[0,1],ylabel='FD')

		# Compute the power spectral density
		# (i.e. how the power of a signal is distributed over the different frequencies)
		#psd=matplotlib.mlab.psd(maskmean,NFFT=128,noverlap=96,Fs=1/TR)

		#plt.clf()
		#fig=plt.figure(figsize=[10,3])
		#fig.subplots_adjust(bottom=0.15)
		#plt.plot(psd[1][2:],numpy.log(psd[0][2:]))
		#plt.title('Log power spectrum of mean signal across mask')
		#plt.xlabel('frequency (secs)')
		#plt.ylabel('log power')
		#plt.savefig(os.path.join(qadir,'meanpsd.png'),bbox_inches='tight')
		#plt.close()

		#plt.clf()
		#plt.imshow(AJKZ,vmin=0,vmax=AJKZ_thresh)
		#plt.xlabel('timepoints')
		#plt.ylabel('slices')
		#plt.title('Spike measure (absolute jackknife Z)')
		#plt.savefig(os.path.join(qadir,'spike.png'),bbox_inches='tight')
		#plt.close()

		#if img.shape[0]<img.shape[1] and img.shape[0]<img.shape[2]:
		#	orientation='saggital'
		#else:
		#	orientation='axial'

		#mk_slice_mosaic(voxmean,os.path.join(qadir,'voxmean.png'),'Image mean (with mask)',contourdata=maskdata)
		#mk_slice_mosaic(voxcv,os.path.join(qadir,'voxcv.png'),'Image coefficient of variation')
		#mk_slice_mosaic(voxsfnr,os.path.join(qadir,'voxsfnr.png'),'Image SFNR')

		# WATCH OUT: PDF report creator excluded
		#mk_report(infile,qadir,datavars)

	# WATCH OUT: excluded as is not being used!
	# def save_vars(infile,qadir,datavars):

	#ToDo: save values in a csv file!?
	#datafile=os.path.join(qadir,'qadata.csv')
	#f=open(datafile,'w')
	#f.write('SNR,%f\n'%numpy.mean(datavars['imgsnr']))
	#f.write('SFNR,%f\n'%datavars['meansfnr'])
	#f.write('drift,%f\n'%datavars['trend'].params[1])
	#f.write('nspikes,%d\n'%len(datavars['spikes']))
	#f.write('nscrub,%d\n'%len(datavars['badvols']))
	#f.close()

	if verbose:
		print '[info] SNR, %f'%numpy.mean(datavars['imgsnr'])
		print '[info] SFNR, %f'%datavars['meansfnr']
		#print '[info] drift, %f'%datavars['trend'].params[1]
		if len(datavars['spikes']) > 0:
		    print '[info] nspikes,%d'%len(datavars['spikes'])

	#if save_sfnr:
	#	sfnrimg=nibabel.Nifti1Image(voxsfnr,img.get_affine())
	#	sfnrimg.to_filename(os.path.join(qadir,'voxsfnr.nii.gz'))

	return datavars