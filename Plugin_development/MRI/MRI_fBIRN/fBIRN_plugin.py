#!/usr/bin/python

# Created 2015-08-20, Jordi Huguet, Dept. Radiology AMC Amsterdam (Z0/3TMRI)

####################################
__author__      = 'Jordi Huguet'  ##
__dateCreated__ = '20150820'      ##
__version__     = '0.1.2'         ##
__versionDate__ = '20151023'      ##
####################################

# - MRI fBIRN Quality Control plugin -
# Initiative for developing a fBIRN-based QC MRI module. Included tests:
#   - SFNR :: signal-to-fluctuation noise ratio

#__version__ = 20151016

# IMPORT FUNCTIONS
#import os
import struct
import numpy as np
#import matplotlib
#if not 'MPLCONFIGDIR' in os.environ:
#    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import fBIRN_lib
from pyWADLib import pydicom_series

try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib


# FUNCTIONS
def logTag():
    return "[fBIRN_plugin] "

def sfnrTest(data,results,params):
    """
    MRI_fBIRN Checks: fBIRN QA
      Signal-to-Fluctuation-Noise Ratio (SFNR)
      ... (more tests to come eventually)

    Workflow:
        1. Read image or sequence
        2. Run test
        3. Build output
    """

    #(1) Reads input file(s) list as a series (scan) single DICOM object
    #    returns DICOM header, raw pizelData object scaled and type of current DICOM object { 2D, 3D, ... }
    #dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data.series_filelist[0],headers_only=False,logTag=logTag())

    #------------------------------------------------------------------
    pixeldataIn = None

    if len(data.series_filelist[0]) > 1:
        # read/load a list of DICOM files
        seriesDataList = pydicom_series.read_files(data.series_filelist[0],showProgress=True, readPixelData=True,skipNonImageFiles=True)

        # check number of series in the array/list
        if len(seriesDataList) != 1:
            raise ValueError("{} ToDo --> lets see how to handle this, actually...shall we?".format(logTag))

        seriesData = seriesDataList[0]

        nTemporalPositions = int(seriesData.info["0020","0105"].value)
        #options['seriesDesc'] = seriesData.info["0008","103E"].value
        #options['protocolName'] = seriesData.info["0018","1030"].value

        # Image pixeldata seems to be transposed when read using wadwrapper_lib methods...
        pixeldataIn = seriesData.get_pixel_array()
        pixeldataIn = np.transpose(pixeldataIn)

        #print len(data)
        #for seriesData in data:
        #    pixeldataIn = seriesData.get_pixel_array()
        #    pixeldataIn = np.transpose(pixeldataIn)
        #    print pixeldataIn.shape
            #print  seriesData.shape
            #print  seriesData.description
            #print  seriesData.info

        #    print seriesData.info["0028","0030"].tag
        #    print seriesData.info["0028","0030"].VR


    elif ( len(data.series_filelist[0]) == 1) and ( wadwrapper_lib.testIfEnhancedDICOM(data.series_filelist[0][0]) ):
         # read/load a single DICOM file
        dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareEnhancedInput(data.series_filelist[0][0],headers_only=False)

        # DICOM keeps NumberOfTemporalPositions nested in sequence items/subitems. And DCM4CHEE seems not to keep them at all (!?).
        # Workaround: Use Philips private tag, should be read as a bit string and converted to string/integer/whatever
        if 'PHILIPS' not in (wadwrapper_lib.readDICOMtag("0008,0070",dcmInfile)).upper():
            raise ValueError("{} Input enhanced dataset type not suitable --> no dynamics/temporal information found".format(logTag))

        nSlicesRaw = wadwrapper_lib.readDICOMtag("2001,1018",dcmInfile)
        try:
            nSlices = struct.unpack("<L",nSlicesRaw)[0]
        except:
            nSlices = nSlicesRaw
        nTemporalPositions = int(pixeldataIn.shape[0])/int(nSlices)

        #options['seriesDesc'] = wadwrapper_lib.readDICOMtag("0008,103E",dcmInfile)   # Check Series Description content
        #options['protocolName'] = wadwrapper_lib.readDICOMtag("0018,1030",dcmInfile) # Check Protocol Name content

        # Image pixeldata seems to be transposed when read using wadwrapper_lib methods...
        pixeldataIn = np.transpose(pixeldataIn)
        new3rdDimension = int(pixeldataIn.shape[2])/nTemporalPositions
        pixeldataIn = np.reshape(pixeldataIn, (pixeldataIn.shape[0],pixeldataIn.shape[1],new3rdDimension,nTemporalPositions))

    else:
        raise ValueError("{} Input dataset type cannot be determined or is not compatible".format(logTag))

    #(2)

    #if 'BIRN' not in options['seriesDesc'] or 'BIRN' not in options['protocolName']:
    #   raise ValueError("{} Input dataset type not suitable".format(logTag))
    # OR
    #   print '[Warning] Not an fBIRN scan --> do nothing'

    #Check if pixeldataIn is actually a numpy array
    if type(pixeldataIn).__module__ != np.__name__ :
        raise ValueError("{} Unable to pull out the pixel data of the incomming DICOM file(s)".format(logTag))

    output = fBIRN_lib.fBIRN_SFNR(pixeldataIn,plot_data=False)

    #(3)

    results.addFloat('mean_SNR', np.mean(output['imgsnr']), quantity='SNR', level=2) #quantity is actually magnitude in the WAD-QC app
    results.addFloat('mean_SFNR', output['meansfnr'], quantity='SFNR', level=2) #quantity is actually magnitude in the WAD-QC app

    #print '[info] SNR, %f'%np.mean(output['imgsnr'])
    #print '[info] SFNR, %f'%output['meansfnr']
    #print '[info] drift, %f'%output['trend'].params[1]
    if len(output['spikes']) > 0:
        print '[info] nspikes,%d'%len(output['spikes'])

    #------------------------------------------------------------------





