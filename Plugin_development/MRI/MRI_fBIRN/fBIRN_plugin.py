#!/usr/bin/python

# Created 2015-08-20, Jordi Huguet, Dept. Radiology AMC Amsterdam (Z0/3TMRI)

####################################
__author__      = 'Jordi Huguet'  ##
__dateCreated__ = '20150820'      ##
__version__     = '0.1.3'         ##
__versionDate__ = '20151127'      ##
####################################

# - MRI fBIRN Quality Control plugin -
# Initiative for developing a fBIRN-based QC MRI module. Included tests:
#   - SFNR :: signal-to-fluctuation noise ratio

#__version__ = 20151016

# IMPORT FUNCTIONS
import os
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

    if ( len(data.series_filelist[0]) > 1 ) or ( len(data.series_filelist[0]) == 1 and os.path.isdir(data.series_filelist[0][0]) ):
        
        fileList = data.series_filelist[0]
        if len(data.series_filelist[0]) == 1 :
            fileList = data.series_filelist[0][0]
        # read/load a list of DICOM files
        seriesDataList = pydicom_series.read_files(fileList,showProgress=True, readPixelData=True,skipNonImageFiles=True)

        # check number of series in the array/list
        if len(seriesDataList) != 1:
            raise ValueError("{} Such test is supposed to apply solely to a single series/scan. Something went wrong...".format(logTag))

        seriesData = seriesDataList[0]

        nTemporalPositions = int(seriesData.info["0020","0105"].value)

        # Image pixeldata seems to be transposed when read using wadwrapper_lib methods...
        pixeldataIn = seriesData.get_pixel_array()
        pixeldataIn = np.transpose(pixeldataIn)
        new3rdDimension = int(pixeldataIn.shape[2])/nTemporalPositions
        print '[Debug] Image dimensions: ' +  str(pixeldataIn.shape)
        print '[Debug] New dimensions: (%s, %s, %s, %s)' %(str(pixeldataIn.shape[0]),str(pixeldataIn.shape[1]),str(new3rdDimension),str(nTemporalPositions))
        pixeldataIn = np.reshape(pixeldataIn, (pixeldataIn.shape[0],pixeldataIn.shape[1],new3rdDimension,nTemporalPositions))

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

        # Image pixeldata seems to be transposed when read using wadwrapper_lib methods...
        pixeldataIn = np.transpose(pixeldataIn)
        new3rdDimension = int(pixeldataIn.shape[2])/nTemporalPositions
        print '[Debug] Image dimensions: ' + str(pixeldataIn.shape)
        print '[Debug] New dimensions: (%s, %s, %s, %s)' %(str(pixeldataIn.shape[0]),str(pixeldataIn.shape[1]),str(new3rdDimension),str(nTemporalPositions))
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
    #if len(output['spikes']) > 0:
        #print '[info] nspikes,%d'%len(output['spikes'])
    if len(output['spikes']) > 0:
        results.addBool('spikes', True, level=2) #Spikes detected
    else :
        results.addBool('spikes', False, level=2) #No spikes detected

    #------------------------------------------------------------------





