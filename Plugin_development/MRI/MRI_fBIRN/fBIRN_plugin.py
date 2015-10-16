#!/usr/bin/python

# Created 2015-08-20, Jordi Huguet, Dept. Radiology AMC Amsterdam (Z0/3TMRI)

####################################
__author__      = 'Jordi Huguet'  ##
__dateCreated__ = '20150820'      ##
#__version__     = '0.1'           ##
__versionDate__ = '20150820'      ##
####################################

# - MRI fBIRN Quality Control plugin -
# Initiative for developing a fBIRN-based QC MRI module. Included tests:
#   - SFNR :: signal-to-fluctuation noise ratio

__version__ = 20151016

# IMPORT FUNCTIONS
#import os
import numpy as np
#import matplotlib
#if not 'MPLCONFIGDIR' in os.environ:
#    os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import fBIRN_lib

try:
    import wadwrapper_lib
except ImportError:
    from pyWADLib import wadwrapper_lib


# FUNCTIONS
def logTag():
    return "[fBIRN_plugin] "

def runTest(data,results,params):
    """
    MRI_fBIRN Checks: fBIRN QA
      Signal-to-Fluctuation-Noise Ratio (SFNR)
      ... (more tests to come eventually)

    Workflow:
        1. Read image or sequence
        2. Run test
        3. Build xml output
    """

    #(1) Reads input file(s) list as a series (scan) single DICOM object
    #    returns DICOM header, raw pizelData object scaled and type of current DICOM object { 2D, 3D, ... }
    dcmInfile,pixeldataIn,dicomMode = wadwrapper_lib.prepareInput(data.series_filelist[0],headers_only=False,logTag=logTag())

    #print 'dcmInfile: %s' %dcmInfile
    #print 'pixeldataIn shape: %s' %str(pixeldataIn.shape)
    #print 'dicomMode: %s' %dicomMode

    #(2)
    seriesDesc = wadwrapper_lib.readDICOMtag("0008,103E",dcmInfile)   # Check Series Description content
    protocolName = wadwrapper_lib.readDICOMtag("0018,1030",dcmInfile) # Check Protocol Name content

    # Check if is an fBIRN phantom image
    if 'FBIRN' in seriesDesc.upper() and 'FBIRN' in protocolName.upper() :
        # Check if is a multiframe image (MRI)
        if dicomMode not in ["3D","Enhanced"] :
            raise ValueError("{} Input dataset type not suitable --> 2D dataset or no PixelData found".format(logTag))

        nTemporalPositions = wadwrapper_lib.readDICOMtag("0020,0105",dcmInfile)

        if nTemporalPositions is '' :
            raise ValueError("{} Input dataset type not suitable --> no dynamics/temporal positions found".format(logTag))
        else:
            transposed_pixeldataIn = pixeldataIn.T
            new3rdDimension = int(transposed_pixeldataIn.shape[2])/int(nTemporalPositions)
            reshaped_pixeldataIn = np.reshape(transposed_pixeldataIn, (transposed_pixeldataIn.shape[0],transposed_pixeldataIn.shape[1],new3rdDimension,nTemporalPositions))

            output = fBIRN_lib.fBIRN_SFNR(reshaped_pixeldataIn,plot_data=False)

        #(3)
        results.addFloat('mean_SNR', np.mean(output['imgsnr']), quantity='SNR')
        results.addFloat('mean_SFNR', output['meansfnr'], quantity='SFNR')

        #print '[info] SNR, %f'%np.mean(output['imgsnr'])
        #print '[info] SFNR, %f'%output['meansfnr']
        #print '[info] drift, %f'%output['trend'].params[1]
        if len(output['spikes']) > 0:
            print '[info] nspikes,%d'%len(output['spikes'])

    else:
        print '[warning] Not an fBIRN scan --> do nothing'



