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

# plot_timeseries.py
#
# TO DO:
# - ...
# - 

import numpy
import matplotlib
# by default matplotlib ships configured to work with a graphical user interface which may require an X11 connection (errors running on background!)
# 	+info:: http://matplotlib.org/faq/howto_faq.html#matplotlib-in-a-web-application-server
matplotlib.use('Agg') 
import matplotlib.pyplot as plot
import statsmodels.api


def plot_timeseries(data,title,outfile,markers=[],markername=[],
                    plottrend=False,ylims=[],plotline=[],xlabel='timepoints',
                    ylabel=[]):

    fig=plot.figure(figsize=[10,3])
    fig.subplots_adjust(bottom=0.15)
    plot.plot(data)
    datarange=numpy.abs(numpy.max(data)-numpy.min(data))
    ntp=len(data)

    if ylims==[]:
        axislims=[0,ntp+1,numpy.min(data) - datarange*0.1,numpy.max(data) + datarange*0.1]
    else:
        axislims=[0,ntp-1,ylims[0],ylims[1]]

    if plottrend:
        X=numpy.vstack((numpy.ones(ntp), numpy.arange(ntp)- numpy.mean(numpy.arange(ntp)),numpy.arange(ntp)**2-numpy.mean(numpy.arange(ntp)**2))).T
        model=statsmodels.api.OLS(data,X)
        results=model.fit()

    plot.axis(axislims)
    fig.hold('on')
    if len(markers)>0:
        for s in markers:
            lp,=plot.plot([s,s],axislims[2:4],linewidth=2,color='red')
        plot.legend([lp],[markername])
    if plottrend:
        plot.plot(numpy.dot(X,results.params),color='black')
    if plotline:
        plot.plot([0,ntp],[plotline,plotline])
        
    plot.title(title)
    plot.xlabel(xlabel)
    if ylabel:
        plot.ylabel(ylabel)
    plot.savefig(outfile,bbox_inches='tight')
    plot.close()

    if plottrend:
        return results
    else:
        return []
