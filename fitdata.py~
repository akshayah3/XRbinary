# -*- coding: utf-8 -*-
"""
This file contains functions concerned with fitting the
synthetic light curves to observational data.
"""

MAXDATAPOINTS =  1001
import sys
from .parmeter import  dataparams

def ReadData(band):
    """
    This function reads a data lightcurve from the file data.filename[band][]
    """
    a = open( str(dataparams.filename[band]), "r" )
    lines = a.readlines()
    if a == None:
        print("Cannot open file %s", dataparams.filename[band]) 
        sys.exit("")

    npoints = 0
    minfields = 10
    maxfields = 0
    while lines != None:
        for line in lines:
            linesplit = line.split()
            nfields = len(linesplit)
            npoints = npoints + 1
            phase, flux, stdev = linesplit[0], linesplit[1], linesplit[2]
            if( npoints > MAXDATAPOINTS ):
                sys.exit("Too many phase points in observed light curve.")
            if( nfields < minfields ):
                minfields = nfields
            if( nfields > maxfields ):
                maxfields = nfields
            if( (phase < -0.5) or (phase > 1.0) ):
                print("   Phase of data point %3ld out of range.\n", npoints)
                sys.exit("")
                dataparams.phase[band][npoints] = phase
                dataparams.flux[band][npoints] = flux
                if( nfields >= 3 ):
                    if( stdev <= 0.0 ):
                        print("   Stand. dev. of data point %3ld out of range.\n", npoints)
                        sys.exit("")
                    dataparams.standdev[band][npoints] = stdev
                else: 
                    dataparams.standdev[band][npoints] = 1.0
    if(npoints <= 0 ):
            sys.exit("No data points in the file containing the observed light curve.")
    dataparams.npoints[band] = npoints
    if( maxfields != minfields ):
        sys.exit("One or more standard deviations missing in data file.")

    """   
    /*******************
    
    Put the points in the observed light curve in order of
    increasing orbital phase.
   
    ********************/
    """
    if( npoints == 1 ):
        return
    for i in range(1, dataparams.npoints[band]):
        for j in range(i+1, dataparams.npoints[band]):
            if( dataparams.phase[band][j] < dataparams.phase[band][i] ):
                x = dataparams.phase[band][i]
                dataparams.phase[band][i] = dataparams.phase[band][j]
                dataparams.phase[band][j] = x
                x = dataparams.flux[band][i]
                dataparams.flux[band][i] = dataparams.flux[band][j]
                dataparams.flux[band][j] = x
                x = dataparams.standdev[band][i]
                dataparams.standdev[band][i] = dataparams.standdev[band][j]
                dataparams.standdev[band][j] = x

    return
