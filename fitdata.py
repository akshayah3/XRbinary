# -*- coding: utf-8 -*-
"""
This file contains functions concerned with fitting the
synthetic light curves to observational data.
"""

def ReadData(band):
    """
    This function reads a data lightcurve from the file data.filename[band][]
    """
    a = open( str(data.filename[band]), "r" )
    if a == Null:
        print("Cannot open file %s", data.filename[band]) 
        Quit("")

    npoints = 0
    minfields = 10
    maxfields = 0
    while( fgets( inputline, 80, a) != NULL ):
        npoints = npoints + 1
        if( npoints > MAXDATAPOINTS ):
            Quit("Too many phase points in observed light curve.")
        nfields = sscanf( inputline, "%lf  %lf  %lf", &phase, &flux, &stdev )
        if( nfields < minfields ):
            minfields = nfields
        if( nfields > maxfields ):
            maxfields = nfields
        if( (phase < -0.5) or (phase > 1.0) ):
            print("   Phase of data point %3ld out of range.\n", npoints)
            Quit("")
        data.phase[band][npoints] = phase
        data.flux[band][npoints] = flux
        if( nfields >= 3 ):
            if( stdev <= 0.0 ):
                print("   Stand. dev. of data point %3ld out of range.\n", npoints)
                Quit("")
            data.standdev[band][npoints] = stdev
            else: 
                data.standdev[band][npoints] = 1.0
    if( npoints <= 0 ):
        Quit("No data points in the file containing the observed light curve.")
    data.npoints[band] = npoints
    if( maxfields != minfields ):
        Quit("One or more standard deviations missing in data file.")

    """   
    /*******************
    
    Put the points in the observed light curve in order of
    increasing orbital phase.
   
    ********************/
    """
    if( npoints == 1 ):
        return
    for i in range(1, data.npoints[band]):
        for j in range(i+1, data.npoints[band]):
            if( data.phase[band][j] < data.phase[band][i] ):
                x = data.phase[band][i]
                data.phase[band][i] = data.phase[band][j]
                data.phase[band][j] = x
                x = data.flux[band][i]
                data.flux[band][i] = data.flux[band][j]
                data.flux[band][j] = x
                x = data.standdev[band][i]
                data.standdev[band][i] = data.standdev[band][j]
                data.standdev[band][j] = x

   return
