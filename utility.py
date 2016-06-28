# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:53:23 2016

"""

"""
Utility functions needed for the project
"""

import math

def BBSquareIntensity(T, lowlambda, highlambda):
    """
    This function returns the mean specific intensity emitted by a black body 
    with temperature T in a square bandpass from lowlambda to highlambda.  
    Specifically, the function returns ( \int_lowlambda^highlambda B(lambda) d lambda ) / deltalambda
    where B(lambda) is in specific intensity form.  It uses the precalculated 
    table Zzeta.dat and does a logarithmic interpolation in the table. The 
    wavelengths must be in centimeters but the intensity is in ergs/cm^2/sec/Angstrom.
    """
    c = 2.99792e+10
    h = 6.626075e-27
    k = 1.38065e-16
    
    if( lowlambda >= highlambda ):
        print("lowlambda must be lt highlambda in BBSquareIntensity().")
    if( (lowlambda == 0) or (highlambda == 0) ):
        print("Both wavelengths must be gt zero in BBSquareIntensity().")

    Zeta2 = ( h * c ) / (  lowlambda * k * T )
    Zeta1 = ( h * c ) / ( highlambda * k * T )
    
    if( (Zeta1 < deltaBBzeta) or (Zeta2 < deltaBBzeta) ):
        print("lambda x T too high in BBSquareIntensity().")

    if( Zeta1 >= BBzetamax ):
        meanI = 0.0
        return meanI
    else:
        logZeta1 = math.log( Zeta1 )
        nlow = Zeta1 / deltaBBzeta
        ZetaLow = nlow * deltaBBzeta
        logZetaLow = math.log( ZetaLow )
        nhigh = nlow + 1
        logZetaHigh = math.log( ZetaLow + deltaBBzeta )
        logZBBzetaLow  = math.log( ZBBzeta[nlow]  )
        logZBBzetaHigh = math.log( ZBBzeta[nhigh] )
        slope = (logZBBzetaHigh - logZBBzetaLow) / (logZetaHigh - logZetaLow)
        logZBBzeta1 = logZBBzetaLow + (logZeta1 - logZetaLow) * slope
        ZBBzeta1 = math.exp( logZBBzeta1 )

    if( Zeta2 >= BBzetamax ):
      ZBBzeta2 = ZBBzeta[maxBBzetaindex]

    else:
        logZeta2 = math.log( Zeta2 )
        nlow = Zeta2 / deltaBBzeta
        ZetaLow = nlow * deltaBBzeta
        logZetaLow = math.log( ZetaLow )
        nhigh = nlow + 1
        logZetaHigh = math.log( ZetaLow + deltaBBzeta )
        logZBBzetaLow  = math.log( ZBBzeta[nlow]  )
        logZBBzetaHigh = math.log( ZBBzeta[nhigh] )
        slope = (logZBBzetaHigh - logZBBzetaLow) / (logZetaHigh - logZetaLow)
        logZBBzeta2 = logZBBzetaLow + (logZeta2 - logZetaLow) * slope
        ZBBzeta2 = math.exp( logZBBzeta2 )

    meanI = ( T*T*T*T ) * (ZBBzeta2 - ZBBzeta1) / (highlambda - lowlambda)
    meanI *= 1.0e-08

    return meanI

    def BBFilterIntensity(T, filter):
        """
        This function returns the mean specific intensity emitted by a black 
        body with temperature T observed through a filter.
        """
        findex = -1;
        for i in range(maxIBBfilterindex):
            if IBBfilterName[i] == filter:
                findex = i
                break

        if( findex == -1 ):
            print("Unknown filter name in BBFilterIntensity.");

        if( (T < IBBT[0]) or (T > IBBT[maxIBBTindex]) ):
            print("T out of range in BBFilterIntensity.")
        if( T == IBBT[maxIBBTindex] ):
            intensity = IBBtable[maxIBBTindex][findex]
            return( intensity )

        minTindex = (T - IBBTmin) / IBBdeltaT
        maxTindex = minTindex + 1
        slope = (IBBtable[maxTindex][findex] - IBBtable[minTindex][findex]) / IBBdeltaT
        intensity = slope * (T - IBBT[minTindex]) + IBBtable[minTindex][findex]

        return intensity
