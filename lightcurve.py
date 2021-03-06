# -*- coding: utf-8 -*-
"""

"""

GRIDXTILES = 401
GRIDZTILES = 401
SIGMA =  5.6704e-5
MAX2TILES  =  40506
MAXDISKTILES = 40506

import math
import sys
from utility import CartDotProd, BBSquareIntensity, BBFilterIntensity
from diagnose import InspectYlimits, InspectEscape, InspectHeat2ByADC, InspectHeatDiskBy1, InspectHeatDiskByID, InspectHeatDiskByADC, InspectHeat2By1, InspectHeat2ByID, InspectHeat2ByDisk
from diskgeom import DiskTopH, DiskBottomH, RhoToA
from diskflux import InnerDiskFlambda, InnerDiskTotFlux, ADCTotFlux
from star1 import Star1
from star2 import Star2
from parmeter import flowcontrol, orbitparams, systemparams, wholediskpars
from parmeter import adcpars, thirdlightparams, XYGrid, dataparams, globalvar, CartVector

def MakeLightCurves():
    """
    This function makes the orbital light curve.
    """
    TotalFlux = [0 for i in range(22)]
    if( flowcontrol.thirdlight == "ON" ):
        if( globalvar.verbose == "ON" ):
            print(" Calculating third light fluxes.\n")
        ThirdLight()
    if( globalvar.verbose == "ON" ): 
        print(" Begin calculating the light curves.\n")

    for iphase in range(0, orbitparams.maxpindex):
        globalvar.LCphase[iphase] = orbitparams.phasemin + iphase * orbitparams.deltaphase
        k = iphase / 10
        k = iphase - 10 * k
        if( k == 0 ):
            if( globalvar.verbose == "ON" ):
                print("    phase number {}   phase =  {}\n").format( iphase, 
                                                         globalvar.LCphase[iphase])
        FluxesAtPhase( globalvar.LCphase[iphase], TotalFlux )
        for band in range(1, orbitparams.nbands):
            globalvar.LCflux[band][iphase] = TotalFlux[band]

    if( orbitparams.normalize == "OFF" ):
        if( globalvar.verbose == "ON" ):
	       print(" Normalizing the light curves.\n")
        Normalize()
    else:
        for band in range(1, orbitparams.nbands):
            for iphase in range(0, orbitparams.maxpindex):
	           globalvar.NormLC[band][iphase] = 0.0

    return

def FluxesAtPhase( phase, TotalFlux ):
    """
    This function calculates the total fluxes from the system at
    a particular orbital phase.  It is the heart of the program.
    """
    calcphase = phase - orbitparams.phaseoffset
    if( calcphase < -0.5 ): 
        calcphase += 1.0
    if( calcphase > 1.0 ):
        calcphase -= 1.0
    sunvector = CartVector()
    sunvector.x = -1.0 * math.sin( systemparams.i ) * math.sin( calcphase * 2*math.pi )
    sunvector.y =        math.cos( systemparams.i )
    sunvector.z = -1.0 * math.sin( systemparams.i ) * math.cos( calcphase * 2*math.pi )

    if( flowcontrol.star1 == "ON" ):
        start = CartVector()        
        start.x = 0.0
        start.y = 0.0
        start.z = systemparams.a
        Star1Escape = EscapeFraction( start, sunvector)

    if( flowcontrol.star2 == "ON"):
        T2mu = [0 for i in range(0, MAX2TILES)]
        T2Escape = [0 for i in range(0, MAX2TILES)]
        
        for itile in range(1, Star2.Ntiles):
            T2mu[itile] = (CartDotProd( sunvector, globalvar.T2normCart[itile] ))
            if( T2mu[itile] > 0.0 ):
                start = CartVector()
                start.x = globalvar.T2x[itile]
                start.y = globalvar.T2y[itile]
                start.z = globalvar.T2z[itile]
                T2Escape[itile] = EscapeFraction( start, sunvector)
            else:
                T2Escape[itile] = 0.0

    if( flowcontrol.disk == "ON"):
        for itile in range(1, wholediskpars.Ntiles):
            TDiskmu = [0 for i in range(0, MAXDISKTILES)]
            TDiskmu[itile] = (CartDotProd( sunvector, globalvar.TDisknormCart[itile]))
            if( TDiskmu[itile] > 0.0 ):
                start = CartVector()
                start.x = globalvar.TDiskx[itile]
                start.y = globalvar.TDisky[itile]
                start.z = globalvar.TDiskz[itile]
                TDiskEscape = []
                TDiskEscape.append(EscapeFraction( start, sunvector))
            else:
                TDiskEscape.append(0.0)
        if( flowcontrol.innerdisk == "ON" ):
            start.x = 0.0
            start.y = 0.0
            start.z = systemparams.a
            InnerDiskEscape = EscapeFraction( start, sunvector)


    for band in range(1, orbitparams.nbands):
        TotalFlux[band] = 0.0

        if( flowcontrol.star1 == "ON" ):
            Star1Emitted = Star1.Star1Flambda( orbitparams.filter[band], 
                                   orbitparams.minlambda[band],
                                   orbitparams.maxlambda[band] )
            star1flux = Star1Emitted * Star1Escape
            TotalFlux[band] += star1flux

        if( flowcontrol.star2 == "ON" ):
            star2flux = 0.0
            T2Emitted = [0 for i in range(0, MAX2TILES)]
            if( orbitparams.filter[band] == "SQUARE"):
                for itile in range(1, Star2.Ntiles):
                    if( T2mu[itile] < 0.0 ):
                       T2Emitted[itile] = 0.0
                    else:
   	                  T2Emitted[itile] =  globalvar.T2I[band][itile]*T2mu[itile] * globalvar.T2dS[itile]
                    star2flux += T2Emitted[itile] * T2Escape[itile]
            else:
                for itile in range(1, Star2.Ntiles):
                    if( T2mu[itile] < 0.0 ):
		             T2Emitted[itile] = 0.0
                    else:
                        if( (globalvar.T2T[itile] > globalvar.IperpT[globalvar.maxIperpTindex]) or (globalvar.T2T[itile] < globalvar.IperpT[0]) ):
                            T2Emitted[itile] = globalvar.T2I[band][itile]*T2mu[itile] * globalvar.T2dS[itile]
                        else:
                            T2Emitted[itile] = globalvar.T2I[band][itile]* Star2.ClaretHmu( globalvar.T2T[itile],globalvar.T2logg[itile],orbitparams.filter[band],T2mu[itile] )* T2mu[itile] * globalvar.T2dS[itile]
                    star2flux += T2Emitted[itile] * T2Escape[itile]
            TotalFlux[band] += star2flux
        if( flowcontrol.disk == "ON" ):
            diskflux = 0.0
            TDiskEmitted = [0 for i in range(0, MAX2TILES)]

            for itile in range(1, wholediskpars.Ntiles):
                if( TDiskmu[itile] < 0.0 ):
	              TDiskEmitted[itile] = 0.0
                else:
   	              TDiskEmitted[itile] = globalvar.TDiskI[band][itile] * TDiskmu[itile] * globalvar.TDiskdS[itile]
                diskflux += TDiskEmitted[itile] * TDiskEscape[itile]
            TotalFlux[band] += diskflux
            if( flowcontrol.innerdisk == "ON" ):
                InnerDiskMu = sunvector.y
                InnerDiskEmitted = InnerDiskFlambda( orbitparams.filter[band], orbitparams.minlambda[band],orbitparams.maxlambda[band] )
                InnerDiskflux = InnerDiskMu * InnerDiskEmitted * InnerDiskEscape;
                TotalFlux[band] += InnerDiskflux

        if( flowcontrol.thirdlight == "ON" ):
            TotalFlux[band] += thirdlightparams.addFlux[band]

        if( flowcontrol.diagnostics == "INSPECTESCAPE"):
            if( flowcontrol.diagnoseband == orbitparams.filter[band] ):
                InspectEscape( sunvector, orbitparams.filter[band],
                           Star1Emitted, Star1Escape,
                           T2mu, T2Emitted, T2Escape,
                           TDiskmu, TDiskEmitted, TDiskEscape);
            sys.exit("Quit after INSPECTESCAPE.")

    return

def ThirdLight():
    """
    This function calculates the amount of third light flux to
    to be added to each bandpass.  The flux is actually added in
    function FluxesAtPhase().  
    """
    TotalFlux = [0 for i in range(22)]
    FluxesAtPhase( thirdlightparams.orbphase, TotalFlux )

    for LCband in range(1, orbitparams.nbands):
        found = 0
        for band in range(1, thirdlightparams.nbands):
            if( orbitparams.filter[LCband] == thirdlightparams.filter[band]):
                if( orbitparams.filter[LCband] == "SQUARE"):
                    if( (orbitparams.minlambda[LCband] == thirdlightparams.minlambda[band]) and (orbitparams.maxlambda[LCband] == thirdlightparams.maxlambda[band]) ):
                        found = 1
                        alpha = thirdlightparams.fraction[band]
                else:
                    found = 1
                    alpha = thirdlightparams.fraction[band]
        if( found == 1 ):
            thirdlightparams.addFlux[LCband] = ( alpha / (1.0 - alpha) ) * TotalFlux[LCband]
        else:
            sys.exit("Thirdlight fraction not specified for a bandpass.")

    return

def Normalize():
    """
 
    This function normalizes the orbital light.  Meaning of the
    normalization types:
        MAXVALUE  normvalue
             The maximum value of each individual light curve is 
             normalized to normvalue.
        FITDATA  normfilter  minlambda maxlambda
             The normalization factor for filter "normfilter" is set
             equal to the value that minimizes the (weighted) mean squared
             squared difference between the observed light curve
             and the calculated light curve.  This same normalization
             factor is then applied to all the calculated light curves.
    """
    TotalFlux = [0 for i in range(22)]
    if( orbitparams.normalize == "MAXVALUE"):
        for band in range(1, orbitparams.nbands):
            maxflux = 0.0
            for iphase in range(0, orbitparams.maxpindex):
                if( globalvar.LCflux[band][iphase] > maxflux ):
                    maxflux = globalvar.LCflux[band][iphase]
            if( maxflux == 0.0 ):
                maxflux = 1.0
            normfactor = orbitparams.normvalue / maxflux
            for iphase in range(0, orbitparams.maxpindex):
                globalvar.NormLC[band][iphase] = normfactor * globalvar.LCflux[band][iphase]
    if( orbitparams.normalize == "FITDATA"):
        calcband = 0
        for band in range(1, orbitparams.nbands):
            if( orbitparams.normfilter == orbitparams.filter[band]):
                if( orbitparams.normfilter == "SQUARE"):
                    if( (orbitparams.normMinlambda == orbitparams.minlambda[band]) and (orbitparams.normMaxlambda == orbitparams.maxlambda[band]) ):
       	             calcband = band
                else:
                    calcband = band
        if( calcband == 0 ):
            sys.exit("Could not find a matching calculated band in Normalize().")
                       
        databand = 0
        for band in range(1, dataparams.nbands):
            if( orbitparams.normfilter == dataparams.filter[band]):
                if( orbitparams.normfilter == "SQUARE"):
                    if( (orbitparams.normMinlambda == dataparams.minlambda[band]) and (orbitparams.normMaxlambda == dataparams.maxlambda[band]) ):
                        databand = band
                else:
                    databand = band
        if( databand == 0 ):
            sys.exit("Could not find a matching data band in Normalize().")

        if( globalvar.verbose == "ON"):
            if( orbitparams.normfilter == "SQUARE"):
	          print("    Normalizing in the SQUARE {} {} bandpass.\n").format(
                                  orbitparams.normMinlambda, orbitparams.normMaxlambda)
            else:
	          print("    Normalizing in the {} filter.\n").format( orbitparams.normfilter)
             
        sum1 = 0.0
        sum2 = 0.0
        for i in range(1, dataparams.npoints[databand]):
            k = i / 10
            k = i - 10 * k
            if( k == 0 ):
                if( globalvar.verbose == "ON" ):
                    print("    phase number {}   phase =  {}\n").format(
                                     i, dataparams.phase[databand][i])
            FluxesAtPhase( dataparams.phase[databand][i], TotalFlux )
            weight = 1.0 / ( dataparams.standdev[databand][i] 
                       * dataparams.standdev[databand][i] )
            sum1 += weight * dataparams.flux[databand][i] * TotalFlux[calcband]
            sum2 += weight * TotalFlux[calcband] * TotalFlux[calcband]
        normfactor = sum1 / sum2
        for band in range(1, orbitparams.nbands):
            for iphase in range(0, orbitparams.maxpindex):
                globalvar.NormLC[band][iphase] = normfactor * globalvar.LCflux[band][iphase]

        dataparams.chisquare = 0.0
        for i in range(1, dataparams.npoints[databand]):
            FluxesAtPhase( dataparams.phase[databand][i], TotalFlux )
            weight = 1.0 / ( dataparams.standdev[databand][i] 
                       * dataparams.standdev[databand][i] )
            error = dataparams.flux[databand][i] - normfactor * TotalFlux[calcband]
            dataparams.chisquare += weight * error * error
        if( globalvar.verbose == "ON"): 
            print(" chi^2({}) = {}\n").format(dataparams.npoints[databand], 
	       dataparams.chisquare)
    else:
        sys.exit("Unrecognized normalization type in function Normalize().")

    return

def Irradiate():
    """
    This function calculates the heating due to irradiation.  Note:
       -- Disk irradiation is calculated first and includes 
             irradiation from star 1, the ADC, and the inner disk.
       -- Star 2 irradiation is calculated second and includes
             irradiation from star 1, the ADC, the inner disk, and the 
             (possibly heated) disk.
 
    In the current version the disk does not heat itself except
    for the inner disk.
 
    Also in the current version the ADC acts as if all its flux
    comes from two points located at y = +/1 adc.height above
    and below the disk; 1/2 of adc.L from each point.
 
    The code sacrificies efficiency for clarity in several places.
    """
    TDiskTold = [0 for i in range(0, MAXDISKTILES)]
    muA1toD = [0 for i in range(0, MAXDISKTILES)]
    transmit1toD = [0 for i in range(0, MAXDISKTILES)]
    DeltaT41toD = [0 for i in range(0, MAXDISKTILES)]
    muAidtoD = [0 for i in range(0, MAXDISKTILES)]
    transmitidtoD = [0 for i in range(0, MAXDISKTILES)]
    DeltaT4idtoD = [0 for i in range(0, MAXDISKTILES)]
    DeltaT4ADCtoD = [0 for i in range(0, MAXDISKTILES)]
    transmitADCtoD = [0 for i in range(0, MAXDISKTILES)]
    muAADCtoD = [0 for i in range(0, MAXDISKTILES)]
    T2Told = [0 for i in range(0, MAXDISKTILES)]
    muA1to2 = [0 for i in range(0, MAXDISKTILES)]
    transmit1to2 = [0 for i in range(0, MAXDISKTILES)]
    DeltaT41to2 = [0 for i in range(0, MAXDISKTILES)]
    muAidto2 = [0 for i in range(0, MAXDISKTILES)]
    transmitidto2 = [0 for i in range(0, MAXDISKTILES)]
    DeltaT4idto2 = [0 for i in range(0, MAXDISKTILES)]
    DeltaT4Dto2 = [0 for i in range(0, MAXDISKTILES)]    
    muAADCto2 = [0 for i in range(0, MAXDISKTILES)]    
    DeltaT4ADCto2 = [0 for i in range(0, MAXDISKTILES)]    
    transmitADCto2 = [0 for i in range(0, MAXDISKTILES)]    


    if( globalvar.verbose == "ON" ): 
        print(" Begin heating by irradiation.\n")
   
    if( flowcontrol.disk == "ON" ):

        if (flowcontrol.star1 == "ON"):
            if( globalvar.verbose == "ON" ): 
                print("   Begin heating the outer disk by star 1.\n")

            for iDisktile in range(1, wholediskpars.Ntiles):
                TDiskTold[iDisktile] = globalvar.TDiskT[iDisktile]
            start = CartVector()            
            start.x = 0.0
            start.y = 0.0
            start.z = systemparams.a
            for iDisktile in range(1, wholediskpars.Ntiles):
                end = CartVector()                
                end.x = globalvar.TDiskx[iDisktile]
                end.y = globalvar.TDisky[iDisktile]
                end.z = globalvar.TDiskz[iDisktile]
                delta = CartVector()
                delta.x = end.x - start.x
                delta.y = end.y - start.y
                delta.z = end.z - start.z
                d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                direction = CartVector()                
                direction.x = delta.x / d
                direction.y = delta.y / d
                direction.z = delta.z / d
                muA1toD[iDisktile] = CartDotProd( direction, 
                                               globalvar.TDisknormCart[iDisktile])
                if( muA1toD[iDisktile] < 0.0 ):
                    DeltaT41toD[iDisktile] = math.abs( Star1.Star1TotFlux( d )
                                            * muA1toD[iDisktile]
                                            / SIGMA )
                    transmit1toD[iDisktile] = Transmission( start, end )
                else:
                    DeltaT41toD[iDisktile] = 0.0
                    transmit1toD[iDisktile] = 0.0
                globalvar.TDiskT4[iDisktile] = globalvar.TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT41toD[iDisktile] * transmit1toD[iDisktile];
                globalvar.TDiskT[iDisktile] = pow( globalvar.TDiskT4[iDisktile], 0.25 )
            if( flowcontrol.diagnostics == "INSPECTHEATING"):
                InspectHeatDiskBy1(TDiskTold, muA1toD, DeltaT41toD, transmit1toD )
        if (flowcontrol.innerdisk == "ON"):
            if( globalvar.verbose == "ON" ): 
                print("   Begin heating the outer disk by the inner disk.\n")

            for iDisktile in range(1, wholediskpars.Ntiles):
                TDiskTold[iDisktile] = globalvar.TDiskT[iDisktile]
                start.x = 0.0
                start.y = 0.0
                start.z = systemparams.a
                for iDisktile in range(1, wholediskpars.Ntiles):
                    end.x = globalvar.TDiskx[iDisktile]
                    end.y = globalvar.TDisky[iDisktile]
                    end.z = globalvar.TDiskz[iDisktile]
                    delta.x = end.x - start.x
                    delta.y = end.y - start.y
                    delta.z = end.z - start.z
                    d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                    direction.x = delta.x / d
                    direction.y = delta.y / d
                    direction.z = delta.z / d
                    muAidtoD[iDisktile] = CartDotProd( direction, 
                                               globalvar.TDisknormCart[iDisktile])
                    if( muAidtoD[iDisktile] < 0.0 ):
                        DeltaT4idtoD[iDisktile] = math.abs( InnerDiskTotFlux( d ) 
                                         * direction.y
                                         * muAidtoD[iDisktile] 
                                         / SIGMA )
                        transmitidtoD[iDisktile] = Transmission( start, end )
                    else:
                        DeltaT4idtoD[iDisktile] = 0.0
                        transmitidtoD[iDisktile] = 0.0
                    globalvar.TDiskT4[iDisktile] = globalvar.TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4idtoD[iDisktile] * transmitidtoD[iDisktile]
                    globalvar.TDiskT[iDisktile] = pow( globalvar.TDiskT4[iDisktile], 0.25 )
                if( flowcontrol.diagnostics == "INSPECTHEATING"):
                    InspectHeatDiskByID(TDiskTold, muAidtoD, DeltaT4idtoD, 
                                   transmitidtoD )

            if (flowcontrol.adc == "ON"):
                 if( globalvar.verbose == "ON" ): 
                     print("   Begin heating the outer disk by the ADC.\n")

                 for iDisktile in range(1, wholediskpars.Ntiles):
                     TDiskTold[iDisktile] = globalvar.TDiskT[iDisktile]
                 
                 start.x = 0.0
                 start.y = adcpars.height
                 start.z = systemparams.a
                 for iDisktile in range(1, wholediskpars.Ntiles):
                     end.x = globalvar.TDiskx[iDisktile]
                     end.y = globalvar.TDisky[iDisktile]
                     end.z = globalvar.TDiskz[iDisktile]
                     delta.x = end.x - start.x
                     delta.y = end.y - start.y
                     delta.z = end.z - start.z
                     d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                     direction.x = delta.x / d
                     direction.y = delta.y / d
                     direction.z = delta.z / d
                     muAADCtoD[iDisktile] = CartDotProd( direction, 
                                               globalvar.TDisknormCart[iDisktile])
                     if( muAADCtoD[iDisktile] < 0.0 ):
                         DeltaT4ADCtoD[iDisktile] = math.abs( ADCTotFlux( d ) 
                                         * muAADCtoD[iDisktile] 
                                         / SIGMA )
                         transmitADCtoD[iDisktile] = Transmission( start, end )
                     else:
                         DeltaT4ADCtoD[iDisktile] = 0.0
                         transmitADCtoD[iDisktile] = 0.0
                     globalvar.TDiskT4[iDisktile] = globalvar.TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4ADCtoD[iDisktile] * transmitADCtoD[iDisktile]
                     globalvar.TDiskT[iDisktile] = pow( globalvar.TDiskT4[iDisktile], 0.25 )
                 if( flowcontrol.diagnostics == "INSPECTHEATING"):
                     InspectHeatDiskByADC( "TOP", TDiskTold, muAADCtoD, 
                                  DeltaT4ADCtoD,  transmitADCtoD )
                 start.x = 0.0
                 start.y = -adcpars.height
                 start.z = systemparams.a
                 for iDisktile in range(1, wholediskpars.Ntiles):
                     end.x = globalvar.TDiskx[iDisktile]
                     end.y = globalvar.TDisky[iDisktile]
                     end.z = globalvar.TDiskz[iDisktile]
                     delta.x = end.x - start.x
                     delta.y = end.y - start.y
                     delta.z = end.z - start.z
                     d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                     direction.x = delta.x / d
                     direction.y = delta.y / d
                     direction.z = delta.z / d
                     muAADCtoD[iDisktile] = CartDotProd( direction, 
                                               globalvar.TDisknormCart[iDisktile])
                     if( muAADCtoD[iDisktile] < 0.0 ):
                         DeltaT4ADCtoD[iDisktile] = math.abs( ADCTotFlux( d ) 
                                         * muAADCtoD[iDisktile] 
                                         / SIGMA )
                         transmitADCtoD[iDisktile] = Transmission( start, end )
                     else:
                         DeltaT4ADCtoD[iDisktile] = 0.0
                         transmitADCtoD[iDisktile] = 0.0
                     globalvar.TDiskT4[iDisktile] = globalvar.TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4ADCtoD[iDisktile] * transmitADCtoD[iDisktile]
                     globalvar.TDiskT[iDisktile] = pow( globalvar.TDiskT4[iDisktile], 0.25 )

                 if( flowcontrol.diagnostics == "INSPECTHEATING"):
                     InspectHeatDiskByADC( "BOTTOM", TDiskTold, muAADCtoD, 
                                  DeltaT4ADCtoD,  transmitADCtoD )

            for band in range(1, orbitparams.nbands):
                if( orbitparams.filter[band] == "SQUARE"):
                    for iDisktile in range(1, wholediskpars.ntiles):
                        globalvar.TDiskI[band][iDisktile] = BBSquareIntensity( globalvar.TDiskT[iDisktile], 
                                     orbitparams.minlambda[band],
                                     orbitparams.maxlambda[band])
                else:
                    for iDisktile in range(1, wholediskpars.Ntiles):
                        globalvar.TDiskI[band][iDisktile] = BBFilterIntensity( 
                                     globalvar.TDiskT[iDisktile], orbitparams.filter[band])

        if( flowcontrol.star2 == "ON" ):

            if( flowcontrol.star1 == "ON" ):
                if( globalvar.verbose == "ON" ): 
                    print("   Begin heating star 2 by star 1.\n")

                for i2tile in range(1, Star2.Ntiles):
                    T2Told[i2tile] = globalvar.T2T[i2tile]
                start.x = 0.0
                start.y = 0.0
                start.z = systemparams.a
                for i2tile in range(1, Star2.Ntiles):
                    end.x = globalvar.T2x[i2tile]
                    end.y = globalvar.T2y[i2tile]
                    end.z = globalvar.T2z[i2tile]
                    delta.x = end.x - start.x
                    delta.y = end.y - start.y
                    delta.z = end.z - start.z
                    d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                    direction.x = delta.x / d
                    direction.y = delta.y / d
                    direction.z = delta.z / d
                    muA1to2[i2tile] = CartDotProd( direction, globalvar.T2normCart[i2tile] )
                    if( muA1to2[i2tile] < 0.0 ):
                        DeltaT41to2[i2tile] = math.abs( Star1.Star1TotFlux( d )
                                           * muA1to2[i2tile] 
                                           / SIGMA )
                        transmit1to2[i2tile] = Transmission( start, end)
                    else:
                        DeltaT41to2[i2tile] = 0.0
                        transmit1to2[i2tile] = 0.0
                    summation = pow( globalvar.T2T[i2tile], 4.0 ) + Star2.albedo * DeltaT41to2[i2tile] * transmit1to2[i2tile]
                    globalvar.T2T[i2tile] = pow( summation, 0.25 )

                if( flowcontrol.diagnostics == "INSPECTHEATING"):
                    InspectHeat2By1(T2Told, muA1to2, DeltaT41to2, transmit1to2 )
 
            if( flowcontrol.disk == "ON" ):
                if( flowcontrol.innerdisk == "ON"):
                    if( globalvar.verbose == "ON" ): 
                        print("   Begin heating star 2 by the inner disk.\n")
              
                    for i2tile in range(1, Star2.Ntiles):
                        T2Told[i2tile] = globalvar.T2T[i2tile]
                    start.x = 0.0
                    start.y = 0.0
                    start.z = systemparams.a
                    for i2tile in range(1, Star2.Ntiles):
                        end.x = globalvar.T2x[i2tile]
                        end.y = globalvar.T2y[i2tile]
                        end.z = globalvar.T2z[i2tile]
                        delta.x = end.x - start.x
                        delta.y = end.y - start.y
                        delta.z = end.z - start.z
                        d =  math.sqrt( delta.x*delta.x + delta.y*delta.y 
                                                   + delta.z*delta.z )
                        direction.x = delta.x / d
                        direction.y = delta.y / d
                        direction.z = delta.z / d
                        muAidto2[i2tile] = CartDotProd( direction, globalvar.T2normCart[i2tile] )
                        if( muAidto2[i2tile] < 0.0 ):
                            DeltaT4idto2[i2tile] = math.abs( InnerDiskTotFlux( d )
                                               * direction.y
                                               * muAidto2[i2tile]
                                               / SIGMA )
                            transmitidto2[i2tile] = Transmission( start, end)
                        else:
                           DeltaT4idto2[i2tile] = 0.0
                           transmitidto2[i2tile] = 0.0
                        summation = pow( globalvar.T2T[i2tile], 4.0 ) + Star2.albedo * DeltaT4idto2[i2tile] * transmitidto2[i2tile]
                        globalvar.T2T[i2tile] = pow( summation, 0.25 )
                    if( flowcontrol.diagnostics == "INSPECTHEATING"):
                        InspectHeat2ByID(T2Told, muAidto2, DeltaT4idto2, 
                                       transmitidto2 )

                if( globalvar.verbose == "ON" ): 
	              print("   Begin heating star 2 by the outer disk.\n")

                for i2tile in range(1, Star2.Ntiles):
                    T2Told[i2tile] = globalvar.T2T[i2tile]
                for i2tile in range(1, Star2.Ntiles):
                    if( globalvar.verbose == "ON"):
                        k = i2tile / 100
                        k = i2tile - 100 * k
                        if( k == 0 ):
                            print("      heating star 2 tile number %5ld\n", i2tile)
                    DeltaT4Dto2[i2tile] = 0.0
                    end.x = globalvar.T2x[i2tile]
                    end.y = globalvar.T2y[i2tile]
                    end.z = globalvar.T2z[i2tile]
                    for iDisktile in range(1, wholediskpars.Ntiles):
                        start.x = globalvar.TDiskx[iDisktile]
                        start.y = globalvar.TDisky[iDisktile]
                        start.z = globalvar.TDiskz[iDisktile]
                        vectord = CartVector()                        
                        vectord.x = end.x - start.x
                        vectord.y = end.y - start.y
                        vectord.z = end.z - start.z
                        muEprime = CartDotProd( vectord, globalvar.TDisknormCart[iDisktile] )
                        if( muEprime > 0.0 ):
                            muAprime = CartDotProd( vectord, globalvar.T2normCart[i2tile] )
                            if( muAprime < 0.0 ):
                                dsquare =   vectord.x * vectord.x + vectord.y * vectord.y + vectord.z * vectord.z
                                d = math.sqrt( dsquare )
                                muE = muEprime / d
                                muA = muAprime / d
                                delT4disk = -( globalvar.TDiskT4[iDisktile] / math.pi) * (muE * muA / dsquare) * globalvar.TDiskdS[iDisktile]
                                delT4disk *= Transmission( start, end )
                                DeltaT4Dto2[i2tile] += delT4disk 
                            summation = pow( globalvar.T2T[i2tile], 4.0) + Star2.albedo * DeltaT4Dto2[i2tile]
                            globalvar.T2T[i2tile] = pow( sum, 0.25 )
                        if( flowcontrol.diagnostics == "INSPECTHEATING"):
                            InspectHeat2ByDisk( T2Told, DeltaT4Dto2, globalvar.T2T)

                    if (flowcontrol.adc == "ON"):
                        if( globalvar.verbose == "ON" ): 
                            print("   Begin heating star 2 by the ADC.\n")
                        for i2tile in range(1, Star2.Ntiles):
                            T2Told[i2tile] = globalvar.T2T[i2tile]

                        start.x = 0.0
                        start.y = adcpars.height
                        start.z = systemparams.a
                        for i2tile in range(1, Star2.Ntiles):
                            end.x = globalvar.T2x[i2tile]
                            end.y = globalvar.T2y[i2tile]
                            end.z = globalvar.T2z[i2tile]
                            delta.x = end.x - start.x
                            delta.y = end.y - start.y
                            delta.z = end.z - start.z
                            d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                            direction.x = delta.x / d
                            direction.y = delta.y / d
                            direction.z = delta.z / d
                            muAADCto2[i2tile] = CartDotProd( direction, globalvar.T2normCart[i2tile] )
                            if( muAADCto2[i2tile] < 0.0 ):
                                DeltaT4ADCto2[i2tile] = math.abs( ADCTotFlux( d ) 
                                             * muAADCto2[i2tile] 
                                             / SIGMA )
                                transmitADCto2[i2tile] = Transmission( start, end )
                            else:
                                DeltaT4ADCto2[i2tile] = 0.0
                                transmitADCto2[i2tile] = 0.0
                                summation = pow( globalvar.T2T[i2tile], 4.0) + Star2.albedo * DeltaT4ADCto2[i2tile] * transmitADCto2[i2tile]
                            globalvar.T2T[i2tile] = pow( sum, 0.25 )
                            if( flowcontrol.diagnostics == "INSPECTHEATING"):
                                InspectHeat2ByADC( "TOP", T2Told, muAADCto2, 
                                  DeltaT4ADCto2,  transmitADCto2 )
                        start.x = 0.0
                        start.y = -adcpars.height
                        start.z = systemparams.a
                        for i2tile in range(1, Star2.Ntiles):
                            end.x = globalvar.T2x[i2tile]
                            end.y = globalvar.T2y[i2tile]
                            end.z = globalvar.T2z[i2tile]
                            delta.x = end.x - start.x
                            delta.y = end.y - start.y
                            delta.z = end.z - start.z
                            d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                            direction.x = delta.x / d
                            direction.y = delta.y / d
                            direction.z = delta.z / d
                            muAADCto2[i2tile] = CartDotProd( direction, globalvar.T2normCart[i2tile] )
                            if( muAADCto2[i2tile] < 0.0 ):
                                DeltaT4ADCto2[i2tile] = math.abs( ADCTotFlux( d ) 
                                             * muAADCto2[i2tile] 
                                             / SIGMA )
                                transmitADCto2[i2tile] = Transmission( start, end )
                            else:
                                DeltaT4ADCto2[i2tile] = 0.0
                                transmitADCto2[i2tile] = 0.0
                            summation = pow( globalvar.T2T[i2tile], 4.0)+ Star2.albedo * DeltaT4ADCto2[i2tile] * transmitADCto2[i2tile]
                        globalvar.T2T[i2tile] = pow( sum, 0.25 )

                    if( flowcontrol.diagnostics == "INSPECTHEATING"):
                        InspectHeat2ByADC( "BOTTOM", T2Told, muAADCto2, 
                                  DeltaT4ADCto2,  transmitADCto2 )

                for band in range(1, orbitparams.nbands):
                    if( orbitparams.filter[band] == "SQUARE"):
                        for i2tile in range(1, Star2.Ntiles):
                            globalvar.T2I[band][i2tile] = BBSquareIntensity( globalvar.T2T[i2tile], 
                                                 orbitparams.minlambda[band], 
                                                 orbitparams.maxlambda[band])
                    else:
                        for i2tile in range(1, Star2.Ntiles):
                            if( (globalvar.T2T[i2tile] > globalvar.IperpT[globalvar.maxIperpTindex]) or (globalvar.T2T[i2tile] < globalvar.IperpT[0]) ):
                                globalvar.T2I[band][i2tile] = BBFilterIntensity( globalvar.T2T[i2tile], 
                                                     orbitparams.filter[band])
                            else:
                                globalvar.T2I[band][i2tile] = Star2.GetIperp( globalvar.T2T[i2tile], globalvar.T2logg[i2tile], 
                                           orbitparams.filter[band])

                if( flowcontrol.diagnostics == "INSPECTHEATING"):
                    sys.exit("Quit after INSPECTHEATING.")

                return

def EscapeFraction( start, direction):
    """
    This function traces light rays through the binary system to
    determine if the light ray escapes the system and contributes
    to the light curve.  It also determines the fraction of the 
    intensity at the beginning of the ray that survives to escape
    from the binary.  In this version the escape fraction is 0 or 1.
 
    Note that the calculation starts after the ray has already
    traversed a few steplengths.  This avoids some of the
    problems introduced by pixelization.
    """
    transmit = 1.0
    steplength = 0.8 * XYGrid.deltal
    deltaray = CartVector()
    deltaray.x = steplength * direction.x
    deltaray.y = steplength * direction.y
    deltaray.z = steplength * direction.z

    ray = CartVector()
    ray.x = start.x + 3.0 * deltaray.x
    ray.y = start.y + 3.0 * deltaray.y
    ray.z = start.z + 3.0 * deltaray.z
    while(True):
        ray.x = ray.x + deltaray.x
        if( ray.x >= XYGrid.xmax ):
	       break
        if( ray.x <= XYGrid.xmin ):
             break
        ray.y = ray.y + deltaray.y
        if( ray.y >= XYGrid.ymax ):
	       break
        if( ray.y <= XYGrid.ymin ):
	       break
        ray.z = ray.z + deltaray.z
        if( ray.z >= XYGrid.zmax ):
            break
        if( ray.z <= XYGrid.zmin ):
            break
        ix = 1.5 + ( ray.x - XYGrid.xmin ) / XYGrid.deltax
        iz = 1.5 + ( ray.z - XYGrid.zmin ) / XYGrid.deltaz
        if( ray.y < XYGrid.Topy[ix][iz] ):
            if( ray.y > XYGrid.Bottomy[ix][iz] ):
                transmit = 0.0
                return( transmit )

    return( transmit )

def Transmission( start, end):
    """
    This function traces light rays between two points in the binary
    system.  It returns the fraction of the intensity at the beginning
    of the ray that survives to the end of the ray.  In this first
    version the transmission is either 1.0 or 0.0.
 
    Note that the function jumps the ray a few stepsizes when if
    first begins to propagate the ray, and the ray is assumed to have 
    arrived at the end of its path if it gets to within 3 pixels of
    the end point.  This is reasonable since the long dimension of a
    typical tile is several pixels long.  These measures ameliorate 
    but do not entirely  eliminate the nastier effects of pixelization.
 
    NOTE:  Program XRbinary spends the vast majority of its time
           in this function.
    """
    ixEnd = 1.5 + ( end.x - XYGrid.xmin ) / XYGrid.deltax
    izEnd = 1.5 + ( end.z - XYGrid.zmin ) / XYGrid.deltaz

    delta = CartVector()
    delta.x = end.x - start.x
    delta.y = end.y - start.y
    delta.z = end.z - start.z
    length = math.sqrt( delta.x * delta.x + delta.y * delta.y + delta.z * delta.z )
    direction = CartVector()
    direction.x = delta.x / length
    direction.y = delta.y / length
    direction.z = delta.z / length
    nsteps = length / XYGrid.deltal
    if( nsteps <= 3 ):
        transmit = 1.0
        return( transmit )
    stepsize = length / nsteps
    deltaray = CartVector()
    deltaray.x = stepsize * direction.x
    deltaray.y = stepsize * direction.y
    deltaray.z = stepsize * direction.z

    ray = CartVector()
    ray.x = start.x + 3.0 * deltaray.x
    ray.y = start.y + 3.0 * deltaray.y
    ray.z = start.z + 3.0 * deltaray.z

    InverseDeltax = 1.0 / XYGrid.deltax
    InverseDeltaz = 1.0 / XYGrid.deltaz
    while(True):
        ray.x += deltaray.x
        ray.y += deltaray.y
        ray.z += deltaray.z
        ix = 1.5 + ( ray.x - XYGrid.xmin ) * InverseDeltax
        iz = 1.5 + ( ray.z - XYGrid.zmin ) * InverseDeltaz
        if( abs( ixEnd - ix ) <= 3 ):
            if( abs( izEnd - iz ) <= 3 ):
                transmit = 1.0
                break
        if( ray.y < XYGrid.Topy[ix][iz] ):
            if( ray.y > XYGrid.Bottomy[ix][iz] ):
                transmit = 0.0
                break

    return( transmit )

def MakeYlimits():
    """
    This function finds the surfaces Grid.Topy[i][j] and Grid.Bottomy[i][j]
    defined by the highest and lowest values of y of the components 
    in the system.  The function also finds maximum and minimum
    values of x, y, z needed so the edges of the grid cover
    the stars and disk.
    """
    if( globalvar.verbose == "ON"):
        print(" Begin making Ylimits grid.\n")
    margin = 0.04 * systemparams.a

    XYGrid.xmin = 0.0
    XYGrid.xmax = 0.0
    XYGrid.ymin = 0.0
    XYGrid.ymax = 0.0
    XYGrid.zmin = 0.0
    XYGrid.zmax = 0.0
    if(  flowcontrol.star1 == "ON" ) or (flowcontrol.innerdisk == "ON" ) or (flowcontrol.adc == "ON" ):
        XYGrid.zmin = systemparams.a
        XYGrid.zmax = systemparams.a
    if( flowcontrol.adc == "ON" ):
        XYGrid.ymin = -adcpars.height
        XYGrid.ymax =  adcpars.height

    if( flowcontrol.star2 == "ON" ):
        for itile in range(1, Star2.Ntiles):
            if( globalvar.T2x[itile] < XYGrid.xmin ): 
                XYGrid.xmin = globalvar.T2x[itile]
            if( globalvar.T2x[itile] > XYGrid.xmax ):
                XYGrid.xmax = globalvar.T2x[itile]
            if( globalvar.T2y[itile] < XYGrid.ymin ):
                XYGrid.ymin = globalvar.T2y[itile]
            if( globalvar.T2y[itile] > XYGrid.ymax ):
                XYGrid.ymax = globalvar.T2y[itile]
            if( globalvar.T2z[itile] < XYGrid.zmin ):
                XYGrid.zmin = globalvar.T2z[itile]
            if( globalvar.T2z[itile] > XYGrid.zmax ):
                XYGrid.zmax = globalvar.T2z[itile]
    if( flowcontrol.disk == "ON" ):
        for itile in range(1, wholediskpars.Ntiles):
            if( globalvar.TDiskx[itile] < XYGrid.xmin ):
                XYGrid.xmin = globalvar.TDiskx[itile]
            if( globalvar.TDiskx[itile] > XYGrid.xmax ):
                XYGrid.xmax = globalvar.TDiskx[itile]
            if( globalvar.TDisky[itile] < XYGrid.ymin ):
                XYGrid.ymin = globalvar.TDisky[itile]
            if( globalvar.TDisky[itile] > XYGrid.ymax ):
                XYGrid.ymax = globalvar.TDisky[itile]
            if( globalvar.TDiskz[itile] < XYGrid.zmin ):
                XYGrid.zmin = globalvar.TDiskz[itile]
            if( globalvar.TDiskz[itile] > XYGrid.zmax ):
                XYGrid.zmax = globalvar.TDiskz[itile]

    XYGrid.xmin -= margin
    XYGrid.xmax += margin
    XYGrid.ymin -= margin
    XYGrid.ymax += margin
    XYGrid.zmin -= margin
    XYGrid.zmax += margin

    XYGrid.Nxtiles = GRIDXTILES - 1
    XYGrid.Nztiles = GRIDZTILES - 1
    XYGrid.deltax = (XYGrid.xmax - XYGrid.xmin) / ( XYGrid.Nxtiles - 1 )
    XYGrid.deltaz = (XYGrid.zmax - XYGrid.zmin) / ( XYGrid.Nztiles - 1 )
    if( (XYGrid.deltax <=0.0) or (XYGrid.deltaz <= 0.0) ):
        sys.exit("Either Grid.deltax or Grid.deltaz equals zero.")
    XYGrid.deltal = math.sqrt( XYGrid.deltax*XYGrid.deltax + XYGrid.deltaz*XYGrid.deltaz)
  
    for ix in range(1, XYGrid.Nxtiles):
        for iz in range(1, XYGrid.Nztiles):
            XYGrid.Topy[ix][iz]    = 0.0
            XYGrid.Bottomy[ix][iz] = 0.0
            x = XYGrid.xmin + (ix - 1) * XYGrid.deltax
            z = XYGrid.zmin + (iz - 1) * XYGrid.deltaz
            if( z < systemparams.rL1 ):
                if( flowcontrol.star2 == "ON" ):
	               XYGrid.Topy[ix][iz] = Star2.Star2TopY( x, z )
	               XYGrid.Bottomy[ix][iz] = -XYGrid.Topy[ix][iz]
            else:
                if( flowcontrol.disk == "ON" ):
                    rho = math.sqrt( x*x + (z - systemparams.a)*(z - systemparams.a) )
                    coszeta = (z - systemparams.a) / rho
                    zeta = math.acos( coszeta )
                    if( x < 0.0 ):
                        zeta = math.pi*2 - zeta
                    a = RhoToA( rho, zeta)
                    XYGrid.Topy[ix][iz] = DiskTopH( a, zeta)
                    XYGrid.Bottomy[ix][iz] = DiskBottomH( a, zeta)

    if( flowcontrol.diagnostics == "INSPECTYLIMITS"):
        InspectYlimits()
        sys.exit("Quit in MakeYlimits after INSPECTYLIMITS.")

    return
