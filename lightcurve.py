# -*- coding: utf-8 -*-
"""

"""

import math
import sys
from .diskflux import maindisk
from .star1 import Star1
from .star2 import Star2
from .parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, XYGrid, dataparams, ReadInput

def MakeLightCurves():
    """
    This function makes the orbital light curve.
    """
    if( flowcontrol.thirdlight == "ON" ):
        if( verbose == "ON" ):
            print(" Calculating third light fluxes.\n")
        ThirdLight()
    if( verbose == "ON" ): 
        print(" Begin calculating the light curves.\n")

    for iphase in range(0, orbitparams.maxpindex):
        LCphase[iphase] = orbitparams.phasemin + iphase * orbitparams.deltaphase
        k = iphase / 10
        k = iphase - 10 * k
        if( k == 0 ):
            if( verbose == "ON" ):
                print("    phase number {}   phase =  {}\n").format( iphase, 
                                                         LCphase[iphase])
        FluxesAtPhase( LCphase[iphase], TotalFlux )
        for band in range(1, orbit.nbands):
            LCflux[band][iphase] = TotalFlux[band]

    if( orbitparams.normalize == "OFF" ):
        if( verbose == "ON" ):
	       print(" Normalizing the light curves.\n")
        Normalize()
    else:
        for band in range(1, orbitparams.nbands):
            for iphase in range(0, orbitparams.maxpindex):
	           NormLC[band][iphase] = 0.0

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
    sunvector.x = -1.0 * math.sin( systemparams.i ) * math.sin( calcphase * 2*math.pi )
    sunvector.y =        math.cos( systemparams.i )
    sunvector.z = -1.0 * math.sin( systemparams.i ) * math.cos( calcphase * 2*math.pi )

    if( flowcontrol.star1 == "ON" ):
        start.x = 0.0
        start.y = 0.0
        start.z = systemparams.a
        Star1Escape = EscapeFraction( start, sunvector)

    if( flowcontrol.star2 == "ON"):
        for itile in range(1, Star2.Ntiles):
            T2mu[itile] = CartDotProd( sunvector, T2normCart[itile] )
            if( T2mu[itile] > 0.0 ):
                start.x = T2x[itile]
                start.y = T2y[itile]
                start.z = T2z[itile]
                T2Escape[itile] = EscapeFraction( start, sunvector)
            else:
                T2Escape[itile] = 0.0

    if( flowcontrol.disk == "ON"):
        for itile in range(1, wholediskpars.Ntiles):
            TDiskmu[itile] = CartDotProd( sunvector, TDisknormCart[itile] )
            if( TDiskmu[itile] > 0.0 ):
                start.x = TDiskx[itile]
                start.y = TDisky[itile]
                start.z = TDiskz[itile]
                TDiskEscape[itile] = EscapeFraction( start, sunvector)
            else:
                TDiskEscape[itile] = 0.0
        if( flowcontrol.innerdisk == "ON" ):
            start.x = 0.0
            start.y = 0.0
            start.z = syspars.a
            InnerDiskEscape = EscapeFraction( start, sunvector)


    for band in range(1, orbitparams.nbands):
        TotalFlux[band] = 0.0

        if( flowcontrol.star1 == "ON" ):
            Star1Emitted = Star1Flambda( orbit.filter[band], 
                                   orbitparams.minlambda[band],
                                   orbitparams.maxlambda[band] )
            star1flux = Star1Emitted * Star1Escape
            TotalFlux[band] += star1flux

        if( flowcontrol.star2 == "ON" ):
            star2flux = 0.0
            if( orbitparams.filter[band] == "SQUARE"):
                for itile in range(1, Star2.Ntiles):
                    if( T2mu[itile] < 0.0 ):
                       T2Emitted[itile] = 0.0
                    else:
   	                  T2Emitted[itile] =  T2I[band][itile]*T2mu[itile] * T2dS[itile]
                    star2flux += T2Emitted[itile] * T2Escape[itile]
            else:
                for itile in range(1, Star2.Ntiles):
                    if( T2mu[itile] < 0.0 ):
		             T2Emitted[itile] = 0.0
                    else:
                        if( (T2T[itile] > IperpT[maxIperpTindex]) or (T2T[itile] < IperpT[0]) ):
                            T2Emitted[itile] = T2I[band][itile]*T2mu[itile] * T2dS[itile]
                        else:
                            T2Emitted[itile] = T2I[band][itile]* ClaretHmu( T2T[itile],T2logg[itile],orbit.filter[band],T2mu[itile] )* T2mu[itile] * T2dS[itile]
                    star2flux += T2Emitted[itile] * T2Escape[itile];
            TotalFlux[band] += star2flux
        if( flowcontrol.disk == "ON" ):
            diskflux = 0.0
            for itile in range(1, wholediskpars.Ntiles):
                if( TDiskmu[itile] < 0.0 ):
	              TDiskEmitted[itile] = 0.0
                else:
   	              TDiskEmitted[itile] = TDiskI[band][itile] * TDiskmu[itile] * TDiskdS[itile]
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
    if( orbitparams.normalize == "MAXVALUE"):
        for band in range(1, orbitparams.nbands):
            maxflux = 0.0
            for iphase in range(0, orbitparams.maxpindex):
                if( LCflux[band][iphase] > maxflux ):
                    maxflux = LCflux[band][iphase]
            if( maxflux == 0.0 ):
                maxflux = 1.0
            normfactor = orbitparams.normvalue / maxflux
            for iphase in range(0, orbitparams.maxpindex):
                NormLC[band][iphase] = normfactor * LCflux[band][iphase]
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

        if( verbose == "ON"):
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
                if( verbose == "ON" ):
                    print("    phase number {}   phase =  {}\n").format(
                                     i, dataparams.phase[databand][i])
            FluxesAtPhase( dataparams.phase[databand][i], TotalFlux )
            weight = 1.0 / ( dataparams.standdev[databand][i] 
                       * dataparams.standdev[databand][i] )
            sum1 += weight * dataparams.flux[databand][i] * TotalFlux[calcband]
            sum2 += weight * TotalFlux[calcband] * TotalFlux[calcband]
        normfactor = sum1 / sum2
        for band in range(1, orbitparams.nbands):
            for iphase in range(0, maxpindex):
                NormLC[band][iphase] = normfactor * LCflux[band][iphase]

        dataparams.chisquare = 0.0
        for i in range(1, dataparams.npoints[databand]):
            FluxesAtPhase( dataparams.phase[databand][i], TotalFlux )
            weight = 1.0 / ( dataparams.standdev[databand][i] 
                       * dataparams.standdev[databand][i] )
            error = dataparams.flux[databand][i] - normfactor * TotalFlux[calcband]
            data.chisquare += weight * error * error
        if( verbose == "ON"): 
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
    if( verbose == "ON" ): 
        print(" Begin heating by irradiation.\n")
   
    if( flowcontrol.disk == "ON" ):

        if (flowcontrol.star1 == "ON"):
            if( verbose == "ON" ): 
                print("   Begin heating the outer disk by star 1.\n")

            for iDisktile in range(1, wholediskpars.Ntiles):
                TDiskTold[iDisktile] = TDiskT[iDisktile]
            start.x = 0.0
            start.y = 0.0
            start.z = systemparams.a
            for iDisktile in range(1, wholediskpars.Ntiles):
                end.x = TDiskx[iDisktile]
                end.y = TDisky[iDisktile]
                end.z = TDiskz[iDisktile]
                delta.x = end.x - start.x
                delta.y = end.y - start.y
                delta.z = end.z - start.z
                d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                direction.x = delta.x / d
                direction.y = delta.y / d
                direction.z = delta.z / d
                muA1toD[iDisktile] = CartDotProd( direction, 
                                               TDisknormCart[iDisktile])
                if( muA1toD[iDisktile] < 0.0 ):
                    DeltaT41toD[iDisktile] = math.abs( Star1TotFlux( d )
                                            * muA1toD[iDisktile]
                                            / SIGMA )
                    transmit1toD[iDisktile] = Transmission( start, end )
                else:
                    DeltaT41toD[iDisktile] = 0.0
                    transmit1toD[iDisktile] = 0.0
                TDiskT4[iDisktile] = TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT41toD[iDisktile] * transmit1toD[iDisktile];
                TDiskT[iDisktile] = pow( TDiskT4[iDisktile], 0.25 )
            if( flowcontrol.diagnostics == "INSPECTHEATING"):
                InspectHeatDiskBy1(TDiskTold, muA1toD, DeltaT41toD, transmit1toD )
        if (flowcontrol.innerdisk == "ON"):
            if( verbose == "ON" ): 
                print("   Begin heating the outer disk by the inner disk.\n")

            for iDisktile in range(1, Ntiles):
                TDiskTold[iDisktile] = TDiskT[iDisktile]
                start.x = 0.0
                start.y = 0.0
                start.z = systemparams.a
                for iDisktile in range(1, wholediskpars.Ntiles):
                    end.x = TDiskx[iDisktile]
                    end.y = TDisky[iDisktile]
                    end.z = TDiskz[iDisktile]
                    delta.x = end.x - start.x
                    delta.y = end.y - start.y
                    delta.z = end.z - start.z
                    d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                    direction.x = delta.x / d
                    direction.y = delta.y / d
                    direction.z = delta.z / d
                    muAidtoD[iDisktile] = CartDotProd( direction, 
                                               TDisknormCart[iDisktile])
                    if( muAidtoD[iDisktile] < 0.0 ):
                        DeltaT4idtoD[iDisktile] = math.abs( InnerDiskTotFlux( d ) 
                                         * direction.y
                                         * muAidtoD[iDisktile] 
                                         / SIGMA )
                        transmitidtoD[iDisktile] = Transmission( start, end )
                    else:
                        DeltaT4idtoD[iDisktile] = 0.0
                        transmitidtoD[iDisktile] = 0.0
                    TDiskT4[iDisktile] = TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4idtoD[iDisktile] * transmitidtoD[iDisktile]
                    TDiskT[iDisktile] = pow( TDiskT4[iDisktile], 0.25 )
                if( flowcontrol.diagnostics == "INSPECTHEATING"):
                    InspectHeatDiskByID(TDiskTold, muAidtoD, DeltaT4idtoD, 
                                   transmitidtoD )

            if (flowcontrol.adc == "ON"):
                 if( verbose == "ON" ): 
                     print("   Begin heating the outer disk by the ADC.\n")

                 for iDisktile in range(1, wholediskpars.Ntiles):
                     TDiskTold[iDisktile] = TDiskT[iDisktile]
                 
                 start.x = 0.0
                 start.y = adc.height
                 start.z = systemparams.a
                 for iDisktile in range(1, Ntiles):
                     end.x = TDiskx[iDisktile]
                     end.y = TDisky[iDisktile]
                     end.z = TDiskz[iDisktile]
                     delta.x = end.x - start.x
                     delta.y = end.y - start.y
                     delta.z = end.z - start.z
                     d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                     direction.x = delta.x / d
                     direction.y = delta.y / d
                     direction.z = delta.z / d
                     muAADCtoD[iDisktile] = CartDotProd( direction, 
                                               TDisknormCart[iDisktile])
                     if( muAADCtoD[iDisktile] < 0.0 ):
                         DeltaT4ADCtoD[iDisktile] = math.abs( ADCTotFlux( d ) 
                                         * muAADCtoD[iDisktile] 
                                         / SIGMA )
                         transmitADCtoD[iDisktile] = Transmission( start, end )
                     else:
                         DeltaT4ADCtoD[iDisktile] = 0.0
                         transmitADCtoD[iDisktile] = 0.0
                     TDiskT4[iDisktile] = TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4ADCtoD[iDisktile] * transmitADCtoD[iDisktile]
                     TDiskT[iDisktile] = pow( TDiskT4[iDisktile], 0.25 )
                 if( flowcontrol.diagnostics == "INSPECTHEATING"):
                     InspectHeatDiskByADC( "TOP", TDiskTold, muAADCtoD, 
                                  DeltaT4ADCtoD,  transmitADCtoD )
                 start.x = 0.0
                 start.y = -adcpars.height
                 start.z = systemparams.a
                 for iDisktile in range(1, wholediskpars.Ntiles):
                     end.x = TDiskx[iDisktile]
                     end.y = TDisky[iDisktile]
                     end.z = TDiskz[iDisktile]
                     delta.x = end.x - start.x
                     delta.y = end.y - start.y
                     delta.z = end.z - start.z
                     d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                     direction.x = delta.x / d
                     direction.y = delta.y / d
                     direction.z = delta.z / d
                     muAADCtoD[iDisktile] = CartDotProd( direction, 
                                               TDisknormCart[iDisktile])
                     if( muAADCtoD[iDisktile] < 0.0 ):
                         DeltaT4ADCtoD[iDisktile] = math.abs( ADCTotFlux( d ) 
                                         * muAADCtoD[iDisktile] 
                                         / SIGMA )
                         transmitADCtoD[iDisktile] = Transmission( start, end )
                     else:
                         DeltaT4ADCtoD[iDisktile] = 0.0
                         transmitADCtoD[iDisktile] = 0.0
                     TDiskT4[iDisktile] = TDiskT4[iDisktile]+ wholediskpars.albedo * DeltaT4ADCtoD[iDisktile] * transmitADCtoD[iDisktile]
                     TDiskT[iDisktile] = pow( TDiskT4[iDisktile], 0.25 )

                 if( flowcontrol.diagnostics == "INSPECTHEATING"):
                     InspectHeatDiskByADC( "BOTTOM", TDiskTold, muAADCtoD, 
                                  DeltaT4ADCtoD,  transmitADCtoD )

            for band in range(1, orbitparams.nbands):
                if( orbitparams.filter[band] == "SQUARE"):
                    for iDisktile in range(1, wholediskpars.ntiles):
                        TDiskI[band][iDisktile] = BBSquareIntensity( TDiskT[iDisktile], 
                                     orbitparams.minlambda[band],
                                     orbitparams.maxlambda[band])
                else:
                    for iDisktile in range(1, wholediskpars.Ntiles):
                        TDiskI[band][iDisktile] = BBFilterIntensity( 
                                     TDiskT[iDisktile], orbitparams.filter[band])

        if( flowcontrol.star2 == "ON" ):

            if( flowcontrol.star1 == "ON" ):
                if( verbose == "ON" ): 
                    print("   Begin heating star 2 by star 1.\n")

                for i2tile in range(1, Star2.Ntiles):
                    T2Told[i2tile] = T2T[i2tile]
                start.x = 0.0
                start.y = 0.0
                start.z = systemparams.a
                for i2tile in range(1, Star2.Ntiles):
                    end.x = T2x[i2tile]
                    end.y = T2y[i2tile]
                    end.z = T2z[i2tile]
                    delta.x = end.x - start.x
                    delta.y = end.y - start.y
                    delta.z = end.z - start.z
                    d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                    direction.x = delta.x / d
                    direction.y = delta.y / d
                    direction.z = delta.z / d
                    muA1to2[i2tile] = CartDotProd( direction, T2normCart[i2tile] )
                    if( muA1to2[i2tile] < 0.0 ):
                        DeltaT41to2[i2tile] = math.abs( Star1TotFlux( d )
                                           * muA1to2[i2tile] 
                                           / SIGMA )
                        transmit1to2[i2tile] = Transmission( start, end)
                    else:
                        DeltaT41to2[i2tile] = 0.0
                        transmit1to2[i2tile] = 0.0
                    summation = pow( T2T[i2tile], 4.0 ) + Star2.albedo * DeltaT41to2[i2tile] * transmit1to2[i2tile]
                    T2T[i2tile] = pow( summation, 0.25 )

                if( flowcontrol.diagnostics == "INSPECTHEATING"):
                    InspectHeat2By1(T2Told, muA1to2, DeltaT41to2, transmit1to2 )
 
            if( flowcontrol.disk == "ON" ):
                if( flowcontrol.innerdisk == "ON"):
                    if( verbose == "ON" ): 
                        print("   Begin heating star 2 by the inner disk.\n")
              
                    for i2tile in range(1, Star2.Ntiles):
                        T2Told[i2tile] = T2T[i2tile]
                    start.x = 0.0
                    start.y = 0.0
                    start.z = systemparams.a
                    for i2tile in range(1, Star2.Ntiles):
                        end.x = T2x[i2tile]
                        end.y = T2y[i2tile]
                        end.z = T2z[i2tile]
                        delta.x = end.x - start.x
                        delta.y = end.y - start.y
                        delta.z = end.z - start.z
                        d =  math.sqrt( delta.x*delta.x + delta.y*delta.y 
                                                   + delta.z*delta.z )
                        direction.x = delta.x / d
                        direction.y = delta.y / d
                        direction.z = delta.z / d
                        muAidto2[i2tile] = CartDotProd( direction, T2normCart[i2tile] )
                        if( muAidto2[i2tile] < 0.0 ):
                            DeltaT4idto2[i2tile] = math.abs( InnerDiskTotFlux( d )
                                               * direction.y
                                               * muAidto2[i2tile]
                                               / SIGMA )
                            transmitidto2[i2tile] = Transmission( start, end)
                        else:
                           DeltaT4idto2[i2tile] = 0.0
                           transmitidto2[i2tile] = 0.0
                        summation = pow( T2T[i2tile], 4.0 ) + Star2.albedo * DeltaT4idto2[i2tile] * transmitidto2[i2tile]
                        T2T[i2tile] = pow( summation, 0.25 )
                    if( flowcontrol.diagnostics == "INSPECTHEATING"):
                        InspectHeat2ByID(T2Told, muAidto2, DeltaT4idto2, 
                                       transmitidto2 )

                if( verbose == "ON" ): 
	              print("   Begin heating star 2 by the outer disk.\n")

                for i2tile in range(1, Star2.Ntiles):
                    T2Told[i2tile] = T2T[i2tile]
                for i2tile in range(1, Star2.Ntiles):
                    if( verbose == "ON"):
                        k = i2tile / 100
                        k = i2tile - 100 * k
                        if( k == 0 ):
                            print("      heating star 2 tile number %5ld\n", i2tile)
                    DeltaT4Dto2[i2tile] = 0.0
                    end.x = T2x[i2tile]
                    end.y = T2y[i2tile]
                    end.z = T2z[i2tile]
                    for iDisktile in range(1, wholediskpars.Ntiles):
                        start.x = TDiskx[iDisktile]
                        start.y = TDisky[iDisktile]
                        start.z = TDiskz[iDisktile]
                        vectord.x = end.x - start.x
                        vectord.y = end.y - start.y
                        vectord.z = end.z - start.z
                        muEprime = CartDotProd( vectord, TDisknormCart[iDisktile] )
                        if( muEprime > 0.0 ):
                            muAprime = CartDotProd( vectord, T2normCart[i2tile] )
                            if( muAprime < 0.0 ):
                                dsquare =   vectord.x * vectord.x + vectord.y * vectord.y + vectord.z * vectord.z
                                d = math.sqrt( dsquare )
                                muE = muEprime / d
                                muA = muAprime / d
                                delT4disk = -( TDiskT4[iDisktile] / math.pi) * (muE * muA / dsquare) * TDiskdS[iDisktile]
                                delT4disk *= Transmission( start, end )
                                DeltaT4Dto2[i2tile] += delT4disk 
                            summation = pow( T2T[i2tile], 4.0) + Star2.albedo * DeltaT4Dto2[i2tile]
                            T2T[i2tile] = pow( sum, 0.25 )
                        if( flowcontrol.diagnostics == "INSPECTHEATING"):
                            InspectHeat2ByDisk( T2Told, DeltaT4Dto2, T2T)

                    if (flowcontrol.adc == "ON"):
                        if( verbose == "ON" ): 
                            print("   Begin heating star 2 by the ADC.\n")
                        for i2tile in range(1, Star2.Ntiles):
                            T2Told[i2tile] = T2T[i2tile]

                        start.x = 0.0
                        start.y = adcpars.height
                        start.z = systemparams.a
                        for i2tile in range(1, Star2.Ntiles):
                            end.x = T2x[i2tile]
                            end.y = T2y[i2tile]
                            end.z = T2z[i2tile]
                            delta.x = end.x - start.x
                            delta.y = end.y - start.y
                            delta.z = end.z - start.z
                            d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                            direction.x = delta.x / d
                            direction.y = delta.y / d
                            direction.z = delta.z / d
                            muAADCto2[i2tile] = CartDotProd( direction, T2normCart[i2tile] )
                            if( muAADCto2[i2tile] < 0.0 ):
                                DeltaT4ADCto2[i2tile] = math.abs( ADCTotFlux( d ) 
                                             * muAADCto2[i2tile] 
                                             / SIGMA )
                                transmitADCto2[i2tile] = Transmission( start, end )
                            else:
                                DeltaT4ADCto2[i2tile] = 0.0
                                transmitADCto2[i2tile] = 0.0
                                summation = pow( T2T[i2tile], 4.0) + Star2.albedo * DeltaT4ADCto2[i2tile] * transmitADCto2[i2tile]
                            T2T[i2tile] = pow( sum, 0.25 )
                            if( flowcontrol.diagnostics == "INSPECTHEATING"):
                                InspectHeat2ByADC( "TOP", T2Told, muAADCto2, 
                                  DeltaT4ADCto2,  transmitADCto2 )
                        start.x = 0.0
                        start.y = -adcpars.height
                        start.z = systemparams.a
                        for i2tile in range(1, Star2.Ntiles):
                            end.x = T2x[i2tile]
                            end.y = T2y[i2tile]
                            end.z = T2z[i2tile]
                            delta.x = end.x - start.x
                            delta.y = end.y - start.y
                            delta.z = end.z - start.z
                            d =  math.sqrt( delta.x*delta.x + delta.y*delta.y + delta.z*delta.z )
                            direction.x = delta.x / d
                            direction.y = delta.y / d
                            direction.z = delta.z / d
                            muAADCto2[i2tile] = CartDotProd( direction, T2normCart[i2tile] )
                            if( muAADCto2[i2tile] < 0.0 ):
                                DeltaT4ADCto2[i2tile] = math.abs( ADCTotFlux( d ) 
                                             * muAADCto2[i2tile] 
                                             / SIGMA )
                                transmitADCto2[i2tile] = Transmission( start, end )
                            else:
                                DeltaT4ADCto2[i2tile] = 0.0
                                transmitADCto2[i2tile] = 0.0
                            summation = pow( T2T[i2tile], 4.0)+ Star2.albedo * DeltaT4ADCto2[i2tile] * transmitADCto2[i2tile]
                        T2T[i2tile] = pow( sum, 0.25 )

                    if( flowcontrol.diagnostics == "INSPECTHEATING"):
                        InspectHeat2ByADC( "BOTTOM", T2Told, muAADCto2, 
                                  DeltaT4ADCto2,  transmitADCto2 )

                for band in range(1, orbitparams.nbands):
                    if( orbitparams.filter[band] == "SQUARE"):
                        for i2tile in range(1, Star2.Ntiles):
                            T2I[band][i2tile] = BBSquareIntensity( T2T[i2tile], 
                                                 orbitparams.minlambda[band], 
                                                 orbitparams.maxlambda[band])
                    else:
                        for i2tile in range(1, Star2.Ntiles):
                            if( (T2T[i2tile] > IperpT[maxIperpTindex]) or (T2T[i2tile] < IperpT[0]) ):
                                T2I[band][i2tile] = BBFilterIntensity( T2T[i2tile], 
                                                     orbitparams.filter[band])
                            else:
                                T2I[band][i2tile] = GetIperp( T2T[i2tile], T2logg[i2tile], 
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
    deltaray.x = steplength * direction.x
    deltaray.y = steplength * direction.y
    deltaray.z = steplength * direction.z

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
        ix = 1.5 + ( ray.x - Grid.xmin ) / XYGrid.deltax
        iz = 1.5 + ( ray.z - Grid.zmin ) / XYGrid.deltaz
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

    delta.x = end.x - start.x
    delta.y = end.y - start.y
    delta.z = end.z - start.z
    length = sqrt( delta.x * delta.x + delta.y * delta.y + delta.z * delta.z )
    direction.x = delta.x / length
    direction.y = delta.y / length
    direction.z = delta.z / length
    nsteps = length / XYGrid.deltal
    if( nsteps <= 3 ):
        transmit = 1.0
        return( transmit )
    stepsize = length / nsteps
    deltaray.x = stepsize * direction.x
    deltaray.y = stepsize * direction.y
    deltaray.z = stepsize * direction.z

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
    if( verbose == "ON"):
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
            if( T2x[itile] < XYGrid.xmin ): 
                XYGrid.xmin = T2x[itile]
            if( T2x[itile] > XYGrid.xmax ):
                XYGrid.xmax = T2x[itile]
            if( T2y[itile] < XYGrid.ymin ):
                XYGrid.ymin = T2y[itile]
            if( T2y[itile] > XYGrid.ymax ):
                XYGrid.ymax = T2y[itile]
            if( T2z[itile] < XYGrid.zmin ):
                XYGrid.zmin = T2z[itile]
            if( T2z[itile] > XYGrid.zmax ):
                XYGrid.zmax = T2z[itile]
    if( flowcontrol.disk == "ON" ):
        for itile in range(1, wholediskpars.Ntiles):
            if( TDiskx[itile] < XYGrid.xmin ):
                XYGrid.xmin = TDiskx[itile]
            if( TDiskx[itile] > XYGrid.xmax ):
                XYGrid.xmax = TDiskx[itile]
            if( TDisky[itile] < XYGrid.ymin ):
                XYGrid.ymin = TDisky[itile]
            if( TDisky[itile] > XYGrid.ymax ):
                XYGrid.ymax = TDisky[itile]
            if( TDiskz[itile] < XYGrid.zmin ):
                XYGrid.zmin = TDiskz[itile]
            if( TDiskz[itile] > XYGrid.zmax ):
                XYGrid.zmax = TDiskz[itile]

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
	               XYGrid.Topy[ix][iz] = Star2TopY( x, z )
	               XYGrid.Bottomy[ix][iz] = -XYGrid.Topy[ix][iz]
            else:
                if( flowcontrol.disk == "ON" ):
                    rho = math.sqrt( x*x + (z - systemparams.a)*(z - systemparams.a) )
                    coszeta = (z - systemparams.a) / rho
                    zeta = math.acos( coszeta )
                    if( x < 0.0 ):
                        zeta = math.pi*2 - zeta
                    a = RhoToA( rho, zeta)
                    Grid.Topy[ix][iz] = DiskTopH( a, zeta)
                    Grid.Bottomy[ix][iz] = DiskBottomH( a, zeta)

    if( flowcontrol.diagnostics == "INSPECTYLIMITS"):
        InspectYlimits()
        sys.exit("Quit in MakeYlimits after INSPECTYLIMITS.")

    return