# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:05:51 2016

"""

import math
from .star1 import Star1
from .star2 import Star2
from .parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, XYGrid, dataparams, ReadInput

sigma = 5.6704e-5


class maindisk:
    amin = None
    amax = None
    da = None
    Hmax = None
    Hpow = None
    L = None
    Tpow = None
    Tamax = None
    Tamin = None
    
    def __init__(self, amin, amax, da, Hmax, Hpow, L, Tpow, Tamax, Tamin):
        self.amin = amin
        self.amax = amax
        self.da = da
        self.Hmax = Hmax
        self.Hpow = Hpow
        self.L = L
        self.Tpow = Tpow
        self.Tamax = Tamax
        self.Tamin = Tamin


    def DiskTopT(self, a, zeta):
        """
        This function returns the temperature of the top of the disk. The 
        temperature of the rim overrides the temperature of the disk if its 
        temperature is greater than the temperature of the underlying main disk.
        """
        
        if( (a > self.amax) or (a < self.amin) ):
            temperature = 0.0
            return temperature
        temperature = self.MainDiskT( a, zeta )

        if( flowcontrol.disktorus == "ON"):
            if( (a <= (disktorusparams.azero + 0.5 * disktorusparams.awidth))
               & (a >= (disktorusparams.azero - 0.5 * disktorusparams.awidth)) ):
	          torusT = self.DiskTorusT( a, zeta)
            if( torusT > temperature ):
                temperature = torusT
            return temperature
        if( flowcontrol.diskrim == "ON"):
            if( a >= (self.amax - diskrimpars.awidth) ):
                rimT = self.DiskRimT( a, zeta)
            if( rimT > temperature ):
                temperature = rimT
            return temperature
        return temperature

    def MainDiskT(self, a, zeta):
        """
        This function returns the _mean_ temperature of tiles on the main disk.
        The function returns temperature = 0 if a is not in the range
        (amin + 0.5*da < a < (amax - 0.5*da)
        The initial temperature distribution is axisymmetric and can be that of
        either a steady-state viscous disk or a power law of the form
        T = Tzero * r^Tpow
        The temperatures are normalized to give the correct main disk luminosity,
        which is an input parameter. The function then checks whether (a,zeta) 
        falls within a disk spot and then if so, calculates the revised T from 
        T = T * diskspot.spotToverT 
        Note that if there is more than one spot and the spots overlap,the 
        fraction by which the spot temperature is reduced are multiplied in the overlap region.
        """
        if( self.TConstant <= 0.0 ):
            if( self.Ttype == "VISCOUS"):
                ymin = 1.0 / 3.0
                ymax = 1.0 - (2.0 /3.0) * math.sqrt(self.amin / self.amax)
                x = (1.0 / self.amin) * ymin - (1.0 / self.amax) * ymax
                self.TConstant = self.L / ( 2.0 * math.pi * sigma * x)
            elif( self.Ttype == "POWERLAW"):
                beta = 4.0 * self.Tpow + 2.0
                if( self.Tpow == -0.5 ):
                    x = self.L / (4.0 * math.pi * sigma )
                    self.TConstant = x / math.log(self.amax / self.amin)
                else:
                    x = (beta * self.L) / (4.0 * math.pi * sigma )
                    self.TConstant = x / ( math.pow(self.amax, beta) - math.pow(self.amin, beta))
                    if(self.TConstant < 0.0):
                        self.TConstant = -self.TConstant
            else:
                sys.exit("Unrecognized temperature distribution in function MainDiskT")

        if( (a < (self.amin + 0.5*self.da)) or (a > (self.amax - 0.5*self.da)) ):
            temperature = 0.0
        else:
            if( self.Ttype == "VISCOUS"):
                r1 = a - self.da / 2.0
                r2 = a + self.da / 2.0
                ymin = 1.0 - (2.0 / 3.0) * math.sqrt( self.amin / r1)
                ymax = 1.0 - (2.0 / 3.0) * math.sqrt( self.amin / r2)
                x = (1.0 / r1) * ymin - (1.0 /r2 ) * ymax 
                T4 = (self.TConstant * x) / ( r2*r2 - r1 * r1)
                temperature = math.pow( T4, 0.25)
            elif( self.Ttype == "POWERLAW"):
                r1 = a - self.da / 2.0
                r2 = a + self.da / 2.0
                if( self.Tpow == -0.5 ):
                    T4 = 2.0 * self.TConstant
                    T4 = T4 * math.log( r2 / r1 ) / (r2*r2 - r1*r1)
                else:
                    beta = 4.0 * self.Tpow + 2.0
                    T4 = 2.0 * self.TConstant / beta
                    T4 = T4 * ( math.pow( r2, beta) - math.pow( r1, beta)) / (r2*r2 - r1*r1)
                    temperature = math.pow ( T4, 0.25)
            if( flowcontrol.diskspots == "ON"):
                for i in range(1, diskspotpars.nspots):
                    if( (zeta >= diskspotpars.zetamin[i]) 
                             & (zeta <= diskspotpars.zetamax[i]) ):
                        if( (a >= diskspotpars.amin[i]) 
                             & (a <= diskspotpars.amax[i]) ):
                            temperature *= diskspotpars.spotToverT[i]
        return temperature

    def DiskTorusT(self, a, zeta):
        """
        The temperature of the torus varies with zeta but not a. There are 
        currently two possibilities for the zeta dependence
        1) SINUSOID
        T = 0.5 * (Tmax + Tmin)
        + 0.5* (Tmax - Tmin) * cos( zeta - zetaTmax );

        2) POINT
        The disk torus height and temperature is defined by a set of points, 
        one point per DISKTORUSPARS= line in the parameter file:
        DISKTORUSPARS=  POINT   Zeta1   H1   T1
        DISKTORUSPARS=  POINT   Zeta2   H2   T2    
            .        .       .     .    .
            .        .       .     .    .
        The temperatures are linearly interpolated between the specified 
        points.

        The Zetas must be in increasing order and disktorus.PointZeta[1]
        must be 0 degrees (this avoids messy computer code).
        """
        if( zeta > 2*math.pi ):
            sys.exit("zeta greater than TWOPI in DiskTorusT.")
        if( zeta < 0.0 ):
            sys.exit("zeta less than zero in DiskTorusT.")

        if( a < (disktorusparams.azero - 0.5 * disktorusparams.awidth) ):
            temperature = 0.0
            return temperature
        if( a > (disktorusparams.azero + 0.5 * disktorusparams.awidth) ):
            temperature = 0.0
            return temperature

        if( disktorusparams.type == "SINUSOID"):
            temperature =    0.5 * ( disktorusparams.Tmax + disktorusparams.Tmin ) + 0.5 * ( disktorusparams.Tmax - disktorusparams.Tmin ) * math.cos( zeta - disktorusparams.ZetaTmax )
        elif( disktorusparams.type == "POINT" ):
            if( disktorusparams.points == 1 ):
                temperature = disktorusparams.PointT[1]
            else:
                for i in range(diskrimpars.points):
                    zetalow = disktorusparams.PointZeta[i]
                    Tlow = disktorusparams.PointT[i]
            if( i < disktorusparams.points ):
                zetahigh = disktorusparams.PointZeta[i+1]
                Thigh = disktorusparams.PointT[i+1]
            else:
                zetahigh = 2*math.pi
                Thigh = disktorusparams.PointT[1]
            if( (zeta >= zetalow) and (zeta < zetahigh) ):
                slope = (Thigh - Tlow) / (zetahigh - zetalow)
                temperature = Tlow + slope * (zeta - zetalow)
                break
        else:
            sys.exit("Unrecognized disk torus type in DiskTorusT.")

        return temperature

    def DiskTorusT(self, a, zeta):
        """
        The temperature of the torus varies with zeta but not a.
        There are currently two possibilities for the zeta dependence

        1) SINUSOID
           T = 0.5 * (Tmax + Tmin)
               + 0.5* (Tmax - Tmin) * cos( zeta - zetaTmax );

        2) POINT
           The disk torus height and temperature is defined by a set of
           points, one point per DISKTORUSPARS= line in the parameter file:
           DISKTORUSPARS=  POINT   Zeta1   H1   T1
           DISKTORUSPARS=  POINT   Zeta2   H2   T2    
            .        .       .     .    .
            .        .       .     .    .
           The temperatures are linearly interpolated between the specified 
           points.

        The Zetas must be in increasing order and disktorus.PointZeta[1]
        must be 0 degrees (this avoids messy computer code).
        """
        if( zeta > 2*math.pi ):
            sys.exit("zeta greater than TWOPI in DiskTorusT.")
        if( zeta < 0.0 ):
            sys.exit("zeta less than zero in DiskTorusT.")

        if( a < (disktorusparams.azero - 0.5 * disktorusparams.awidth) ):
            temperature = 0.0
            return temperature
        if( a > (disktorusparams.azero + 0.5 * disktorusparams.awidth) ):
            temperature = 0.0
            return temperature
        if( disktorusparams.type == "SINUSOID"):
            temperature =    0.5 * ( disktorusparams.Tmax + disktorusparams.Tmin ) + 0.5 * ( disktorusparams.Tmax - disktorusparams.Tmin ) * math.cos( zeta - disktorusparams.ZetaTmax )
        elif( disktorusparams.type == "POINT"):
            if( disktorusparams.points == 1 ):
                temperature = disktorusparams.PointT[1]
            else:
	          for i in range(disktorusparams.points):
                    zetalow = disktorusparams.PointZeta[i]
                    Tlow = disktorusparams.PointT[i]
                    if( i < disktorusparams.points ):
                        zetahigh = disktorusparams.PointZeta[i+1]
                        Thigh = disktorusparams.PointT[i+1]
                    else:
                        zetahigh = 2*math.pi
                        Thigh = disktorusparams.PointT[1]
                    if( (zeta >= zetalow) and (zeta < zetahigh) ):
                        slope = (Thigh - Tlow) / (zetahigh - zetalow)
                        temperature = Tlow + slope * (zeta - zetalow)
                        break
        else:
            sys.quit("Unrecognized disk torus type in DiskTorusT.")

        return temperature

    def DiskEdgeT(self, zeta):
        """
        
        """
        temperature = diskedgepars.T;
        if( diskedgepars.Tspot > diskedgepars.T ):
            z1 = diskedgepars.ZetaMid - diskedgepars.ZetaWidth / 2.0
            z2 = diskedgepars.ZetaMid + diskedgepars.ZetaWidth / 2.0
            if( z1 < 0.0 ):
                if( zeta >= (z1 + 2*math.pi) ):
                    temperature = diskedgepars.Tspot
                if( zeta <= z2 ):
                    temperature = diskedgepars.Tspot
            elif( z2 > 2*math.pi ):
                if( zeta < (z2 - 2*math.pi) ):
	              temperature = diskedgepars.Tspot
                if( zeta >= z1 ):
 	              temperature = diskedgepars.Tspot
            elif( (zeta >= z1) and (zeta <= z2) ):
	          temperature = diskedgepars.Tspot
        return temperature

    def DiskL(self):
        """
        Calculate the luminosity of the disk by adding up the fluxes from all 
        its tiles.This function should not be used until after heating by 
        irradiation has been calculated.  Note that the luminosity of the disk 
        will be different from maindiskL if there are spots, edges, torii, etc.
        on the disk.Also find the maximum and minimum temperatures of the disk
        tiles on the top and bottom (not the edge) of the disk.
        """
        luminosity = 0.0;
        wholediskpars.TopTmax = 0.0;
        wholediskpars.TopTmin = 1.0e12;
        for i in range(wholediskpars.ntiles):
            luminosity += sigma * TDiskT4[i] * TDiskdS[i]
            if( TDiska[itile] < 0.999 * self.amax ):
                if( TDiskT[itile] > wholediskpars.TopTmax ):
                    wholediskpars.TopTmax = TDiskT[itile]
                if( TDiskT[itile] < wholediskpars.TopTmin ):
                    wholediskpars.TopTmin = TDiskT[itile]
        return luminosity

    def InnerDiskFlambda( Filter, minlambda, maxlambda):
        """
        This function returns the contribution of the inner disk to the
        observed spectrum at wavelength lambda.  Note that the wavelengths
        must be in centimeters but the returned flux is in
        ergs/sec/cm^2/Angstrom.  The returned quantity must be 
        multiplied by the geometric projection factor
                    mu = cos( theta ) 
        to get the observed quantity.
        """
        if( Filter == "SQUARE"):
            flux = ( innerdiskpars.L / (2.0 * innerdiskpars.sigmaT4) )* BBSquareIntensity( innerdiskpars.T, minlambda, maxlambda )
        else:
            flux = ( innerdiskpars.L / (2.0 * innerdiskpars.sigmaT4) )* BBFilterIntensity( innerdiskpars.T, Filter )

        return( flux )


    def ADCTotFlux(distance):
        """
        This function returns the integrated flux from ONE of
        the ADC points:
        The integrated flux is just adc.L/2.0 diluted by the
        area of the sphere around the ADC point.
        """
        totalflux = 0.5 * adcpars.L / (4*math.pi * distance * distance )
        return( totalflux )

    def InnerDiskTotFlux(distance):
        """
        This function returns the integrated flux from the inner disk.
        The flux has been integrated over wavelength and over 
        the surface of the disk.  The returned quantity must by
        multiplied by the geometric projection factor 
             mu = cos( theta ) 
        to get the irradiating flux.  The factor is TWOPI
        instead of FOURPI because of the mu factor.
        """
        totalflux = innerdiskpars.L / (2*math.pi * distance * distance )
        return( totalflux )