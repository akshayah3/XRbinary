# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 16:53:23 2016

"""

"""
Utility functions needed for the project
"""

import math
import sys
from parmeter import globalvar, CylVector, CartVector

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
    
    if( (Zeta1 < globalvar.deltaBBzeta) or (Zeta2 < globalvar.deltaBBzeta) ):
        print("lambda x T too high in BBSquareIntensity().")

    if( Zeta1 >= globalvar.BBzetamax ):
        meanI = 0.0
        return meanI
    else:
        logZeta1 = math.log( Zeta1 )
        nlow = Zeta1 / globalvar.deltaBBzeta
        ZetaLow = nlow * globalvar.deltaBBzeta
        logZetaLow = math.log( ZetaLow )
        nhigh = nlow + 1
        logZetaHigh = math.log( ZetaLow + globalvar.deltaBBzeta )
        logZBBzetaLow  = math.log( globalvar.ZBBzeta[nlow]  )
        logZBBzetaHigh = math.log( globalvar.ZBBzeta[nhigh] )
        slope = (logZBBzetaHigh - logZBBzetaLow) / (logZetaHigh - logZetaLow)
        logZBBzeta1 = logZBBzetaLow + (logZeta1 - logZetaLow) * slope
        ZBBzeta1 = math.exp( logZBBzeta1 )

    if( Zeta2 >= globalvar.BBzetamax ):
      ZBBzeta2 = globalvar.ZBBzeta[globalvar.maxBBzetaindex]

    else:
        logZeta2 = math.log( Zeta2 )
        nlow = Zeta2 / globalvar.deltaBBzeta
        ZetaLow = nlow * globalvar.deltaBBzeta
        logZetaLow = math.log( ZetaLow )
        nhigh = nlow + 1
        logZetaHigh = math.log( ZetaLow + globalvar.deltaBBzeta )
        logZBBzetaLow  = math.log( globalvar.ZBBzeta[nlow]  )
        logZBBzetaHigh = math.log( globalvar.ZBBzeta[nhigh] )
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
    for i in range(globalvar.maxIBBfilterindex):
        if globalvar.IBBfilterName[i] == filter:
            findex = i
            break

    if( findex == -1 ):
        print("Unknown filter name in BBFilterIntensity.");

    if( (T < globalvar.IBBT[0]) or (T > globalvar.IBBT[globalvar.maxIBBTindex]) ):
        print("T out of range in BBFilterIntensity.")
    if( T == globalvar.IBBT[globalvar.maxIBBTindex] ):
        intensity = globalvar.IBBtable[globalvar.maxIBBTindex][findex]
        return( intensity )

    minTindex = (T - globalvar.IBBTmin) / globalvar.IBBdeltaT
    maxTindex = minTindex + 1
    slope = (globalvar.IBBtable[maxTindex][findex] - globalvar.IBBtable[minTindex][findex]) / globalvar.IBBdeltaT
    intensity = slope * (T - globalvar.IBBT[minTindex]) + globalvar.IBBtable[minTindex][findex]

    return intensity

def AngleDistance( theta1, phi1, theta2, phi2):
    """
    Calculates the angular distance between two directions, where
    the directions are given by their coordinates (theta,phi) in
    the spherical polar coordinate system.  All angles in radians.
    """
    cosa = math.cos(theta1) * math.cos(theta2) + math.sin(theta1) * math.sin(theta2) * math.cos(phi2 - phi1)
    a = math.acos( cosa )

    return( a )

def CartDotProd( A, B):
    """
    Calculates the dot product of two vectors, both of which
    are in Cartesian coordinates.
    """
    prod = A.x * B.x + A.y * B.y + A.z * B.z

    return( prod )

def Cart2Sphere( Acart, theta, phi ):
    """
    This converts a vector from Cartesian coordinates to 
    spherical polar coordinates.  More specifically, it calculates
    the components of the vector in spherical polar coordinates
    given its components in Cartesian coordinates.
    """
    sint = math.sin( theta )
    cost = math.cos( theta )
    sinp = math.sin( phi )
    cosp = math.cos( phi )

    m11 =   sint * cosp
    m21 =   cost * cosp
    m31 = - sinp

    m12 =   sint * sinp
    m22 =   cost * sinp
    m32 =   cosp

    m13 =   cost
    m23 = - sint
    m33 =   0.0

    Asphere = CylVector()
    Asphere.r     = m11 * Acart.x + m12 * Acart.y + m13 * Acart.z
    Asphere.theta = m21 * Acart.x + m22 * Acart.y + m23 * Acart.z
    Asphere.phi   = m31 * Acart.x + m32 * Acart.y + m33 * Acart.z

    return( Asphere )

def Sphere2Cart( Asphere, theta, phi ):
    """
    This converts a vector from spherical polar coordinates to 
    Cartesian coordinates.  More specifically, it calculates
    the components of the vector in Cartesian coordinates
    given its components in Spherical polar coordinates.
    """
    sint = math.sin( theta )
    cost = math.cos( theta )
    sinp = math.sin( phi )
    cosp = math.cos( phi )

    m11 =   sint * cosp
    m21 =   sint * sinp
    m31 =   cost

    m12 =   cost * cosp
    m22 =   cost * sinp
    m32 = - sint

    m13 = - sinp
    m23 =   cosp
    m33 =   0.0

    Acart = CartVector()
    Acart.x = m11 * Asphere.r + m12 * Asphere.theta + m13 * Asphere.phi
    Acart.y = m21 * Asphere.r + m22 * Asphere.theta + m23 * Asphere.phi
    Acart.z = m31 * Asphere.r + m32 * Asphere.theta + m33 * Asphere.phi

    return( Acart )

def Cyl2Cart( Asphere, zeta ):
    """
    This converts a vector from cylindrical coordinates to 
    Cartesian coordinates.  More specifically, it calculates
    the components of the vector in Cartesian coordinates
    given its components in cylindrical coordinates.
    """
    sinzeta = math.sin( zeta )
    coszeta = math.cos( zeta )

    m11 =   sinzeta
    m21 =   0.0
    m31 =   coszeta

    m12 =   coszeta
    m22 =   0.0
    m32 = - sinzeta

    m13 =   0.0
    m23 =   1
    m33 =   0.0

    Acart = CartVector()
    Acart.x = m11 * Asphere.rho + m12 * Asphere.zeta + m13 * Asphere.h
    Acart.y = m21 * Asphere.rho + m22 * Asphere.zeta + m23 * Asphere.h
    Acart.z = m31 * Asphere.rho + m32 * Asphere.zeta + m33 * Asphere.h

    return( Acart )

def Planck( mode, temperature, lambdanu ):
    """
    Returns the specific intensity emitted by the surface of a black body:
 
      Fnu     = ( 2*h*nu^3 / c^2)      / ( exp[ h*nu / k*T ] - 1 )

       Flambda = ( 2*h*c^2 / lambda^5 ) / ( exp[ h*c / lambda*k*T ] -1 )

     Input data:
      mode         "NU"     returns Fnu     in units of erg/sec/cm^2/Hz
                   "LAMBDA" returns Flambda in units of erg/sec/cm^2/cm
      temperature  in degrees Kelvin
      lambdanu     wavelength in cm if mode=LAMBDA
                   frequency  in Hz if mode=NU
 
    Note that the Planck function is normalized such that 
      (integral over frequency) = ( sigma / pi ) * T**4
    Thus, it is the monochromatic specific intensity per unit area.
    """
    h = 6.62608e-27
    c=  2.99792e+10 
    k = 1.38066e-16 
    if( temperature <= 1.0 ):
        sys.exit("Temperature out of range in function Planck.")
    if( mode == "NU" ):
        nu = lambdanu
        x = ( 2.0 * h * nu * nu * nu ) / ( c * c )
        y = ( h * nu ) / ( k * temperature )
        Fnu = x / ( math.exp( y ) - 1.0 )
        return( Fnu )
    elif( mode == "LAMBDA"):
        Lambda = lambdanu
        x = ( 2.0 * h * c * c ) / pow( Lambda, 5.0 )
        y = ( h * c ) / ( Lambda * k * temperature )
        Flambda = x / ( math.exp( y ) - 1.0 )
        return( Flambda )
    else:
        sys.exit("Unrecognized mode in function Planck().")
      
    return(-1.0)
