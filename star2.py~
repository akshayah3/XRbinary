# -*- coding: utf-8 -*-
"""
All functions concerned with the secondary star (the lobe-filling star or 
"star 2") of the binary system are in this file.
"""

G  = 6.6742e-8
SIGMA = 5.6704e-5
import sys
import math
from .parmeter import systemparams
from .parmeter import XYGrid, globalvar

class star2:
    targetNtiles = None
    meanT = None
    albedo = None
    L = None
    volume = None
    meang = None
    meanr = None
    logg = None
    beta = None
    frontradius = None
    backradius = None
    sideradius = None
    poleradius = None
    """
    Parameters
    ==========
    """
    
    def __init__(self, targetNtiles, meanT, albedo, L):
        self.targetNtiles = targetNtiles
        self.albedo = albedo
        self.L = L
        self.meanT = meanT

    def MakeStar2Tiles(self):
        """
        This function distributes the tiles over the surface of the secondary 
        star and calculates various properties of the tiles.
        """

    def V(self, r, theta, phi):
        """
        This function calculates the potential of the zero velocity surface at the 
        point (r, theta, phi).
        """
        l = math.sin(theta) * math.cos(phi);
        n = math.cos(theta);
        r1 = math.sqrt( systemparams.a * systemparams.a  +  r * r  - 2.0 * systemparams.a  * r * n )
        r2 = r
        if (r1 <= 0.0) or (r2 <= 0.0):
            print("Failure in function V.")
        rho = math.sqrt( (r*l)*(r*l) + (r*n - systemparams.zcm)*(r*n - systemparams.zcm) )
        potential =   systemparams.a / r1 + systemparams.q * systemparams.a / r2 + 0.5 * (1 + systemparams.q) * (rho / systemparams.a) * (rho / systemparams.a)
        potential = -(G * systemparams.M1 / systemparams.a) * potential
        return potential

    def gradV(self, r, theta, phi):
        """
        This function calculates the gradient of the zero velocity potential
        at the point (r, theta, phi).  The gradient is returned as a vector in 
        spherical polar coordinates.
        """  
        epsilon    = 1.0e-7
        deltar     = epsilon * systemparams.a
        deltatheta = epsilon * math.pi
        deltaphi   = epsilon * (math.pi*2)

        gradVr = (self.V( r + deltar, theta, phi) 
            - self.V( r - deltar, theta, phi) ) / (2.0 * deltar)
        if((theta <= deltatheta) or (theta >= (math.pi - deltatheta))):
            gradVt = 0.0
            gradVp = 0.0
        else:
            dVdtheta = (self.V( r, theta + deltatheta, phi) 
                   - self.V( r, theta - deltatheta, phi) ) / (2.0 * deltatheta)
            gradVt   = ( 1.0 / r ) * dVdtheta;
            dVdphi   = (self.V( r, theta, phi + deltaphi)
                   - self.V( r, theta, phi - deltaphi) )   / (2.0 * deltaphi)
            gradVp   = ( 1.0 / (r * math.sin(theta) ) ) * dVdphi

        return gradVr, gradVt, gradVp

    def Normal(gradVr, gradVt, gradVp):
        """
        This function calculates the normal vector by taking normal = delV / |delV|
        where delV is the gradient vector.  The vectors are assumed to be 
        expressed in spherical polar coordinates.
        """    
        denom = math.sqrt(gradVr*gradVr + gradVt*gradVt + gradVp*gradVp)
        if denom > 0.0:
            normVr = gradVr / denom
            normVt = gradVt / denom
            normVp = gradVp / denom
        else:
            normVr = 1.0
            normVt = 0.0
            normVp = 0.0
   
        return normVr, normVt, normVp

    def FindL1(self):
        """
        This function finds the position of the inner Lagrangian point.
        """
        epsilon = 1.0e-7;
        epsilonr = epsilon * systemparams.a;
        rOld = 0.5 * systemparams.a;
        for i in range(500):
            delVplus, _, _  = self.gradV( rOld + epsilonr, 0.0, 0.0)
            delVminus, _, _ = self.gradV( rOld - epsilonr, 0.0, 0.0)
            slope = (delVplus - delVminus) / (2.0 * epsilonr)
            deltar = 0.5 * (delVplus + delVminus) / slope
            r = rOld - deltar
            if (r <= 0.0) or (r >= systemparams.a):
                print("FindL1 failed.")
            if i == 500:
                print("Too many iterations in function FindL1.")
            x = math.abs( (rOld - r) / rOld )
            if ( x < epsilon ):
                break
            rOld = r
        return r

    def findR(self, Vtarget, theta, phi):
        """
        This function finds the radius at which the zero velocity potential equals 
        targetV in the direction (theta, phi).  The radius is given as a 
        fraction of a. NOTE: This function only works if the target potential 
        is less than the potential at the Roche lobe.
        """
        if( Vtarget > systemparams.VL1 ): 
            print("Attempted to find a point outside the Roche lobe in findR).")

        epsilonr = 1.0e-7 * systemparams.a
        r = 0.8 * systemparams.rL1
        for i in range(100):
            if i ==100:
                print("Too many iterations in findR.")
            Vplus  = self.V( r + epsilonr, theta, phi)
            Vminus = self.V( r - epsilonr, theta, phi)
            Vzero  = self.V( r,           theta, phi)
            slope = (Vplus - Vminus) / (2.0 * epsilonr)
            if slope == 0.0:
                slope = 100.0 * Vzero / systemparams.a
            deltar = (Vtarget - Vzero) / slope
            if( math.abs(deltar) > (0.05 * r) ): 
                deltar = (0.05 * r) * (deltar / math.abs(deltar))
            r = r + deltar 
            if( r >= systemparams.rL1 ):
                r = systemparams.rL1
            x = math.abs( (Vtarget - self.V(r, theta, phi)) / Vtarget )
            if( x <= 1.0e-6 ):
                break
            return r

    def TileArea(self, r, theta, normVr, dtheta, dphi):
        """
        This function calculates the surface area of the tiles.  A tile has 
        sides with lengths dtheta and dphi.  The area of the tile is set equal 
        to the area of that part of the curved equipotential surface represented
        by the flat tile.  r is the radius at the center of the segment. The 
        vector "normals" must be in spherical polar coords.
        """

        dS =   2.0 * r * r * dphi * math.sin( theta ) * math.sin( dtheta / 2.0 )
        dS = dS / normVr

        return dS

    def Star2TopY(self, x, z):
        """
        This function returns the value of y at the surface of star 2 for (x,z)
        """
        converge = 1.0e-6
        Vtarget = systemparams.VL1

        if( z >= systemparams.rL1 ):
            topy = 0.0
            return topy
        r = math.sqrt( x*x + z*z )
        if r > 0.0:
            costheta = z / r
            theta = math.acos( costheta )
            phi = 0.0
        if x < 0.0:
            phi = math.pi
            Vzero = self.V(r, theta, phi)
        if Vzero > Vtarget:
            topy = 0.0
            return topy

        epsilony = 1.0e-7 * systemparams.a
        y = 0.5 * systemparams.rL1
        for i in range(150):
            if i == 149:
	           print( x, z, y)
            print("Too many iterations in Star2Topy.")

        r = math.sqrt( x*x + y*y + z*z)
        costheta = z / r
        theta = math.acos( costheta )
        cosphi = x / math.sqrt( x*x + y*y )
        phi = math.acos( cosphi )
        Vzero  = self.V( r, theta, phi)
        change = (Vzero - Vtarget) / Vtarget
        if( math.abs( change ) <= converge ):
            break

        yplus = y + epsilony
        r = math.sqrt(x*x + yplus*yplus + z*z)
        costheta = z / r
        theta = math.acos(costheta);
        cosphi = x / math.sqrt(x*x + yplus*yplus)
        phi = math.acos( cosphi )
        Vplus  = self.V( r, theta, phi)

        yminus = y - epsilony
        if( yminus < 0.0 ):
	       yminus = 0.0
        r = math.sqrt( x*x + yminus*yminus + z*z)
        costheta = z / r
        theta = math.acos(costheta)
        cosphi = x / math.sqrt( x*x + yminus*yminus )
        phi = math.acos(cosphi)
        Vminus = self.V(r, theta, phi)

        slope = (Vplus - Vminus) / (yplus - yminus)
        if(slope == 0.0):
            slope = 100.0 * Vzero / systemparams.a
        deltay = (Vtarget - Vzero) / slope
        if( math.abs(deltay) > (0.05 * y) ): 
            deltay = (0.05 * y) * ( deltay / math.abs( deltay ) )
        y = y + deltay 
        if( y > (1.01 * XYGrid.ymax) ):
            y = 1.01 * XYGrid.ymax
        if( y < 0.0 ):
	       y = 0.0
        topy = y 
    
        return topy

    def GetGDbeta():
        """
        This function extracts the appropriate value of the gravity
        darkening beta parameter from the fourbeta[] array, where
        gravity darkening law is:
        Teff = <Teff> * ( g / <g> )^beta
        """
        if( globalvar.maxGDindex == 1):
            beta = 0.25 * globalvar.fourbeta[1]
            return( beta )
        if( star2.meanT <= globalvar.GDT[0] ):
             beta = 0.25 * globalvar.fourbeta[1]
             return( beta )
        if( star2.meanT >= globalvar.GDT[globalvar.maxGDindex] ):
            beta = 0.25 * globalvar.fourbeta[globalvar.maxGDindex]
            return( beta )      
        for i in range(0, globalvar.maxGDindex):
            if( (star2.meanT >= globalvar.GDT[i]) and (star2.meanT < globalvar.GDT[i+1]) ):
                ilow = i
                ihigh = i + 1
                break

        if( globalvar.GDT[ihigh] == globalvar.GDT[ilow] ):
            beta = 0.25 * globalvar.fourbeta[ilow]
            return( beta )
        slope = ( globalvar.fourbeta[ihigh] - globalvar.fourbeta[ilow] ) / ( globalvar.GDT[ihigh] - globalvar.GDT[ilow] )
        beta = globalvar.fourbeta[ilow] + slope * ( star2.meanT - globalvar.GDT[ilow] )
        beta = 0.25 * beta;

        return( beta )


    def ClaretHmu( T, logg, Filter,mu):
        """
        This is the Claret (2000, AA, 363, 1081) limb-darkening function.
        It does not interpolate in logg and T but instead just goes to 
        the nearest grid point.
        It does interpolate in mu.
        """
        findex = -1
        for i in range(0, globalvar.maxLDfilterindex):
            if( globalvar.LDfilterName[i] == Filter):
                findex = i
                break
        if( findex == -1 ):
            sys.exit("Unknown filter name in ClaretHmu.")

        if( T <= globalvar.LDT[0] ):
            Tindex = 0
        elif( T >= globalvar.LDT[globalvar.maxLDTindex] ):
            Tindex = globalvar.maxLDTindex
        else:
            deltaT = globalvar.LDT[1] - globalvar.LDT[0]
            Tindex = 0.5 + ( T - globalvar.LDT[0] ) / deltaT
        if( logg <= globalvar.LDlogg[0] ):
            gindex = 0
        elif( logg >= globalvar.LDlogg[globalvar.maxLDgindex] ):
            gindex = globalvar.maxLDgindex
        else:
            deltag =  globalvar.LDlogg[1] - globalvar.LDlogg[0]
            gindex = 0.5 + ( logg - globalvar.LDlogg[0] ) / deltag

        a1 = globalvar.LDtable[gindex][Tindex][findex][1]
        a2 = globalvar.LDtable[gindex][Tindex][findex][2]
        a3 = globalvar.LDtable[gindex][Tindex][findex][3]
        a4 = globalvar.LDtable[gindex][Tindex][findex][4]
        if( mu < 0.0 ):
            sys.exit("mu less than zero in ClaretHmu.")
        sqrtmu = math.sqrt( mu )
        h = 1.0 - a1 * (1.0 - sqrtmu)    - a2 * (1.0 - mu) - a3 * (1.0 - mu*sqrtmu) - a4 * (1.0 - mu*mu)
        return ( h )


    def GetIperp( T, logg, Filter ):
        """
        This function returns the specific intensity at zero zenith
        angle by linearly interpolated Iperptable[][][] in temperature
        and log(g).
        """
        findex = -1
        for i in range(0, globalvar.maxIperpfilterindex):
            if( globalvar.IperpfilterName[i] == Filter):
                findex = i
                break
        if( findex == -1 ):
            sys.exit("Unknown filter name in GetIperp.")

        if( T <= globalvar.IperpT[0] ):
            sys.exit("T out of range (too low) in GetIperp.")
        elif( T >= globalvar.IperpT[globalvar.maxIperpTindex] ):
            sys.exit("T out of range (too high) in GetIperp.")
        else:
            deltaT = globalvar.IperpT[1] - globalvar.IperpT[0]
            Tindex1 = ( T - globalvar.IperpT[0] + 0.1 ) / deltaT
            Tindex2 = Tindex1 + 1
            weightT2 = ( T - globalvar.IperpT[Tindex1] ) / deltaT
            weightT1 = 1.0 - weightT2
        if( logg <= globalvar.Iperplogg[0] ):
            gindex1 = 0
            gindex2 = 0
            weightg1 = 1.0
            weightg2 = 0.0
        elif( logg >= globalvar.Iperplogg[globalvar.maxIperpgindex] ):
            gindex1 = globalvar.maxIperpgindex
            gindex2 = globalvar.maxIperpgindex
            weightg1 = 0.0
            weightg2 = 1.0
        else:
            deltag =  globalvar.Iperplogg[1] - globalvar.Iperplogg[0]
            gindex1 = ( logg - globalvar.Iperplogg[0] + 0.1 ) / deltag
            gindex2 = gindex1 + 1
            weightg2 = ( logg - globalvar.Iperplogg[gindex1] ) / deltag
            weightg1 = 1.0 - weightg2
        intensity =  weightg1 * weightT1 * globalvar.Iperptable[gindex1][Tindex1][findex] + weightg1 * weightT2 * globalvar.Iperptable[gindex1][Tindex2][findex]+ weightg2 * weightT1 * globalvar.Iperptable[gindex2][Tindex1][findex]+ weightg2 * weightT2 * globalvar.Iperptable[gindex2][Tindex2][findex]

        return( intensity )


    def Star2L():
        """
        Calculate the luminosity of star 2 by adding up the fluxes
        from all the tiles.  This function should not be used until
        after heating by irradiation has been calculated.
        """
        luminosity = 0.0
        for itile in range(1, star2.Ntiles):
            luminosity += SIGMA * pow( globalvar.T2T[itile], 4.0) * globalvar.T2dS[itile]

        return( luminosity )