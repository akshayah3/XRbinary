# -*- coding: utf-8 -*-
"""
All functions concerned with the secondary star (the lobe-filling star or 
"star 2") of the binary system are in this file.
"""

import math

class star2:
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
        r1 = math.sqrt( syspars.a * syspars.a  +  r * r  - 2.0 * syspars.a  * r * n )
        r2 = r
        if (r1 <= 0.0) or (r2 <= 0.0):
            print("Failure in function V.")
        rho = math.sqrt( (r*l)*(r*l) + (r*n - syspars.zcm)*(r*n - syspars.zcm) )
        potential =   syspars.a / r1 + syspars.q * syspars.a / r2 + 0.5 * (1 + syspars.q) * (rho / syspars.a) * (rho / syspars.a)
        potential = -(G * syspars.M1 / syspars.a) * potential
        return potential

    def gradV(self, r, theta, phi):
        """
        This function calculates the gradient of the zero velocity potential
        at the point (r, theta, phi).  The gradient is returned as a vector in 
        spherical polar coordinates.
        """  
        epsilon    = 1.0e-7
        deltar     = epsilon * syspars.a
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
        epsilonr = epsilon * syspars.a;
        rOld = 0.5 * syspars.a;
        for i in range(500):
            delVplus, _, _  = self.gradV( rOld + epsilonr, 0.0, 0.0)
            delVminus, _, _ = self.gradV( rOld - epsilonr, 0.0, 0.0)
            slope = (delVplus - delVminus) / (2.0 * epsilonr)
            deltar = 0.5 * (delVplus + delVminus) / slope
            r = rOld - deltar
            if (r <= 0.0) or (r >= syspars.a):
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
        if( Vtarget > syspars.VL1 ): 
            print("Attempted to find a point outside the Roche lobe in findR).")

        epsilonr = 1.0e-7 * syspars.a
        r = 0.8 * syspars.rL1
        for i in range(100):
            if i ==100:
                print("Too many iterations in findR.")
            Vplus  = self.V( r + epsilonr, theta, phi)
            Vminus = self.V( r - epsilonr, theta, phi)
            Vzero  = self.V( r,           theta, phi)
            slope = (Vplus - Vminus) / (2.0 * epsilonr)
            if slope == 0.0:
                slope = 100.0 * Vzero / syspars.a
            deltar = (Vtarget - Vzero) / slope
            if( math.abs(deltar) > (0.05 * r) ): 
                deltar = (0.05 * r) * (deltar / math.abs(deltar))
            r = r + deltar 
            if( r >= syspars.rL1 ):
                r = syspars.rL1
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
        Vtarget = syspars.VL1

        if( z >= syspars.rL1 ):
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

        epsilony = 1.0e-7 * syspars.a
        y = 0.5 * syspars.rL1
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
            slope = 100.0 * Vzero / syspars.a
        deltay = (Vtarget - Vzero) / slope
        if( math.abs(deltay) > (0.05 * y) ): 
            deltay = (0.05 * y) * ( deltay / math.abs( deltay ) )
        y = y + deltay 
        if( y > (1.01 * Grid.ymax) ):
            y = 1.01 * Grid.ymax
        if( y < 0.0 ):
	       y = 0.0
        topy = y 
    
        return topy

