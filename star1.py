# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 15:47:31 2016

"""
import math
from utility import BBFilterIntensity, BBSquareIntensity

sigma = 5.6704e-5

class Star1:
    L = None
    T = None
    radius = None
    sigmaT4 = None
    """
    Star 1 is a small sphere if it represents neutron star, or an invisible point mass if it
    represents a black hole.  In either case, it is much smaller than all other dimensions in the
    system and is unresolved by any of the grids used in calculating the light curves.
    
    Parameters
    ==========
    
    L: Luminosity of the compact Star.
    
    T: Temperature of the compact Star.
    
    
    """

    def __init__(self, L, T):
        self.L = L
        self.T = T

    @property
    def luminosity(self):
        """
        Returns the luminosity of the compact Star.
        """
        return self.L

    @property
    def temperature(self):
        """
        Returns the temperature of the compact Star.
        """
        return self.T

    @property
    def radius(self):
        """
        Returns the radius of the compact star.
        """
        radius = self.L/((((self.T)**4) * sigma)* math.pi)
        return radius

    def Star1TotFlux(self, distance):
        """
        This function returns the integrated flux from star 1.
        The flux has been integrated over wavelength and over 
        the surface of the star.
        """
        totalflux = (self.L/ 4*math.pi*distance*distance)
        return totalflux
        
    def MeanRocheRadius(self, q):
        """
        This function returns the mean radius of the Roche lobe
        in units of the separation of the two stars, e.g., <R_lobe>/a.
        using Eggleton's (1983, ApJ, 268, 368) formula.
        q = (Mass of star inside the Roche lobe) / (Mass of the other star)
        """
        x = math.pow(q, 0.33)
        radius = 0.49 * x * x / ( 0.6 * x * x + math.log(1 + x) )
        return radius

    def Star1Flambda(self, minlambda, maxlambda, filter = "Square"):
        """
        This function returns the contribution of star 1 to the observed 
        spectrum at wavelength lambda.  Note that the wavelengths must be in 
        centimeters but the returned flux is in ergs/sec/cm^2/Angstrom.
        """        
        
        if filter == "Square":
            flux = ( self.L / (4.0 * (self.T)**4) ) * BBSquareIntensity( self.T, minlambda, maxlambda )
        else:
            flux = ( self.L / (4.0 * (self.T)**4) ) * BBFilterIntensity( self.T, filter )
        return flux