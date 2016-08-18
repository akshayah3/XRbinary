# -*- coding: utf-8 -*-
"""
All the output functions are in this file.
"""
MSOL  =   1.9891e33
import math
from diskflux import maindisk
from star1 import Star1
from star2 import Star2
from parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, dataparams, globalvar

def Quit(outputline):
    print("\n   {}\n", outputline)
    print("   Program terminated.\n")
    exit( 0 )
    return

def WriteResults():
    """
    This function writes out the results of the calculations.
    """
    WriteLightCurves()
    WriteSysPars()

    return

def WriteLightCurves():
    """
    This function writes out the orbital light curve into the 
    file litecurves.out
    """
    out = open(filenames.lightcurves, "w")
    if( out == None):
        Quit("Cannot open lightcurve output file .LC")

    if( orbitparams.normalize == "OFF"):
        out.write(" NORMALIZE=        {}\n\n".format( (orbitparams.normalize)))
    else:
        if( orbitparams.normalize == "MAXVALUE"):
            out.write(" NORMALIZE=   {}  {}\n\n".format( (orbitparams.normalize, orbitparams.normvalue)))
        if( orbitparams.normalize == "FITDATA"):
            if( orbitparams.normfilter == "SQUARE"): 
                lambda1 = 1.0e8 * orbitparams.normMinlambda
                lambda2 = 1.0e8 * orbitparams.normMaxlambda
                out.write(" NORMALIZE=   {}  {}  {}  {}\n".format( orbitparams.normalize, orbitparams.normfilter, lambda1, lambda2))
            else:
	          out.write(" NORMALIZE=  {}  {}\n".format(
		                   orbitparams.normalize, orbitparams.normfilter))
            out.write(" chi^2 = {}\n\n".format( (dataparams.chisquare)))

    outputline = "    Orbital "
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            dummy =  "      %s      " % orbitparams.filter[band]
            outputline += dummy
        else:
            dummy = "                %s              " % orbitparams.filter[band]
            outputline += dummy
    outputline += "\n"
    out.write("%s"% (outputline))

    outputline = "     Phase "
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            lambda1 = 1.0e8 * orbitparams.minlambda[band]
            lambda2 = 1.0e8 * orbitparams.maxlambda[band]
            dummy = "      ({} A - {} A) ".format(lambda1, lambda2)
            outputline += dummy
        else:
            outputline += "                           "
    outputline += "\n"
    out.write("%s" %outputline)

    for i in range(0, orbitparams.maxpindex):
        outputline = "    {} ".format(globalvar.LCphase[i])
        for band in range(1, orbitparams.nbands):
            dummy =  " {}  {}".format( globalvar.LCflux[band][i], globalvar.NormLC[band][i])
            outputline += dummy
        outputline += "\n"
        out.write("%s"% (outputline))
        out.close()

    return    

def WriteSysPars():
    """
    Write various system parameters to a file named syspars.out 
    """
    out = open(filenames.syspars, "w")
    if( out == None):
        Quit("Cannot open .SysPars output file.")

    out.write( "STAR1=            %s\n"% flowcontrol.star1)
    out.write( "STAR2=            %s\n"% flowcontrol.star2)
    out.write( "DISK=             %s\n"% flowcontrol.disk)
    out.write( "DISKRIM=          %s\n"% flowcontrol.diskrim)
    out.write("DISKTORUS=        %s\n"% flowcontrol.disktorus)
    out.write("INNERDISK=        %s\n"% flowcontrol.innerdisk)
    out.write( "DISKSPOTS=        %s\n"% flowcontrol.diskspots);
    out.write( "ADC=              %s\n"% flowcontrol.adc)
    out.write("THIRDLIGHT=       %s\n"% flowcontrol.thirdlight)
    out.write("IRRADIATION=      %s\n"% flowcontrol.irradiation)

    out.write("\n")
    out.write("orbit.phasemin      =  %7.4f\n"%  orbitparams.phasemin)
    out.write("orbit.phasemax      =  %7.4f\n"%  orbitparams.phasemax)
    out.write("orbit.deltaphase    =  %7.4f\n"%  orbitparams.deltaphase)
    out.write("orbit.maxpindex     =   %3ld\n"%  orbitparams.maxpindex)
    out.write("orbit.phaseoffset   =  %7.4f\n"%  orbitparams.phaseoffset)
    out.write("  Bandpass    Type     Min Wavelength  Max Wavelength\n")
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            lambda1 = 1.0e8 * orbitparams.minlambda[band]
            lambda2 = 1.0e8 * orbitparams.maxlambda[band]
            out.write("     %2ld      %s        %5.0f           %5.0f\n"%
                   band, orbitparams.filter[band], lambda1, lambda2)
        else:
	      out.write( "     %2ld         %s\n"% (band, orbitparams.filter[band]))
    if( orbitparams.normalize == "OFF"):
        out.write( "NORMALIZE=        %s\n"% orbitparams.normalize)
    else:
        if( orbitparams.normalize == "MAXVALUE"):
            out.write("NORMALIZE=   %s  %12.4e\n"% (orbitparams.normalize,
	                                         orbitparams.normvalue))
        if( orbitparams.normalize == "FITDATA"):
            if( orbitparams.normfilter == "SQUARE"):
                lambda1 = 1.0e8 * orbitparams.normMinlambda
                lambda2 = 1.0e8 * orbitparams.normMaxlambda
                out.write("NORMALIZE=   %s  %s  %7.1f  %7.1f\n"% (orbitparams.normalize, orbitparams.normfilter, lambda1, lambda2))
            else:
	          out.write( "NORMALIZE=  %s  %s\n"%
		                   (orbitparams.normalize, orbitparams.normfilter))

    out.write( "\n")
    out.write("syspars.p           =  %15.8e\n"% systemparams.p)
    out.write("syspars.omega       =  %15.8e\n"% systemparams.omega)
    out.write("syspars.K2          =  %15.8e\n"% systemparams.K2)
    out.write("syspars.q           = %8.4f\n"%  systemparams.q)
    out.write("syspars.i           =  %7.4f\n"%  systemparams.i)
    out.write("syspars.a           =  %15.8e\n"% systemparams.a)
    out.write("syspars.zcm         =  %15.8e\n"% systemparams.zcm)
    out.write("syspars.M1 (gm)     =  %10.3e\n"% systemparams.M1)
    out.write("syspars.M2 (gm)     =  %10.3e\n"% systemparams.M2)
    x = systemparams.M1 / MSOL
    out.write("syspars.M1 (Msun)   =  %6.3f\n"%  x)
    x = systemparams.M2 / MSOL
    out.write("syspars.M2 (Msun)   =  %6.3f\n"%  x)
    out.write("syspars.rL1         =  %15.8e\n"% systemparams.rL1)
    out.write("syspars.VL1         =  %15.8e\n"% systemparams.VL1)
    out.write("syspars.MeanLobe1Radius =  %7.4f\n"% 
                                        systemparams.MeanLobe1Radius)
    out.write("syspars.MeanLobe2Radius =  %7.4f\n"% 
                                        systemparams.MeanLobe2Radius)

    if( flowcontrol.star1 == "ON"):
        out.write( "\n")
        out.write( "star1.L             =  %10.3e\n"% Star1.L)
        out.write( "star1.T             =  %10.3e\n"% Star1.T)
        out.write( "star1.sigmaT4       =  %12.5e\n"% Star1.sigmaT4)
        out.write( "star1.radius        =  %12.5e\n"% Star1.radius)

    if( flowcontrol.star2 == "ON"):
        out.write( "\n")
        out.write("star2.targetNtiles  =   %5ld\n"%  Star2.targetNtiles)
        out.write( "star2.Ntiles        =   %5ld\n"%  Star2.Ntiles)
        out.write("star2.frontradius   =  %12.5e\n"% Star2.frontradius)
        out.write( "star2.poleradius    =  %12.5e\n"% Star2.poleradius)
        out.write( "star2.sideradius    =  %12.5e\n"% Star2.sideradius)
        out.write( "star2.backradius    =  %12.5e\n"% Star2.backradius)
        out.write( "star2.volume        =  %12.5e\n"% Star2.volume)
        out.write( "star2.meanr         =  %12.5e\n"% Star2.meanr)
        out.write( "star2.meang         =  %10.3e\n"% Star2.meang)
        out.write( "star2.logg          =   %6.3f\n"% Star2.logg)
        out.write( "star2.meanT         =  %5.0f\n"%  Star2.meanT)
        out.write( "star2.beta          =   %5.3f\n"% Star2.beta)
        out.write( "star2.albedo        =   %5.2f\n"% Star2.albedo)
        out.write( "star2.L             =  %12.5e\n"% Star2.L)

    if( flowcontrol.star2spots == "ON"):
       out.write( "\n")
       out.write( "star2spot.nspots = %2ld\n"% star2spotparams.nspots)
       out.write( "   i   theta[i]  phi[i]  radius[i] T-Ratio[i]\n")
       for i in range(1, star2spotparams.nspots):
           x = star2spotparams.theta[i]  * (360.0 / 2*math.pi)
           y = star2spotparams.phi[i]    * (360.0 / 2*math.pi)
           z = star2spotparams.radius[i] * (360.0 / 2*math.pi)
           out.write( "  %2ld   %7.2f  %7.2f  %7.2f     %5.3f\n"% 
                     i, x, y, z, star2spotparams.SpotToverStarT[i] )

    if( flowcontrol.disk == "ON"):
        out.write( "\n")
        out.write( "disk.targetNtiles   =   %5ld\n"%  wholediskpars.targetNtiles)
        out.write( "disk.Ntiles         =   %5ld\n"%  wholediskpars.Ntiles)
        out.write( "disk.e              =   %5.3f\n"% wholediskpars.e)
        x = wholediskpars.zetazero * (360.0/2*math.pi)
        out.write( "disk.zetazero       =   %5.1f\n"% x)
        out.write( "disk.albedo         =   %5.2f\n"% wholediskpars.albedo)
        out.write( "disk.L              =  %12.5e\n"% wholediskpars.L)
        out.write( "disk.TopTmax        =  %10.3e\n"% wholediskpars.TopTmax)
        out.write( "disk.TopTmin        =  %10.3e\n"% wholediskpars.TopTmin)

        out.write( "\n")
        out.write( "maindisk.amin       =  %11.4e\n"% maindisk.amin)
        out.write( "maindisk.amax       =  %11.4e\n"% maindisk.amax)
        out.write( "maindisk.Hmax       =  %11.4e\n"% maindisk.Hmax)
        out.write( "maindisk.Hpow       =  %5.2f\n"%  maindisk.Hpow)
        out.write( "maindisk.Ttype      =  %s\n"%     maindisk.Ttype)
        if( maindisk.Ttype =="POWER"):
            out.write( "maindisk.Tpow       =  %5.2f\n"% maindisk.Tpow)
        out.write( "maindisk.maindiskL  =  %10.3e\n"% maindisk.maindiskL)
        out.write( "maindisk.Tamax      =  %10.3e\n"% maindisk.Tamax)
        out.write( "maindisk.Tamin      =  %10.3e\n"% maindisk.Tamin)

        out.write( "\n")
        out.write( "diskedge.T          =   %9.3e\n"% diskedgepars.T)
        out.write( "diskedge.Tspot      =   %9.3e\n"% diskedgepars.Tspot)
        x = diskedgepars.ZetaMid * (360.0/2*math.pi)
        out.write( "diskedge.ZetaMid    =   %5.1f\n"%  x)
        x = diskedgepars.ZetaWidth * (360.0/2*math.pi)
        out.write( "diskedge.ZetaWidth  =   %5.1f\n"%  x)

    if( flowcontrol.innerdisk == "ON"):
        out.write( "\n")
        out.write( "innerdisk.T         =  %9.2e\n"% innerdiskpars.T)
        out.write( "innerdisk.L         =  %9.2e\n"% innerdiskpars.L)
        out.write( "innerdisk.sigmaT4   =  %9.2e\n"% innerdiskpars.sigmaT4)
        out.write( "innerdisk.radius    =  %9.2e\n"% innerdiskpars.radius)

    if( flowcontrol.diskrim == "ON"):
        out.write( "\n")
        out.write( "diskrim.awidth      =  %11.4e\n"% diskrimpars.awidth)
        out.write( "diskrim.type        =   %s\n"%     diskrimpars.type)
        out.write( "diskrim.Hmax        =  %11.4e\n"% diskrimpars.Hmax)
        out.write( "diskrim.Hmin        =  %11.4e\n"% diskrimpars.Hmin)
        x = diskrimpars.ZetaHmax * (360.0/2*math.pi)
        out.write( "diskrim.ZetaHmax    =   %5.1f\n"%  x)
        out.write( "diskrim.Tmax        =   %9.3e\n"% diskrimpars.Tmax)
        out.write( "diskrim.Tmin        =   %9.3e\n"% diskrimpars.Tmin)
        x = diskrimpars.ZetaTmax * (360.0/2*math.pi)
        out.write( "diskrim.ZetaTmax    =   %5.1f\n"%  x)
        if( diskrimpars.type == "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, diskrimpars.points):
                x = diskrimpars.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%i, x, diskrimpars.PointH[i], diskrimpars.PointT[i])

    if( flowcontrol.disktorus == "ON"):
        out.write( "\n")
        out.write( "disktorus.azero     =  %11.4e\n"% disktorusparams.azero)
        out.write( "disktorus.awidth    =  %11.4e\n"% disktorusparams.awidth)
        out.write( "disktorus.type      =   %s\n"%     disktorusparams.type)
        out.write( "disktorus.Hmax      =  %11.4e\n"% disktorusparams.Hmax)
        out.write( "disktorus.Hmin      =  %11.4e\n"% disktorusparams.Hmin)
        x = disktorusparams.ZetaHmax * (360.0/2*math.pi)
        out.write( "disktorus.ZetaHmax  =   %5.1f\n"%  x)
        out.write( "disktorus.Tmax      =   %9.3e\n"% disktorusparams.Tmax)
        out.write( "disktorus.Tmin      =   %9.3e\n"% disktorusparams.Tmin)
        x = disktorusparams.ZetaTmax * (360.0/2*math.pi)
        out.write( "disktorus.ZetaTmax  =   %5.1f\n"%  x)
        if( disktorusparams.type ==  "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, disktorusparams.points):
                x = disktorusparams.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%i, x, disktorusparams.PointH[i], disktorusparams.PointT[i])

    if( flowcontrol.diskspots == "ON"):
        out.write( "\n")
        out.write("diskspot.npoints = %2ld\n"% diskspotpars.nspots)
        out.write( "   i  ZetaMin[i]  ZetaMax[i]   aMin[i]     aMax[i]  spotToverT[i] \n")
        for i in range(1, diskspotpars.nspots):
            x = diskspotpars.zetamin[i] * (360.0 / 2*math.pi)
            y = diskspotpars.zetamax[i] * (360.0 / 2*math.pi)
            out.write( "  %2ld    %5.1f       %5.1f    %10.4e  %10.4e    %5.2f\n"%i, x, y, diskspotpars.amin[i], diskspotpars.amax[i], diskspotpars.spotToverT[i])

    if( flowcontrol.adc == "ON"):
        out.write( "\n")
        out.write( "adc.L               =  %10.3e\n"% adcpars.L)
        out.write( "adc.height          =  %11.4e\n"% adcpars.height)

    if ( flowcontrol.thirdlight == "ON"):
        out.write( "\n")
        out.write( "thirdlight.orbphase =  %6.3f\n"%  thirdlightparams.orbphase)
        out.write( "Third Light Fractions:\n")
        out.write( "  Bandpass  Type   Min Wavelength  Max Wavelength  Fraction\n")
        for band in range(1, thirdlightparams.nbands):
            if( thirdlightparams.filter[band] == "SQUARE"):
                lambda1 = 1.0e8 * thirdlightparams.minlambda[band]
                lambda2 = 1.0e8 * thirdlightparams.maxlambda[band]
                out.write( "     %2ld    %s       %5.0f           %5.0f        %5.3f\n"%
                      band, thirdlightparams.filter[band], 
                      lambda1, lambda2, thirdlightparams.fraction[band])
            else:
	          out.write( "     %2ld       %s                                      %5.3f\n"% 
                      band, thirdlightparams.filter[band],
                      thirdlightparams.fraction[band])

    if( dataparams.nbands > 0):
        out.write( "\n")
        out.write( "Observed Light Curve Data Files:\n")
        out.write( "                      Minimum     Maximum     Data\n")
        out.write( "  Bandpass   Type    Wavelength  Wavelength  Points    Filename\n")
        for band in range(1, dataparams.nbands):
            if dataparams.filter[band] == "SQUARE":
                lambda1 = 1.0e8 * dataparams.minlambda[band]
                lambda2 = 1.0e8 * dataparams.maxlambda[band]
                out.write( "     %2ld     %s     %5.0f       %5.0f       %3ld   %s\n"%
                      band, dataparams.filter[band], lambda1, lambda2,
                      dataparams.npoints[band], dataparams.filename[band])
            else:
	          out.write( "     %2ld        %s                               %3ld   %s\n"% 
                         band, dataparams.filter[band], 
                         dataparams.npoints[band], dataparams.filename[band])

    out.close()

    return
