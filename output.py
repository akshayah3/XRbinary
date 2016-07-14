# -*- coding: utf-8 -*-
"""
All the output functions are in this file.
"""

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
    out = open(filename.lightcurves, "w")
    if( is_open(filename.lightcurves) and open_status(filename.lightcurves) == "w"):
        Quit("Cannot open lightcurve output file .LC")

    if( orbit.normalize == "OFF"):
        out.write(" NORMALIZE=        %s\n\n" % (orbit.normalize))
    else:
        if( orbit.normalize == "MAXVALUE"):
            out.write(" NORMALIZE=   %s  %12.4e\n\n" % (orbit.normalize, orbit.normvalue))
        if( orbit.normalize == "FITDATA"):
            if( orbit.normfilter == "SQUARE"): 
                lambda1 = 1.0e8 * orbit.normMinlambda
                lambda2 = 1.0e8 * orbit.normMaxlambda
                out.write(" NORMALIZE=   %s  %s  %7.1f  %7.1f\n" % orbit.normalize, orbit.normfilter, lambda1, lambda2)
            else:
	          out.write(" NORMALIZE=  %s  %s\n" %
		                   orbit.normalize, orbit.normfilter)
            out.write(" chi^2 = %12.5e\n\n" % (data.chisquare))

    outputline = "    Orbital "
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            dummy =  "      %s      " % orbit.filter[band]
            outputline += dummy
        else:
            dummy = "                %s              " % orbit.filter[band]
            outputline += dummy
    outputline += "\n"
    out.write("%s"% (outputline))

    outputline = "     Phase "
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            lambda1 = 1.0e8 * orbit.minlambda[band]
            lambda2 = 1.0e8 * orbit.maxlambda[band]
            dummy = "      (%4.0f A - %4.0f A) " % (lambda1, lambda2)
            outputline += dummy
        else:
            outputline += "                           "
    outputline += "\n"
    out.write("%s" %outputline)

    for i in range(0, orbit.maxpindex):
        outputline = "    %7.4f " % (LCphase[i])
        for band in range(1, orbit.nbands):
            dummy =  " %11.4e  %11.4e" %( LCflux[band][i], NormLC[band][i])
            outputline += dummy
        outputline += "\n"
        out.write("%s"% (outputline))
        out.close()

    return    

def WriteSysPars():
    """
    Write various system parameters to a file named syspars.out 
    """
    out = open(filename.syspars, "w")
    if( out == NULL):
        Quit("Cannot open .SysPars output file.")

    out.write( "STAR1=            %s\n"% control.star1)
    out.write( "STAR2=            %s\n"% control.star2)
    out.write( "DISK=             %s\n"% control.disk)
    out.write( "DISKRIM=          %s\n"% control.diskrim)
    out.write("DISKTORUS=        %s\n"% control.disktorus)
    out.write("INNERDISK=        %s\n"% control.innerdisk)
    out.write( "DISKSPOTS=        %s\n"% control.diskspots);
    out.write( "ADC=              %s\n"% control.adc)
    out.write("THIRDLIGHT=       %s\n"% control.thirdlight)
    out.write("IRRADIATION=      %s\n"% control.irradiation)

    out.write("\n")
    out.write("orbit.phasemin      =  %7.4f\n"%  orbit.phasemin)
    out.write("orbit.phasemax      =  %7.4f\n"%  orbit.phasemax)
    out.write("orbit.deltaphase    =  %7.4f\n"%  orbit.deltaphase)
    out.write("orbit.maxpindex     =   %3ld\n"%  orbit.maxpindex)
    out.write("orbit.phaseoffset   =  %7.4f\n"%  orbit.phaseoffset)
    out.write("  Bandpass    Type     Min Wavelength  Max Wavelength\n")
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            lambda1 = 1.0e8 * orbit.minlambda[band]
            lambda2 = 1.0e8 * orbit.maxlambda[band]
            out.write("     %2ld      %s        %5.0f           %5.0f\n"%
                   band, orbit.filter[band], lambda1, lambda2)
        else:
	      out.write( "     %2ld         %s\n"% (band, orbit.filter[band]))
    if( orbit.normalize == "OFF"):
        out.write( "NORMALIZE=        %s\n"% orbit.normalize)
    else:
        if( orbit.normalize == "MAXVALUE"):
            out.write("NORMALIZE=   %s  %12.4e\n"% (orbit.normalize,
	                                         orbit.normvalue))
        if( orbit.normalize == "FITDATA"):
            if( orbit.normfilter == "SQUARE"):
                lambda1 = 1.0e8 * orbit.normMinlambda
                lambda2 = 1.0e8 * orbit.normMaxlambda
                out.write("NORMALIZE=   %s  %s  %7.1f  %7.1f\n"% (orbit.normalize, orbit.normfilter, lambda1, lambda2))
            else:
	          out.write( "NORMALIZE=  %s  %s\n"%
		                   (orbit.normalize, orbit.normfilter))

    out.write( "\n")
    out.write("syspars.p           =  %15.8e\n"% syspars.p)
    out.write("syspars.omega       =  %15.8e\n"% syspars.omega)
    out.write("syspars.K2          =  %15.8e\n"% syspars.K2)
    out.write("syspars.q           = %8.4f\n"%  syspars.q)
    out.write("syspars.i           =  %7.4f\n"%  syspars.i)
    out.write("syspars.a           =  %15.8e\n"% syspars.a)
    out.write("syspars.zcm         =  %15.8e\n"% syspars.zcm)
    out.write("syspars.M1 (gm)     =  %10.3e\n"% syspars.M1)
    out.write("syspars.M2 (gm)     =  %10.3e\n"% syspars.M2)
    x = syspars.M1 / MSOL
    out.write("syspars.M1 (Msun)   =  %6.3f\n"%  x)
    x = syspars.M2 / MSOL
    out.write("syspars.M2 (Msun)   =  %6.3f\n"%  x)
    out.write("syspars.rL1         =  %15.8e\n"% syspars.rL1)
    out.write("syspars.VL1         =  %15.8e\n"% syspars.VL1)
    out.write("syspars.MeanLobe1Radius =  %7.4f\n"% 
                                        syspars.MeanLobe1Radius)
    out.write("syspars.MeanLobe2Radius =  %7.4f\n"% 
                                        syspars.MeanLobe2Radius)

    if( control.star1 == "ON"):
        out.write( "\n")
        out.write( "star1.L             =  %10.3e\n"% star1.L)
        out.write( "star1.T             =  %10.3e\n"% star1.T)
        out.write( "star1.sigmaT4       =  %12.5e\n"% star1.sigmaT4)
        out.write( "star1.radius        =  %12.5e\n"% star1.radius)

    if( control.star2 == "ON"):
        out.write( "\n")
        out.write("star2.targetNtiles  =   %5ld\n"%  star2.targetNtiles)
        out.write( "star2.Ntiles        =   %5ld\n"%  star2.Ntiles)
        out.write("star2.frontradius   =  %12.5e\n"% star2.frontradius)
        out.write( "star2.poleradius    =  %12.5e\n"% star2.poleradius)
        out.write( "star2.sideradius    =  %12.5e\n"% star2.sideradius)
        out.write( "star2.backradius    =  %12.5e\n"% star2.backradius)
        out.write( "star2.volume        =  %12.5e\n"% star2.volume)
        out.write( "star2.meanr         =  %12.5e\n"% star2.meanr)
        out.write( "star2.meang         =  %10.3e\n"% star2.meang)
        out.write( "star2.logg          =   %6.3f\n"% star2.logg)
        out.write( "star2.meanT         =  %5.0f\n"%  star2.meanT)
        out.write( "star2.beta          =   %5.3f\n"% star2.beta)
        out.write( "star2.albedo        =   %5.2f\n"% star2.albedo)
        out.write( "star2.L             =  %12.5e\n"% star2.L)

    if( control.star2spots == "ON"):
       out.write( "\n")
       out.write( "star2spot.nspots = %2ld\n"% star2spot.nspots)
       out.write( "   i   theta[i]  phi[i]  radius[i] T-Ratio[i]\n")
       for i in range(1, star2spot.nspots):
           x = star2spot.theta[i]  * (360.0 / 2*math.pi)
           y = star2spot.phi[i]    * (360.0 / 2*math.pi)
           z = star2spot.radius[i] * (360.0 / 2*math.pi)
           out.write( "  %2ld   %7.2f  %7.2f  %7.2f     %5.3f\n"% 
                     i, x, y, z, star2spot.SpotToverStarT[i] )

    if( control.disk == "ON"):
        out.write( "\n")
        out.write( "disk.targetNtiles   =   %5ld\n"%  disk.targetNtiles)
        out.write( "disk.Ntiles         =   %5ld\n"%  disk.Ntiles)
        out.write( "disk.e              =   %5.3f\n"% disk.e)
        x = disk.zetazero * (360.0/2*math.pi)
        out.write( "disk.zetazero       =   %5.1f\n"% x)
        out.write( "disk.albedo         =   %5.2f\n"% disk.albedo)
        out.write( "disk.L              =  %12.5e\n"% disk.L)
        out.write( "disk.TopTmax        =  %10.3e\n"% disk.TopTmax)
        out.write( "disk.TopTmin        =  %10.3e\n"% disk.TopTmin)

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
        out.write( "diskedge.T          =   %9.3e\n"% diskedge.T)
        out.write( "diskedge.Tspot      =   %9.3e\n"% diskedge.Tspot)
        x = diskedge.ZetaMid * (360.0/2*math.pi)
        out.write( "diskedge.ZetaMid    =   %5.1f\n"%  x)
        x = diskedge.ZetaWidth * (360.0/2*math.pi)
        out.write( "diskedge.ZetaWidth  =   %5.1f\n"%  x)

    if( control.innerdisk == "ON"):
        out.write( "\n")
        out.write( "innerdisk.T         =  %9.2e\n"% innerdisk.T)
        out.write( "innerdisk.L         =  %9.2e\n"% innerdisk.L)
        out.write( "innerdisk.sigmaT4   =  %9.2e\n"% innerdisk.sigmaT4)
        out.write( "innerdisk.radius    =  %9.2e\n"% innerdisk.radius)

    if( control.diskrim == "ON"):
        out.write( "\n")
        out.write( "diskrim.awidth      =  %11.4e\n"% diskrim.awidth)
        out.write( "diskrim.type        =   %s\n"%     diskrim.type)
        out.write( "diskrim.Hmax        =  %11.4e\n"% diskrim.Hmax)
        out.write( "diskrim.Hmin        =  %11.4e\n"% diskrim.Hmin)
        x = diskrim.ZetaHmax * (360.0/2*math.pi)
        out.write( "diskrim.ZetaHmax    =   %5.1f\n"%  x)
        out.write( "diskrim.Tmax        =   %9.3e\n"% diskrim.Tmax)
        out.write( "diskrim.Tmin        =   %9.3e\n"% diskrim.Tmin)
        x = diskrim.ZetaTmax * (360.0/2*math.pi)
        out.write( "diskrim.ZetaTmax    =   %5.1f\n"%  x)
        if( diskrim.type == "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, diskrim.points):
                x = diskrim.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%i, x, diskrim.PointH[i], diskrim.PointT[i])

    if( control.disktorus == "ON"):
        out.write( "\n")
        out.write( "disktorus.azero     =  %11.4e\n"% disktorus.azero)
        out.write( "disktorus.awidth    =  %11.4e\n"% disktorus.awidth)
        out.write( "disktorus.type      =   %s\n"%     diskrim.type)
        out.write( "disktorus.Hmax      =  %11.4e\n"% disktorus.Hmax)
        out.write( "disktorus.Hmin      =  %11.4e\n"% disktorus.Hmin)
        x = disktorus.ZetaHmax * (360.0/2*math.pi)
        out.write( "disktorus.ZetaHmax  =   %5.1f\n"%  x)
        out.write( "disktorus.Tmax      =   %9.3e\n"% disktorus.Tmax)
        out.write( "disktorus.Tmin      =   %9.3e\n"% disktorus.Tmin)
        x = disktorus.ZetaTmax * (360.0/2*math.pi)
        out.write( "disktorus.ZetaTmax  =   %5.1f\n"%  x)
        if( disktorus.type ==  "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, disktorus.points):
                x = disktorus.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%i, x, disktorus.PointH[i], disktorus.PointT[i])

    if( control.diskspots == "ON"):
        out.write( "\n")
        out.write("diskspot.npoints = %2ld\n"% diskspot.nspots)
        out.write( "   i  ZetaMin[i]  ZetaMax[i]   aMin[i]     aMax[i]  spotToverT[i] \n")
        for i in range(1, diskspot.nspots):
            x = diskspot.zetamin[i] * (360.0 / 2*math.pi)
            y = diskspot.zetamax[i] * (360.0 / 2*math.pi)
            out.write( "  %2ld    %5.1f       %5.1f    %10.4e  %10.4e    %5.2f\n"%i, x, y, diskspot.amin[i], diskspot.amax[i], diskspot.spotToverT[i])

    if( control.adc == "ON"):
        out.write( "\n")
        out.write( "adc.L               =  %10.3e\n"% adc.L)
        out.write( "adc.height          =  %11.4e\n"% adc.height)

    if ( control.thirdlight == "ON"):
        out.write( "\n")
        out.write( "thirdlight.orbphase =  %6.3f\n"%  thirdlight.orbphase)
        out.write( "Third Light Fractions:\n")
        out.write( "  Bandpass  Type   Min Wavelength  Max Wavelength  Fraction\n")
        for band in range(1, thirdlight.nbands):
            if( thirdlight.filter[band] == "SQUARE"):
                lambda1 = 1.0e8 * thirdlight.minlambda[band]
                lambda2 = 1.0e8 * thirdlight.maxlambda[band]
                out.write( "     %2ld    %s       %5.0f           %5.0f        %5.3f\n"%
                      band, thirdlight.filter[band], 
                      lambda1, lambda2, thirdlight.fraction[band])
            else:
	          out.write( "     %2ld       %s                                      %5.3f\n"% 
                      band, thirdlight.filter[band],
                      thirdlight.fraction[band])

    if( data.nbands > 0):
        out.write( "\n")
        out.write( "Observed Light Curve Data Files:\n")
        out.write( "                      Minimum     Maximum     Data\n")
        out.write( "  Bandpass   Type    Wavelength  Wavelength  Points    Filename\n")
        for band in range(1, data.nbands):
            if data.filter[band] == "SQUARE":
                lambda1 = 1.0e8 * data.minlambda[band]
                lambda2 = 1.0e8 * data.maxlambda[band]
                out.write( "     %2ld     %s     %5.0f       %5.0f       %3ld   %s\n"%
                      band, data.filter[band], lambda1, lambda2,
                      data.npoints[band], data.filename[band])
            else:
	          out.write( "     %2ld        %s                               %3ld   %s\n"% 
                         band, data.filter[band], 
                         data.npoints[band], data.filename[band])

    out.close()

    return