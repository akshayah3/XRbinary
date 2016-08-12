# -*- coding: utf-8 -*-
"""

"""
SIGMA = 5.6704e-5
import math
import sys
from .diskflux import maindisk
from .star1 import Star1
from .star2 import Star2
from .parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, XYGrid, dataparams, ReadInput

def main():
    if( sys.argv != 2 ):
        sys.exit("Wrong number of command line parameters.")

    Initialize( sys.argv[0] )

    ReadInput()

    CalcSysPars()

    if( flowcontrol.star2 == "ON" ):
        MakeStar2Tiles()

    if( flowcontrol.disk == "ON"):
        MakeDiskTiles()

    MakeYlimits()
    if( flowcontrol.irradiation == "ON"):
        Irradiate()
    star2.L = Star2L()
    disk.L  = DiskL()

    MakeLightCurves()

    WriteResults()

    return

def Initialize(parfilename ):
    """
    /*****************************************************

    Initialize the input parameters to impossible or
    default values so that nothing gets set by accident.

     ***************************************************/
    """
    filenames.parfile = parfilename
    filenames.syspars = parfilename
    filenames.syspars += ".SysPars"
    filenames.lightcurves = parfilename
    filenames.lightcurves += ".LC"

    verbose = "OFF"
    flowcontrol.diagnostics = "OFF"
    flowcontrol.diagnosephase = -1.0

    flowcontrol.star1 = "MISSING"
    flowcontrol.star2 = "MISSING"
    flowcontrol.disk = "MISSING"
    flowcontrol.diskrim = "MISSING"
    flowcontrol.disktorus = "MISSING"
    flowcontrol.innerdisk = "MISSING"
    flowcontrol.diskspots = "MISSING"
    flowcontrol.adc = "MISSING"
    flowcontrol.thirdlight = "MISSING"
    flowcontrol.irradiation = "MISSING"

    maxGDindex          = -1
    maxLDgindex         = -1
    maxLDTindex         = -1
    maxLDfilterindex    = -1
    maxIperpgindex      = -1
    maxIperpTindex      = -1
    maxIperpfilterindex = -1
    maxIBBTindex        = -1
    maxIBBfilterindex   = -1
    IBBTmin             = -1.0
    IBBTmax             = -1.0
    IBBdeltaT           = -1.0
    maxBBzetaindex      = -1

    orbitparams.phasemin      = -1.0
    orbitparams.phasemax      = -1.0
    orbitparams.deltaphase    = -1.0
    orbitparams.maxpindex     = -1
    orbitparams.phaseoffset   = -1.0
    orbitparams.nbands        =  0
    orbitparams.normalize = "MISSING"
    orbitparams.normvalue  = -1.0

    systemparams.p           = -1.0
    systemparams.omega       = -1.0
    systemparams.K2          = -1.0
    systemparams.q           = -1.0
    systemparams.i           = -1.0
    systemparams.a           = -1.0
    systemparams.zcm         = -1.0
    systemparams.M1          = -1.0
    systemparams.M2          = -1.0
    systemparams.rL1         = -1.0
    systemparams.VL1         = -1.0

    Star1.L             = -1.0
    Star1.T             = -1.0
    Star1.sigmaT4       = -1.0
    Star1.radius        = -1.0

    Star2.targetNtiles  = -1
    Star2.Ntiles        = -1
    Star2.volume        = -1.0
    Star2.meanr         = -1.0
    Star2.meang         =  0.0
    Star2.logg          = -1.0
    Star2.meanT         = -1.0
    Star2.beta          = -1.0
    Star2.albedo        = -1.0
    Star2.L             = -1.0
    Star2.frontradius   = -1.0
    Star2.poleradius    = -1.0
    Star2.sideradius    = -1.0
    Star2.backradius    = -1.0

    star2spotparams.nspots    =  0

    wholediskpars.targetNtiles   = -1
    wholediskpars.Ntiles         = -1
    wholediskpars.e              = -1.0
    wholediskpars.zetazero       = -1.0
    wholediskpars.albedo         = -1.0
    wholediskpars.L              = -1.0
    wholediskpars.TopTmax        = -1.0
    wholediskpars.TopTmin        = -1.0

    maindisk.amin       = -1.0
    maindisk.amax       = -1.0
    maindisk.Hmax       = -1.0
    maindisk.Hpow       = -1.0
    maindisk.Ttype = "MISSING"
    maindisk.maindiskL  = -1.0
    maindisk.Tpow       = -10.0
    maindisk.TConstant  = -1.0
    maindisk.Tamax      = -1.0
    maindisk.Tamin      = -1.0

    diskedgepars.T          = -1.0
    diskedgepars.Tspot      = -1.0
    diskedgepars.ZetaMid    = -1.0
    diskedgepars.ZetaWidth  = -1.0

    diskrimpars.type = "MISSING"
    diskrimpars.awidth      = -1.0
    diskrimpars.Hmax        = -1.0
    diskrimpars.Hmin        = -1.0
    diskrimpars.ZetaHmax    = -1.0
    diskrimpars.Tmax        = -1.0
    diskrimpars.Tmin        = -1.0
    diskrimpars.ZetaTmax    = -1.0
    diskrimpars.points      =  0

    disktorusparams.type = "MISSING"
    disktorusparams.azero     = -1.0
    disktorusparams.awidth    = -1.0
    disktorusparams.Hmax      = -1.0
    disktorusparams.Hmin      = -1.0
    disktorusparams.ZetaHmax  = -1.0
    disktorusparams.Tmax      = -1.0
    disktorusparams.Tmin      = -1.0
    disktorusparams.ZetaTmax  = -1.0
    disktorusparams.points    =  0

    innerdiskpars.T         = -1.0
    innerdiskpars.L         = -1.0
    
    diskspotpars.nspots     = 0

    adcpars.L               = -1.0
    adcpars.height          = -1.0

    thirdlightparams.nbands   = 0
    thirdlightparams.orbphase = -100.0

    XYGrid.Nxtiles        = -1
    XYGrid.Nztiles        = -1

    dataparams.nbands         = 0

    return

def CalcSysPars():
    """    
    /****************************************************

    This function uses the various input parameters to calculate
    other useful system parameters.  It also converts some
    other input data to cgs units.

    *********************************************************/
    """
    if flowcontrol.diagnostics != "OFF":
        if( (flowcontrol.diagnosephase < orbitparams.phasemin)
	         or (flowcontrol.diagnosephase > orbitparams.phasemax) ):
            sys.exit("diagnosephase must be ge phasemin and le phasemax.")
        flowcontrol.diagnoseindex = 0.5 + ( flowcontrol.diagnosephase - orbitparams.phasemin ) / orbitparams.deltaphase

    orbitparams.maxpindex    = 1.0e-7 + (orbitparams.phasemax - orbitparams.phasemin) / orbitparams.deltaphase
    for band in range(1, orbitparams.nbands):
        orbitparams.minlambda[band] *= 1.0e-8
        orbitparams.maxlambda[band] *= 1.0e-8
    if( orbitparams.normalize == "FITDATA"):
        if( orbitparams.normfilter == "SQUARE"):
            orbitparams.normMinlambda *= 1.0e-8
            orbitparams.normMaxlambda *= 1.0e-8
    systemparams.p          *= 86400.0;
    systemparams.omega      = math.pi*2 / systemparams.p
    systemparams.i          *= math.pi*2 / 360.0
    if( systemparams.K2 > 0.0 ):
        systemparams.K2         *= 1.0e5;
        V2orb              = systemparams.K2 / math.sin( systemparams.i )
        systemparams.a          = V2orb * (1.0 + systemparams.q) / systemparams.omega
        M1plusM2           = pow( systemparams.omega, 2) * pow( systemparams.a, 3) / G
        systemparams.M1         = M1plusM2 / (1.0 + systemparams.q)
        systemparams.M2         = systemparams.M1 * systemparams.q
    else:
        systemparams.M1      = systemparams.M1 * MSOL
        systemparams.M2      = systemparams.M1 * systemparams.q
        M1plusM2        = systemparams.M1 + systemparams.M2
        a3              = (G * M1plusM2 * systemparams.p * systemparams.p) / (4.0 * math.pi * math.pi )
        systemparams.a       = pow( a3, 0.3333333333333 )
        a_2             = systemparams.a / (1.0 + systemparams.q)
        V2orb           = systemparams.omega * a_2 / 1.0e5
        systemparams.K2      = V2orb * math.sin( systemparams.i )
    systemparams.zcm        = systemparams.a / (1.0 + systemparams.q)
    systemparams.rL1        = FindL1()  
    systemparams.VL1        = V( systemparams.rL1, 0.0, 0.0 )
    systemparams.MeanLobe2Radius = MeanRocheRadius( systemparams.q )
    systemparams.MeanLobe1Radius = MeanRocheRadius( 1.0 / systemparams.q )

    Star1.sigmaT4      = SIGMA * pow( Star1.T, 4.0 )
    r2                 = Star1.L / ( 4.0 * math.pi * Star1.sigmaT4 )
    Star1.radius       = math.sqrt( r2 )

    if( star2spotparams.nspots > 0 ):
        for i in range (1, star2spotparams.nspots):
            star2spotparams.theta[i]  *= 2*math.pi / 360.0
            star2spotparams.phi[i]    *= 2*math.pi / 360.0
            star2spotparams.radius[i] *= 2*math.pi / 360.0
    wholediskpars.zetazero      *= math.pi*2 / 360.0
    maindisk.amin      *= systemparams.a
    maindisk.amax      *= systemparams.a
    maindisk.Hmax      *= systemparams.a

    diskedgepars.ZetaMid   *= 2*math.pi / 360.0
    diskedgepars.ZetaWidth *= 2*math.pi / 360.0

    if( flowcontrol.diskrim == "ON"):
        if( diskrimpars.type == "POINT"):
            MakeDiskRim()
            for i in range(1, diskrimpars.points):
                diskrimpars.PointZeta[i] *= 2*math.pi / 360.0
                diskrimpars.PointH[i]    *= systemparams.a
        diskrimpars.awidth     *= systemparams.a
        diskrimpars.Hmax       *= systemparams.a
        diskrimpars.Hmin       *= systemparams.a
        diskrimpars.ZetaHmax   *= math.pi*2 / 360.0
        diskrimpars.ZetaTmax   *= math.pi*2 / 360.0

    if( flowcontrol.disktorus == "ON"):
        if( disktorusparams.type == "POINT"):
            MakeDiskTorus()
            for i in range(1, disktorusparams.points):
                disktorusparams.PointZeta[i] *= 2*math.pi / 360.0
                disktorusparams.PointH[i]    *= systemparams.a
        disktorusparams.azero    *= systemparams.a
        disktorusparams.awidth   *= systemparams.a
        disktorusparams.Hmax     *= systemparams.a
        disktorusparams.Hmin     *= systemparams.a
        disktorusparams.ZetaHmax *= 2*math.pi / 360.0
        disktorusparams.ZetaTmax *= 2*math.pi / 360.0

    if( flowcontrol.innerdisk == "ON"):
        innerdiskpars.sigmaT4      = SIGMA * pow( innerdiskpars.T, 4.0 )
        r2                     = innerdiskpars.L / ( 2.0 * math.pi * innerdiskpars.sigmaT4 )
        innerdiskpars.radius       = math.sqrt( r2 )

    if( diskspotpars.nspots > 0 ):
      for i in range(1, diskspotpars.nspots):
          if( diskspotpars.zetamin[i] > diskspotpars.zetamax[i] ):
              diskspotpars.nspots += 1
              if( diskspotpars.nspots >= 20 ):
	             sys.exit("diskspot.nspots too large in function CalcSysPars.")
              for j in range(diskspotpars.nspots, i, -1):
                  diskspotpars.zetamin[j] = diskspotpars.zetamin[j-1]
                  diskspotpars.zetamax[j] = diskspotpars.zetamax[j-1]
                  diskspotpars.amin[j] = diskspotpars.amin[j-1]
                  diskspotpars.amax[j] = diskspotpars.amax[j-1]
                  diskspotpars.spotToverT[j] = diskspotpars.spotToverT[j-1]
               
              diskspotpars.zetamax[i+1] = 360.0
              diskspotpars.zetamin[i+1] = diskspotpars.zetamin[i]
              diskspotpars.zetamax[i] = diskspotpars.zetamax[i]
              diskspotpars.zetamin[i] = 0.0
          
          diskspotpars.zetamin[i] *= math.pi*2 / 360.0
          diskspotpars.zetamax[i] *= math.pi*2 / 360.0
          diskspotpars.amin[i]    *= syspars.a
          if( diskspotpars.amin[i] < maindisk.amin ):
              diskspotpars.amin[i] = maindisk.amin
          diskspotpars.amax[i]    *= systemparams.a
          if( diskspotpars.amax[i] > maindisk.amax ):
              diskspotpars.amax[i] = maindisk.amax

    if( flowcontrol.adc == "ON"):
        adcpars.height *= systemparams.a

    if( flowcontrol.thirdlight == "ON"):
        for band in range(1, thirdlightparams.nbands):
            thirdlightparams.minlambda[band] *= 1.0e-8
            thirdlightparams.maxlambda[band] *= 1.0e-8

    for band in range(1, dataparams.nbands):
       dataparams.minlambda[band] *= 1.0e-8
       dataparams.maxlambda[band] *= 1.0e-8

    if( flowcontrol.diagnostics == "INSPECTSYSPARS"):
       WriteSysPars()
       sys.exit("Quit in main.c after INSPECTSYSPARS.")
    return
