# -*- coding: utf-8 -*-
"""

"""

import sys

def main():
    if( sys.argv != 2 ):
        Quit("Wrong number of command line parameters.")

    Initialize( sys.argv[0] )

    ReadInput()

    CalcSysPars()

    if( control.star2 == "ON" ):
        MakeStar2Tiles()

    if( control.disk == "ON"):
        MakeDiskTiles()

    MakeYlimits()
    if( control.irradiation == "ON"):
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
    filename.parfile = parfilename
    filename.syspars = parfilename
    filename.syspars += ".SysPars"
    filename.lightcurves = parfilename
    filename.lightcurves += ".LC"

    verbose = "OFF"
    control.diagnostics = "OFF"
    control.diagnosephase = -1.0

    control.star1 = "MISSING"
    control.star2 = "MISSING"
    control.disk = "MISSING"
    control.diskrim = "MISSING"
    control.disktorus = "MISSING"
    control.innerdisk = "MISSING"
    control.diskspots = "MISSING"
    control.adc = "MISSING"
    control.thirdlight = "MISSING"
    control.irradiation = "MISSING"

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

    orbit.phasemin      = -1.0
    orbit.phasemax      = -1.0
    orbit.deltaphase    = -1.0
    orbit.maxpindex     = -1
    orbit.phaseoffset   = -1.0
    orbit.nbands        =  0
    orbit.normalize = "MISSING"
    orbit.normvalue  = -1.0

    syspars.p           = -1.0
    syspars.omega       = -1.0
    syspars.K2          = -1.0
    syspars.q           = -1.0
    syspars.i           = -1.0
    syspars.a           = -1.0
    syspars.zcm         = -1.0
    syspars.M1          = -1.0
    syspars.M2          = -1.0
    syspars.rL1         = -1.0
    syspars.VL1         = -1.0

    star1.L             = -1.0
    star1.T             = -1.0
    star1.sigmaT4       = -1.0
    star1.radius        = -1.0

    star2.targetNtiles  = -1
    star2.Ntiles        = -1
    star2.volume        = -1.0
    star2.meanr         = -1.0
    star2.meang         =  0.0
    star2.logg          = -1.0
    star2.meanT         = -1.0
    star2.beta          = -1.0
    star2.albedo        = -1.0
    star2.L             = -1.0
    star2.frontradius   = -1.0
    star2.poleradius    = -1.0
    star2.sideradius    = -1.0
    star2.backradius    = -1.0

    star2spot.nspots    =  0

    disk.targetNtiles   = -1
    disk.Ntiles         = -1
    disk.e              = -1.0
    disk.zetazero       = -1.0
    disk.albedo         = -1.0
    disk.L              = -1.0
    disk.TopTmax        = -1.0
    disk.TopTmin        = -1.0

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

    diskedge.T          = -1.0
    diskedge.Tspot      = -1.0
    diskedge.ZetaMid    = -1.0
    diskedge.ZetaWidth  = -1.0

    diskrim.type = "MISSING"
    diskrim.awidth      = -1.0
    diskrim.Hmax        = -1.0
    diskrim.Hmin        = -1.0
    diskrim.ZetaHmax    = -1.0
    diskrim.Tmax        = -1.0
    diskrim.Tmin        = -1.0
    diskrim.ZetaTmax    = -1.0
    diskrim.points      =  0

    disktorus.type = "MISSING"
    disktorus.azero     = -1.0
    disktorus.awidth    = -1.0
    disktorus.Hmax      = -1.0
    disktorus.Hmin      = -1.0
    disktorus.ZetaHmax  = -1.0
    disktorus.Tmax      = -1.0
    disktorus.Tmin      = -1.0
    disktorus.ZetaTmax  = -1.0
    disktorus.points    =  0

    innerdisk.T         = -1.0
    innerdisk.L         = -1.0
    
    diskspot.nspots     = 0

    adc.L               = -1.0
    adc.height          = -1.0

    thirdlight.nbands   = 0
    thirdlight.orbphase = -100.0

    Grid.Nxtiles        = -1
    Grid.Nztiles        = -1

    data.nbands         = 0

    return

def CalcSysPars():
    """    
    /****************************************************

    This function uses the various input parameters to calculate
    other useful system parameters.  It also converts some
    other input data to cgs units.

    *********************************************************/
    """
    if control.diagnostics != "OFF":
        if( (control.diagnosephase < orbit.phasemin)
	         or (control.diagnosephase > orbit.phasemax) ):
            Quit("diagnosephase must be ge phasemin and le phasemax.")
        control.diagnoseindex = 0.5 + ( control.diagnosephase - orbit.phasemin ) / orbit.deltaphase

    orbit.maxpindex    = 1.0e-7 + (orbit.phasemax - orbit.phasemin) / orbit.deltaphase
    for band in range(1, orbit.nbands):
        orbit.minlambda[band] *= 1.0e-8
        orbit.maxlambda[band] *= 1.0e-8
    if( orbit.normalize == "FITDATA"):
        if( orbit.normfilter == "SQUARE"):
            orbit.normMinlambda *= 1.0e-8
            orbit.normMaxlambda *= 1.0e-8
    syspars.p          *= 86400.0;
    syspars.omega      = math.pi*2 / syspars.p
    syspars.i          *= math.pi*2 / 360.0
    if( syspars.K2 > 0.0 ):
        syspars.K2         *= 1.0e5;
        V2orb              = syspars.K2 / math.sin( syspars.i )
        syspars.a          = V2orb * (1.0 + syspars.q) / syspars.omega
        M1plusM2           = pow( syspars.omega, 2) * pow( syspars.a, 3) / G
        syspars.M1         = M1plusM2 / (1.0 + syspars.q)
        syspars.M2         = syspars.M1 * syspars.q
    else:
        syspars.M1      = syspars.M1 * MSOL
        syspars.M2      = syspars.M1 * syspars.q
        M1plusM2        = syspars.M1 + syspars.M2
        a3              = (G * M1plusM2 * syspars.p * syspars.p) / (4.0 * math.pi * math.pi )
        syspars.a       = pow( a3, 0.3333333333333 )
        a_2             = syspars.a / (1.0 + syspars.q)
        V2orb           = syspars.omega * a_2 / 1.0e5
        syspars.K2      = V2orb * math.sin( syspars.i )
    syspars.zcm        = syspars.a / (1.0 + syspars.q)
    syspars.rL1        = FindL1()  
    syspars.VL1        = V( syspars.rL1, 0.0, 0.0 )
    syspars.MeanLobe2Radius = MeanRocheRadius( syspars.q )
    syspars.MeanLobe1Radius = MeanRocheRadius( 1.0 / syspars.q )

    star1.sigmaT4      = SIGMA * pow( star1.T, 4.0 )
    r2                 = star1.L / ( 4.0 * math.pi * star1.sigmaT4 )
    star1.radius       = math.sqrt( r2 )

    if( star2spot.nspots > 0 ):
        for i in range (1, star2spot.nspots):
            star2spot.theta[i]  *= 2*math.pi / 360.0
            star2spot.phi[i]    *= 2*math.pi / 360.0
            star2spot.radius[i] *= 2*math.pi / 360.0
    disk.zetazero      *= math.pi*2 / 360.0
    maindisk.amin      *= syspars.a
    maindisk.amax      *= syspars.a
    maindisk.Hmax      *= syspars.a

    diskedge.ZetaMid   *= 2*math.pi / 360.0
    diskedge.ZetaWidth *= 2*math.pi / 360.0

    if( control.diskrim == "ON"):
        if( diskrim.type == "POINT"):
            MakeDiskRim()
            for i in range(1, diskrim.points):
                diskrim.PointZeta[i] *= 2*math.pi / 360.0
                diskrim.PointH[i]    *= syspars.a
        diskrim.awidth     *= syspars.a
        diskrim.Hmax       *= syspars.a
        diskrim.Hmin       *= syspars.a
        diskrim.ZetaHmax   *= math.pi*2 / 360.0
        diskrim.ZetaTmax   *= math.pi*2 / 360.0

    if( control.disktorus == "ON"):
        if( disktorus.type == "POINT"):
            MakeDiskTorus()
            for i in range(1, disktorus.points):
                disktorus.PointZeta[i] *= 2*math.pi / 360.0
                disktorus.PointH[i]    *= syspars.a
        disktorus.azero    *= syspars.a
        disktorus.awidth   *= syspars.a
        disktorus.Hmax     *= syspars.a
        disktorus.Hmin     *= syspars.a
        disktorus.ZetaHmax *= 2*math.pi / 360.0
        disktorus.ZetaTmax *= 2*math.pi / 360.0

    if( control.innerdisk == "ON"):
        innerdisk.sigmaT4      = SIGMA * pow( innerdisk.T, 4.0 )
        r2                     = innerdisk.L / ( 2.0 * math.pi * innerdisk.sigmaT4 )
        innerdisk.radius       = math.sqrt( r2 )

    if( diskspot.nspots > 0 ):
      for i in range(1, diskspot.nspots):
          if( diskspot.zetamin[i] > diskspot.zetamax[i] ):
              diskspot.nspots += 1
              if( diskspot.nspots >= 20 ):
	             Quit("diskspot.nspots too large in function CalcSysPars.")
              for( j = diskspot.nspots; j > i; j-- ):
	             diskspot.zetamin[j] = diskspot.zetamin[j-1]
                   diskspot.zetamax[j] = diskspot.zetamax[j-1]
                   diskspot.amin[j] = diskspot.amin[j-1]
                   diskspot.amax[j] = diskspot.amax[j-1]
                   diskspot.spotToverT[j] = diskspot.spotToverT[j-1]
               
              diskspot.zetamax[i+1] = 360.0
	        diskspot.zetamin[i+1] = diskspot.zetamin[i]
 	        diskspot.zetamax[i] = diskspot.zetamax[i]
              diskspot.zetamin[i] = 0.0
          
          diskspot.zetamin[i] *= math.pi*2 / 360.0
          diskspot.zetamax[i] *= math.pi*2 / 360.0
          diskspot.amin[i]    *= syspars.a
          if( diskspot.amin[i] < maindisk.amin ):
              diskspot.amin[i] = maindisk.amin
          diskspot.amax[i]    *= syspars.a
         if( diskspot.amax[i] > maindisk.amax ):
              diskspot.amax[i] = maindisk.amax

    if( control.adc == "ON"):
        adc.height *= syspars.a

    if( control.thirdlight == "ON"):
        for( band in range(1, thirdlight.nbands):
            thirdlight.minlambda[band] *= 1.0e-8
            thirdlight.maxlambda[band] *= 1.0e-8

   for band in range(1, data.nbands):
       data.minlambda[band] *= 1.0e-8
       data.maxlambda[band] *= 1.0e-8

   if( control.diagnostics == "INSPECTSYSPARS"):
       WriteSysPars()
       Quit("Quit in main.c after INSPECTSYSPARS.")
   return
