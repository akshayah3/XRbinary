# -*- coding: utf-8 -*-
"""
Functions concerned with the geometry of the accretion disk 
around the compact star are in this file.
"""
import sys
import math
from .diskflux import maindisk
from .parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, XYGrid, dataparams, ReadInput


def MakeDiskTiles():
    """
    This function distributes tiles over the surface of the disk
    and calculates various properties of the tiles. The tiles
    all have roughly the same area.
    """
    if( verbose == "ON"):
     print(" Begin making disk tiles.\n")
    
    """
    /************************************************************
    
       Calculate the typical area and dimensions of the tiles.
    
       Method:  Calculate a rough area for the disk and divide it
           by the target number of tiles to get a target area for
           the tiles.  The target edge length of the tiles is 
               targetda  =  sqrt( target area).
    
    **************************************************************/
    """
    maindiskarea = math.pi * ( maindisk.amax * maindisk.amax
                              - maindisk.amin * maindisk.amin)
    rimheight= 0.0
    if( flowcontrol.diskrim == "ON"):
        rimheight = 0.5 * (diskrimpars.Hmin + diskrimpars.Hmax)
    edgearea = 2.0 * math.pi * maindisk.amax * (maindisk.Hmax + rimheight)
    extraArim = 0.0
    extraRrim = 0.0
    if( flowcontrol.diskrim == "ON"):
        flatArea = math.pi * maindisk.amax * maindisk.amax - math.pi * (maindisk.amax - diskrimpars.awidth) * (maindisk.amax - diskrimpars.awidth)
        extra = 0.5 * (rimheight * rimheight) / (diskrimpars.awidth * diskrimpars.awidth)
        extraArim = flatArea * extra;
        extraRrim = diskrimpars.awidth * extra
    extraAtorus = 0.0
    extraRtorus = 0.0
    if( flowcontrol.disktorus == "ON" ):
       flatArea =   math.pi * (disktorusparams.azero + 0.5 * disktorusparams.awidth) * (disktorusparams.azero + 0.5 * disktorusparams.awidth)- math.pi * (disktorusparams.azero - 0.5 * disktorusparams.awidth) * (disktorusparams.azero - 0.5 * disktorusparams.awidth)
       torusheight = 0.5 * (disktorusparams.Hmax + disktorusparams.Hmin)
       extra =  0.5 * (torusheight * torusheight) / (0.5 * disktorusparams.awidth * 0.5 * disktorusparams.awidth)
       extraAtorus = flatArea * extra
       extraRtorus = disktorusparams.awidth * extra
    rougharea = 2.0 * ( maindiskarea +  edgearea + extraAtorus + extraArim)

    targetarea = rougharea / wholediskpars.targetNtiles
    targetda = math.sqrt( targetarea )
    nringMain = 0.5 + (maindisk.amax - maindisk.amin
                              + extraRtorus + extraRrim) / targetda
    if( nringMain < 2 ):
        sys.exit("disk.targetNtiles too small to cover the disk properly.")
    daMain = (maindisk.amax - maindisk.amin - 1.0e-7*maindisk.amax) / nringMain
    maindisk.da = daMain
    nringEdge = 0.5 + (maindisk.Hmax + rimheight) / targetda
    if( nringEdge < 1 ):
        nringEdge = 1
    """
    /******************************************************
   
      Now lay out the tiles.  Note that the tiles are first
      distributed in (a,zeta) space and then converted to
      (rho,zeta) space.
      Every ring must have at least mintiles tiles per ring.
   
    **********************************************************/

    *********************************************************
   
      Tiles on the top surface of the disk.
   
    ********************************************************/
    """
    mintiles = 30
    tilenumber = 0
    for ring in range(1, nringMain):
        a1 = (1.0 + 1.0e-7) * maindisk.amin + (ring -1) * daMain
        a2 = a1 + daMain
        h1 = MainDiskH( a1 )
        h2 = MainDiskH( a2 )
        a = 0.5 * (a1 + a2)
        if( flowcontrol.disktorus == "ON"):
	      if( (a > (disktorusparams.azero - 0.5 * disktorusparams.awidth))
	              and (a < (disktorusparams.azero + 0.5 * disktorusparams.awidth)) ):
                h1 = DiskTopH( a1, disktorusparams.ZetaHmax )
                h2 = DiskTopH( a2, disktorusparams.ZetaHmax )
        if( flowcontrol.diskrim == "ON"):
	       if( a > (maindisk.amax - diskrimpars.awidth) ):
                 h1 = DiskTopH( a1, diskrimpars.ZetaHmax )
                 h2 = DiskTopH( a2, diskrimpars.ZetaHmax )
        slope = (h2 - h1) / (a2 - a1)
        projection = math.sqrt( 1.0 + slope*slope )
        if( projection > 3.0 ): 
	      projection = 4.0
        ringarea = projection * ( math.pi*a2*a2 - math.pi*a1*a1)
        nzeta = 0.5 + (ringarea / targetarea)
        if( nzeta < mintiles ):
            nzeta = mintiles
        dzeta = 2*math.pi / nzeta
        for i in range(1, nzeta):
            tilenumber += 1
            if( tilenumber >= (MAXDISKTILES / 2) ):
	          sys.exit("Attempted to make too many disk tiles.")
            zeta1 = (i - 1) * dzeta
            zeta2 = zeta1 + dzeta - 1.0e-10
            zeta = 0.5 * (zeta1 + zeta2)
            TDiska[tilenumber] = a
            TDiskZeta[tilenumber] = zeta
            TDiskH[tilenumber] = 0.5 * ( DiskTopH( a1, zeta) 
                                   +  DiskTopH( a2, zeta) )
            TDiskRho[tilenumber] = AToRho( a, zeta )
            TDiskx[tilenumber] = TDiskRho[tilenumber] * math.sin( zeta )
            TDisky[tilenumber] = TDiskH[tilenumber]
            TDiskz[tilenumber] = TDiskRho[tilenumber] * math.cos( zeta ) + systemparams.a
            TDisknormCyl[tilenumber] = DiskNormal( "TOP", a1, a2, zeta1, zeta2)
            TDisknormCart[tilenumber] = Cyl2Cart( TDisknormCyl[tilenumber], 
                                          zeta )
            TDiskdS[tilenumber] = TopTileArea( tilenumber, a1, a2, zeta1, zeta2)
            TDiskT[tilenumber] = DiskTopT( a, zeta )
            TDiskT4[tilenumber] = math.pow( TDiskT[tilenumber], 4.0)
        if( ring == 1 ):
            maindisk.Tamin = maindisk.MainDiskT( a, 0.0)
        if( ring == nringMain ):
            maindisk.Tamax = maindisk.MainDiskT( a, 0.0)
    maintiles = tilenumber


    """
    /**************************************************
   
      Now put the tiles on the edge of the disk.  The tiles
      are in a different order than on the surface.
   
    **************************************************/
    """
    nzeta = 0.5 + edgearea / (nringEdge * targetarea)
    if( nzeta < mintiles ):
        nzeta = mintiles
    dzeta = math.pi*2 / nzeta
    for i in range(1, nzeta):
        zeta1 = (i - 1) * dzeta
        zeta2 = zeta1 + dzeta - 1.0e-10;
        zeta = 0.5 * (zeta1 + zeta2)
        height = maindisk.Hmax + DiskRimH( maindisk.amax, zeta )
        dh = height / nringEdge
        for ring in range(1, nringEdge):
            tilenumber += 1
            if( tilenumber >= (MAXDISKTILES / 2) ):
	           sys.quit("Attempted to make too many disk tiles.")
            TDiska[tilenumber] = maindisk.amax
            h = (ring - 0.5) * dh
            TDiskRho[tilenumber] = AToRho( maindisk.amax, zeta )
            TDiskZeta[tilenumber] = zeta
            TDiskH[tilenumber] = h
            TDiskx[tilenumber] = TDiskRho[tilenumber] * math.sin( zeta )
            TDisky[tilenumber] = TDiskH[tilenumber]
            TDiskz[tilenumber] = TDiskRho[tilenumber] * math.cos( zeta ) + systemparams.a
            TDisknormCyl[tilenumber] = DiskNormal( "EDGE", a1, a2, zeta1, zeta2)
            TDisknormCart[tilenumber] = Cyl2Cart( TDisknormCyl[tilenumber],
                                                      zeta )
            TDiskdS[tilenumber] = EdgeTileArea( tilenumber, zeta1, zeta2, dh)
            TDiskT[tilenumber] = DiskEdgeT( zeta )
            TDiskT4[tilenumber] = math.pow( TDiskT[tilenumber], 4.0)
    edgetiles = tilenumber - maintiles
   
    """
    /************************************************************
   
      Duplicate the tiles on the bottom of the disk.  This is 
      necessary because the bottom tiles can irradiate and heat 
      parts of the secondary star that are visible.
   
    ************************************************************/
    """
    for i in range(1, tilenumber):
        if( (i+tilenumber) >= MAXDISKTILES ):
	      sys.quit("Attempted to make too many disk tiles.")
        TDiska[i+tilenumber]    =  TDiska[i]
        TDiskRho[i+tilenumber]  =  TDiskRho[i]
        TDiskZeta[i+tilenumber] =  TDiskZeta[i]
        TDiskH[i+tilenumber]    = -TDiskH[i]
        TDiskx[i+tilenumber]    =  TDiskx[i]
        TDisky[i+tilenumber]    = -TDisky[i]
        TDiskz[i+tilenumber]    =  TDiskz[i]
        TDiskdS[i+tilenumber]   =  TDiskdS[i]
        TDisknormCyl[i+tilenumber].rho  =  TDisknormCyl[i].rho
        TDisknormCyl[i+tilenumber].zeta =  TDisknormCyl[i].zeta
        TDisknormCyl[i+tilenumber].h    = -TDisknormCyl[i].h
        TDisknormCart[i+tilenumber].x   =  TDisknormCart[i].x
        TDisknormCart[i+tilenumber].y   = -TDisknormCart[i].y
        TDisknormCart[i+tilenumber].z   =  TDisknormCart[i].z
        TDiskT[i+tilenumber]  = TDiskT[i]
        TDiskT4[i+tilenumber] = TDiskT4[i]
    wholediskpars.Ntiles = 2 * tilenumber

    """   
    /************************************************************
   
      The last step:  Calculate the specific intensity emitted
      by each disk tile.  The disk is assumed to emit angle-independent
      black body intensities.  The mean intensities depend, then,
      on the bandpass:
   
    ************************************************************/
    """
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            for i in range(1, wholediskpars.Ntiles):
                TDiskI[band][i] = BBSquareIntensity( TDiskT[i], 
                                  orbitparams.minlambda[band],
                                  orbitparams.maxlambda[band])
        else:
            for i in range(1, wholediskpars.Ntiles):
                TDiskI[band][i] = BBFilterIntensity( 
                                  TDiskT[i], orbitparams.filter[band])

    if( flowcontrol.diagnostics == "INSPECTDISKTILES"):
        InspectDiskTiles( targetarea, nringMain, maintiles, 
                            nringEdge, edgetiles)
        sys.exit("Quit after INSPECTDISKTILES.")

    return

    def DiskTopH(a, zeta):
        """        
        This function returns the upper boundary of the disk
        at position (rho, zeta).
        """
        hdisk = MainDiskH( a )
        hrim  = DiskRimH( a, zeta )
        htorus = DiskTorusH( a, zeta )
        height = hdisk + hrim + htorus

        return height

    def DiskBottomH(a, zeta):
        """
        This function returns the lower boundary of the disk at position (rho, zeta).
        The current version assumes that the disk is symmetric
        about the orbital plane.
        """
        hdisk = MainDiskH( a )
        hrim  = DiskRimH( a, zeta )
        htorus = DiskTorusH( a, zeta )
        height = -( hdisk + hrim + htorus )

        return height

    def MainDiskH(a):
        """
        This function returns the height of the main, axi-symmetric
        part of the disk.  The height has the functional form
        h = Hmax * ( (a - amin) / (amax - amin) )^Hpow
        for
        amin < a < amax
        and
        h = 0 otherwise.
        Note:  To make a disk with constant thickness equal to
           Hmax, just set Hpow = 0.

        The weird factor (1.0 + 1.0e-7) is protection against 
        roundoff error at the edge of the disk.
        """
        if( a < maindisk.amin ):
            height = 0.0
        elif ( a > ( (1.0 + 1.0e-7) * maindisk.amax ) ):
            height = 0.0
        else:
            x = (a - maindisk.amin) / (maindisk.amax - maindisk.amin)
            if( x == 0.0 ):
	           height = maindisk.Hmax
            else:
                 height = maindisk.Hmax * math.pow( x, maindisk.Hpow )
        return height

def DiskRimH( a, zeta ):
    """
    This function returns the height of the disk rim at 
    position (a, zeta) .  The cross section of the rim is 
    modeled as a truncated ellipse added on top of the main 
    disk with
       vertical semi-major axis   = H
       horizontal semi-minor axis = diskrim.awidth
 
        (h/H)^2 + ( (maindisk.amax -a) / awidth )^2 = 1
 
    The rim height H is a function of zeta but its width is not.
    There are currently two possibilities for the zeta dependence
 
    1) SINUSOID
       H = 0.5 * (Hmax + Hmin)
             + 0.5* (Hmax - Hmin) * cos( zeta - zetaHmax );
 
    2) POINT
       The disk rim height and temperature is defined by a set 
       of points, one point per DISKRIMPARS= line in the parameter file:
          DISKRIMPARS=  POINT   Zeta1   H1   T1
          DISKRIMPARS=  POINT   Zeta2   H2   T2    
             .        .       .     .    .
             .        .       .     .    .
       The Zetas need not be in order.
       The heights are linearly interpolated between the specified points.
 
    The weird factor (1.0 + 1.0e-7) is protection against 
    roundoff error at the edge of the disk.
    """
    if( zeta > 2*math.pi ):
        sys.exit("zeta greater than TWOPI in DiskRimH.")
    if( zeta < 0.0 ):
        sys.exit("zeta less than zero in DiskRimH.")

    if( a < (maindisk.amax - diskrimpars.awidth) ):
        height = 0.0
        return( height )
    if( a > ( (1.0 + 1.0e-7) * maindisk.amax ) ):
        height = 0.0
        return( height )

    if( diskrimpars.type == "SINUSOID" ):
        Hzeta =    0.5 * ( diskrimpars.Hmax + diskrimpars.Hmin ) + 0.5 * ( diskrimpars.Hmax - diskrimpars.Hmin ) * math.cos( zeta - diskrimpars.ZetaHmax )
    elif( diskrimpars.type == "POINT" ):
        if( diskrimpars.points == 1 ):
            Hzeta = diskrimpars.PointH[1]
        else:
	      for i in range(1, diskrimpars.points):
               zetalow = diskrimpars.PointZeta[i]
               Hlow = diskrimpars.PointH[i]
               if( i < diskrimpars.points ):
                   zetahigh = diskrimpars.PointZeta[i+1];
                   Hhigh = diskrimpars.PointH[i+1]
               else:
                   zetahigh = 2*math.pi
                   Hhigh = diskrimpars.PointH[1]
               if( (zeta >= zetalow) and (zeta < zetahigh) ):
                  slope = (Hhigh - Hlow) / (zetahigh - zetalow)
                  Hzeta = Hlow + slope * (zeta - zetalow)
                  break
    else:
        sys.exit("Unrecognized disk rim type in DiskRimH.")

    x = (maindisk.amax - a) / diskrimpars.awidth
    y = math.sqrt( 1.0 - x*x)
    height = Hzeta * y

    return( height )


def DiskTorusH( a, zeta ):
    """
    This function returns the height of the disk torus
    at position (a, zeta).   The cross section of the torus is 
    modeled as an ellipse added on top of the main disk with
       ellipse center             = disktorus.azero
       vertical semi-major axis   = H
       horizontal semi-minor axis = 0.5 * disktorus.awidth
 
        (h/H)^2 + ( (disktorus.azero -a) / awidth )^2 = 1
 
    The torus height H is a function of zeta but its width is not. 
    There are currently two possibilities for the zeta dependence
 
    1) SINUSOID
       H = 0.5 * (Hmax + Hmin)
             + 0.5* (Hmax - Hmin) * cos( zeta - zetaHmax );
 
    2) POINT
       The disk torus height and temperature is defined by a set of
       points, one point per DISKTORUSPARS= line in the parameter file:
          DISKTORUSPARS=  POINT   Zeta1   H1   T1
          DISKTORUSPARS=  POINT   Zeta2   H2   T2    
             .        .       .     .    .
             .        .       .     .    .
       The heights are linearly interpolated between the specified points.
 
       The Zetas must be in increasing order and disktorus.PointZeta[1]
       must be 0 degrees (this avoids messy computer code).
       The heights are linearly interpolated between the specified points.
    """
    if( zeta > 2*math.pi ):
        sys.exit("zeta greater than TWOPI in DiskTorusH.")
    if( zeta < 0.0 ):
        sys.exit("zeta less than zero in DiskTorusH.")

    if( a < (disktorusparams.azero - 0.5 * disktorusparams.awidth) ):
        height = 0.0
        return( height )
    if( a > (disktorusparams.azero + 0.5 * disktorusparams.awidth) ):
        height = 0.0
        return( height )

    if( disktorusparams.type == "SINUSOID" ):
        Hzeta =    0.5 * ( disktorusparams.Hmax + disktorusparams.Hmin ) + 0.5 * ( disktorusparams.Hmax - disktorusparams.Hmin ) * math.cos( zeta - disktorusparams.ZetaHmax )
    elif( disktorusparams.type == "POINT" ):
        if( disktorusparams.points == 1 ):
            Hzeta = disktorusparams.PointH[1]
        else:
	      for i in range(1, disktorusparams.points):
               zetalow = disktorusparams.PointZeta[i]
               Hlow = disktorusparams.PointH[i]
               if( i < disktorusparams.points ):
                   zetahigh = disktorusparams.PointZeta[i+1]
                   Hhigh = disktorusparams.PointH[i+1]
               else:
                   zetahigh = 2*math.pi
                   Hhigh = disktorusparams.PointH[1]
               if( (zeta >= zetalow) and (zeta < zetahigh) ):
                   slope = (Hhigh - Hlow) / (zetahigh - zetalow)
                   Hzeta = Hlow + slope * (zeta - zetalow);
                   break
    else:
        sys.exit("Unrecognized disk torus type in DiskTorusH.");

    x = (disktorusparams.azero - a) / (0.5 * disktorusparams.awidth)
    y = math.sqrt( 1.0 - x*x)
    height = Hzeta * y

    return( height )

def DiskNormal( where, a1, a2, zeta1, zeta2):
    """
    This function calculates the unit vector that is
    normal to the surface of a disk tile.  Note: the vector is
    normal to the tile, not normal to the surface of the disk
    underlying the tile.
    """ 
    if( (a2 - a1) <= 0.0 ):
        sys.exit("a2 le a1 in DiskNormal().")
    if( (zeta2 - zeta1) <= 0.0 ):
        sys.exit("zeta2 le zeta1 in DiskNormal().")
    a = 0.5 * (a1 + a2)
    zeta = 0.5 * (zeta1 + zeta2)
    if( where == "TOP"):
        rho1 = AToRho( a1, zeta)
        rho2 = AToRho( a2, zeta)
        rho =  AToRho( a,  zeta)
        h1 = DiskTopH( a1, zeta)
        h2 = DiskTopH( a2, zeta)
        DhDrho = (h2 - h1) / (rho2 - rho1)
        h1 = DiskTopH( a, zeta1)
        h2 = DiskTopH( a, zeta2)
        DhDzeta = (h2 - h1) / (zeta2 - zeta1)
        norm.rho  = -DhDrho
        norm.zeta = -(1.0 / rho) * DhDzeta
        norm.h    =  1.0
    elif( where == "EDGE"):
        rho1 = AToRho( a, zeta1)
        rho2 = AToRho( a, zeta2)
        rho =  AToRho( a, zeta)
        DrhoDzeta = (rho2 - rho1) / (zeta2 - zeta1)
        norm.rho  = 1.0
        norm.zeta = -(1.0 / rho) * DrhoDzeta
        norm.h    = 0.0
    else:
      sys.exit("Unrecognized location in DiskNormal.")

    beta = 1.0 / math.sqrt( norm.rho*norm.rho + norm.zeta*norm.zeta + norm.h*norm.h)
    norm.rho  = beta * norm.rho
    norm.zeta = beta * norm.zeta
    norm.h    = beta * norm.h

    return( norm )


def TopTileArea( itile, a1, a2, zeta1, zeta2):
    """
    This function calculates the area of tiles on the top surface
    of the disk extending from a1 to a2, and from zeta1 to zeta2.
    Note that the tile is finite in size, so the expression
    for flatarea is not just  d S = rho d rho d zeta.
    """
    if( zeta1 < 0.0 ):
        sys.exit("zeta less than zero in TopTileArea.")
    if( zeta2 > 2*math.pi ):
        sys.exit("zeta greater than TWOPI in TopTileArea.")
    zeta = 0.5 * (zeta1 + zeta2)
    dzeta = zeta2 - zeta1
    if( dzeta < 0.0 ):
        sys.exit("dzeta is less than zero in TopTileArea.")

    rho1 = AToRho( a1, zeta)
    rho2 = AToRho( a2, zeta)
    flatarea = math.pi * (rho2*rho2 - rho1*rho1) * ( dzeta / 2*math.pi )

    if( TDisknormCyl[itile].h == 0.0 ):
        sys.exit("TDisknormCyl.h equals zero in TopTileArea.")
    Area = flatarea / TDisknormCyl[itile].h

    return( Area )

def EdgeTileArea( itile, zeta1, zeta2, dh):
    """
    This function calculates the area of edge tiles.
    """
    if( zeta1 < 0.0 ):
        sys.exit("zeta less than zero in EdgeTileArea.")
    if( zeta2 > 2*math.pi ):
        sys.exit("zeta greater than TWOPI in EdgeTileArea.")
    zeta = 0.5 * (zeta1 + zeta2)
    dzeta = zeta2 - zeta1
    if( dzeta < 0.0 ):
        sys.exit("dzeta less than zero in EdgeTileArea.")

    flatarea = TDiskRho[itile] * dzeta * dh

    if( TDisknormCyl[itile].rho == 0.0 ):
        sys.exit("TDisknormCyl.rho equals zero in EdgeTileArea.")
    Area = flatarea / TDisknormCyl[itile].rho

    return( Area )

def RhoToA( rho, zeta ):
    """
    Calculates a from (rho, zeta).
    """
    if( wholediskpars.e == 0.0 ):
        a = rho
    else:
        a = rho * (1.0 + wholediskpars.e * math.cos(zeta - wholediskpars.zetazero) ) / (1.0 - wholediskpars.e * wholediskpars.e)

    return( a )

def AToRho( a, zeta ):
    """
    Calculates rho from ( a, zeta).
    """
    if( wholediskpars.e == 0.0 ):
        rho = a
    else:
        rho =  ( a * (1.0 - wholediskpars.e * wholediskpars.e) ) / (1.0 + wholediskpars.e * cos(zeta - wholediskpars.zetazero) )
                             
    return( rho )

def MakeDiskRim():
    """
    This function is invoked if the disk rim is ON and the rim
    is specified by a set of discrete points along its rim.  
    It does the following:
       1)  Sort the rim points in order of increasing zeta.
       2)  Check for duplicate points.
       3)  Create a point a zeta = 0 if one does not already
             exist by linearly interpolating in zeta.
       4)  Finds maximum and minimum values for H and T.
 
    Note:  The function assumes that zeta is given in degrees.
    """
    if( diskrimpars.points <= 0 ):
        sys.exit("SortRimPoints:  diskrim.points le 0.")
    elif( diskrimpars.points > 1 ):
        for i in range(1, diskrimpars.points):
            for j in range(i+1, diskrimpars.points):
                if( diskrimpars.PointZeta[i] > diskrimpars.PointZeta[j] ):
                    dummy = diskrimpars.PointZeta[i]
                    diskrimpars.PointZeta[i] = diskrimpars.PointZeta[j]
                    diskrimpars.PointZeta[j] = dummy
                    dummy = diskrimpars.PointH[i]
                    diskrimpars.PointH[i] = diskrimpars.PointH[j]
                    diskrimpars.PointH[j] = dummy
                    dummy = diskrimpars.PointT[i]
                    diskrimpars.PointT[i] = diskrimpars.PointT[j]
                    diskrimpars.PointT[j] = dummy

    if( diskrimpars.points > 1 ):
        for i in range(1, diskrimpars.points):
            if( diskrimpars.PointZeta[i] == diskrimpars.PointZeta[i+1] ):
                sys.exit("DISKRIMPARS: At least two points have the same PointZeta.")
    if( diskrimpars.PointZeta[1] != 0.0 ):
        diskrimpars.points += 1
        for i in range(diskrimpars.points, 2, -1):
             diskrimpars.PointZeta[i] = diskrimpars.PointZeta[i-1]
             diskrimpars.PointH[i]    = diskrimpars.PointH[i-1]
             diskrimpars.PointT[i]    = diskrimpars.PointT[i-1]
        diskrimpars.PointZeta[1] = 0.0
        zeta2 = diskrimpars.PointZeta[2];
        zetaN = 360.0 - diskrimpars.PointZeta[diskrimpars.points];
        diskrimpars.PointH[1] = ( zetaN * diskrimpars.PointH[2] + zeta2 * diskrimpars.PointH[diskrimpars.points]) / (zeta2 + zetaN)
        diskrimpars.PointT[1] = ( zetaN * diskrimpars.PointT[2] + zeta2 * diskrimpars.PointT[diskrimpars.points]) / (zeta2 + zetaN)
    diskrimpars.Hmax = diskrimpars.PointH[1]
    diskrimpars.Hmin = diskrimpars.PointH[1]
    diskrimpars.ZetaHmax = diskrimpars.PointZeta[1]
    diskrimpars.Tmin = diskrimpars.PointT[1]
    diskrimpars.Tmax = diskrimpars.PointT[1]
    diskrimpars.ZetaTmax = diskrimpars.PointZeta[1]
    for i in range(1, diskrimpars.points):
        if( diskrimpars.PointH[i] > diskrimpars.Hmax ):
            diskrimpars.Hmax = diskrimpars.PointH[i]
            diskrimpars.ZetaHmax = diskrimpars.PointZeta[i]
        if( diskrimpars.PointH[i] < diskrimpars.Hmin ):
            diskrimpars.Hmin = diskrimpars.PointH[i]
        if( diskrimpars.PointT[i] > diskrimpars.Tmax ):
            diskrimpars.Tmax = diskrimpars.PointT[i]
            diskrimpars.ZetaTmax = diskrimpars.PointZeta[i]
        if( diskrimpars.PointT[i] < diskrimpars.Tmin ):
            diskrimpars.Tmin = diskrimpars.PointT[i]
    return

def MakeDiskTorus():
    """
    This function is invoked if the disk torus is ON and the torus
    is specified by a set of discrete points around the disk.  
    It does the following:
       1)  Sort the torus points in order of increasing zeta.
       2)  Check for duplicate points.
       3)  Create a point a zeta = 0 if one does not already
             exist by linearly interpolating in zeta.
       4)  Finds maximum and minimum values for H and T.
 
    Note:  The function assumes that zeta is given in degrees.
    """
    if( disktorusparams.points <= 0 ):
        sys.exit("SortTorusPoints:  disktorus.points le 0.")
    elif( disktorusparams.points > 1 ):
        for i in range(1, disktorusparams.points):
            for j in range(i+1, disktorusparams.points):
                if( disktorusparams.PointZeta[i] > disktorusparams.PointZeta[j] ):
                    dummy = disktorusparams.PointZeta[i]
                    disktorusparams.PointZeta[i] = disktorusparams.PointZeta[j]
                    disktorusparams.PointZeta[j] = dummy
                    dummy = disktorusparams.PointH[i]
                    disktorusparams.PointH[i] = disktorusparams.PointH[j]
                    disktorusparams.PointH[j] = dummy
                    dummy = disktorusparams.PointT[i]
                    disktorusparams.PointT[i] = disktorusparams.PointT[j]
                    disktorusparams.PointT[j] = dummy
    if( disktorusparams.points > 1 ):
        for i in range(1, disktorusparams.points):
            if( disktorusparams.PointZeta[i] == disktorusparams.PointZeta[i+1] ):
                sys.exit("DISKTORUSPARS: At least two points have the same PointZeta.")

    if( disktorusparams.PointZeta[1] != 0.0 ):
        disktorusparams.points += 1
        for i in range(disktorusparams.points,2,-1):
            disktorusparams.PointZeta[i] = disktorusparams.PointZeta[i-1]
            disktorusparams.PointH[i]    = disktorusparams.PointH[i-1]
            disktorusparams.PointT[i]    = disktorusparams.PointT[i-1]
        disktorusparams.PointZeta[1] = 0.0
        zeta2 = disktorusparams.PointZeta[2]
        zetaN = 360.0 - disktorusparams.PointZeta[disktorusparams.points]
        disktorusparams.PointH[1] = ( zetaN * disktorusparams.PointH[2] + zeta2 * disktorusparams.PointH[disktorusparams.points]) / (zeta2 + zetaN)
        disktorusparams.PointT[1] = ( zetaN * disktorusparams.PointT[2] + zeta2 * disktorusparams.PointT[disktorusparams.points]) / (zeta2 + zetaN)
    disktorusparams.Hmax = disktorusparams.PointH[1]
    disktorusparams.Hmin = disktorusparams.PointH[1]
    disktorusparams.ZetaHmax = disktorusparams.PointZeta[1]
    disktorusparams.Tmin = disktorusparams.PointT[1]
    disktorusparams.Tmax = disktorusparams.PointT[1]
    disktorusparams.ZetaTmax = disktorusparams.PointZeta[1]
    for i in range(1, disktorusparams.points):
        if( disktorusparams.PointH[i] > disktorusparams.Hmax ):
            disktorusparams.Hmax = disktorusparams.PointH[i]
            disktorusparams.ZetaHmax = disktorusparams.PointZeta[i]
        if( disktorusparams.PointH[i] < disktorusparams.Hmin ):
            disktorusparams.Hmin = disktorusparams.PointH[i]
        if( disktorusparams.PointT[i] > disktorusparams.Tmax ):
            disktorusparams.Tmax = disktorusparams.PointT[i]
            disktorusparams.ZetaTmax = disktorusparams.PointZeta[i]
        if( disktorusparams.PointT[i] < disktorusparams.Tmin ):
            disktorusparams.Tmin = disktorusparams.PointT[i]
    return
    