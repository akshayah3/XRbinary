# -*- coding: utf-8 -*-
"""
Functions concerned with the geometry of the accretion disk 
around the compact star are in this file.
"""
import math
from diskflux import Disk


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
    maindiskarea = math.pi * ( Disk.amax * Disk.amax
                              - Disk.amin * Disk.amin)
    rimheight= 0.0
    if( control.diskrim == "ON"):
        rimheight = 0.5 * (diskrim.Hmin + diskrim.Hmax)
    edgearea = 2.0 * math.pi * Disk.amax * (Disk.Hmax + rimheight)
    extraArim = 0.0
    extraRrim = 0.0
    if( control.diskrim == "ON"):
        flatArea = math.pi * Disk.amax * Disk.amax - math.pi * (Disk.amax - diskrim.awidth) * (Disk.amax - diskrim.awidth)
        extra = 0.5 * (rimheight * rimheight) / (diskrim.awidth * diskrim.awidth)
        extraArim = flatArea * extra;
        extraRrim = diskrim.awidth * extra
    extraAtorus = 0.0
    extraRtorus = 0.0
    if( control.disktorus == "ON" ):
       flatArea =   math.pi * (disktorus.azero + 0.5 * disktorus.awidth) * (disktorus.azero + 0.5 * disktorus.awidth)- math.pi * (disktorus.azero - 0.5 * disktorus.awidth) * (disktorus.azero - 0.5 * disktorus.awidth)
       torusheight = 0.5 * (disktorus.Hmax + disktorus.Hmin)
       extra =  0.5 * (torusheight * torusheight) / (0.5 * disktorus.awidth * 0.5 * disktorus.awidth);
       extraAtorus = flatArea * extra
       extraRtorus = disktorus.awidth * extra
    rougharea = 2.0 * ( maindiskarea +  edgearea + extraAtorus + extraArim)

    targetarea = rougharea / disk.targetNtiles
    targetda = math.sqrt( targetarea )
    nringMain = 0.5 + (Disk.amax - Disk.amin
                              + extraRtorus + extraRrim) / targetda
    if( nringMain < 2 ):
        Quit("disk.targetNtiles too small to cover the disk properly.")
    daMain = (Disk.amax - Disk.amin - 1.0e-7*Disk.amax) / nringMain
    Disk.da = daMain
    nringEdge = 0.5 + (Disk.Hmax + rimheight) / targetda
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
        a1 = (1.0 + 1.0e-7) * Disk.amin + (ring -1) * daMain
        a2 = a1 + daMain
        h1 = MainDiskH( a1 )
        h2 = MainDiskH( a2 )
        a = 0.5 * (a1 + a2)
        if( control.disktorus == "ON"):
	      if( (a > (disktorus.azero - 0.5 * disktorus.awidth))
	              and (a < (disktorus.azero + 0.5 * disktorus.awidth)) ):
                h1 = DiskTopH( a1, disktorus.ZetaHmax )
                h2 = DiskTopH( a2, disktorus.ZetaHmax )
        if( control.diskrim == "ON"):
	       if( a > (Disk.amax - diskrim.awidth) ):
                 h1 = DiskTopH( a1, diskrim.ZetaHmax )
                 h2 = DiskTopH( a2, diskrim.ZetaHmax )
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
	          Quit("Attempted to make too many disk tiles.")
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
            TDiskz[tilenumber] = TDiskRho[tilenumber] * math.cos( zeta ) + syspars.a
            TDisknormCyl[tilenumber] = DiskNormal( "TOP", a1, a2, zeta1, zeta2)
            TDisknormCart[tilenumber] = Cyl2Cart( TDisknormCyl[tilenumber], 
                                          zeta )
            TDiskdS[tilenumber] = TopTileArea( tilenumber, a1, a2, zeta1, zeta2)
            TDiskT[tilenumber] = DiskTopT( a, zeta )
            TDiskT4[tilenumber] = math.pow( TDiskT[tilenumber], 4.0)
        if( ring == 1 ):
            Disk.Tamin = Disk.MainDiskT( a, 0.0)
        if( ring == nringMain ):
            Disk.Tamax = Disk.MainDiskT( a, 0.0)
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
        height = Disk.Hmax + DiskRimH( Disk.amax, zeta )
        dh = height / nringEdge
        for ring in range(1, nringEdge):
            tilenumber += 1
            if( tilenumber >= (MAXDISKTILES / 2) ):
	           Quit("Attempted to make too many disk tiles.")
            TDiska[tilenumber] = maindisk.amax
            h = (ring - 0.5) * dh
            TDiskRho[tilenumber] = AToRho( maindisk.amax, zeta )
            TDiskZeta[tilenumber] = zeta
            TDiskH[tilenumber] = h
            TDiskx[tilenumber] = TDiskRho[tilenumber] * math.sin( zeta )
            TDisky[tilenumber] = TDiskH[tilenumber]
            TDiskz[tilenumber] = TDiskRho[tilenumber] * math.cos( zeta ) + syspars.a
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
	      Quit("Attempted to make too many disk tiles.")
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
    disk.Ntiles = 2 * tilenumber

    """   
    /************************************************************
   
      The last step:  Calculate the specific intensity emitted
      by each disk tile.  The disk is assumed to emit angle-independent
      black body intensities.  The mean intensities depend, then,
      on the bandpass:
   
    ************************************************************/
    """
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            for i in range(1, disk.Ntiles):
                TDiskI[band][i] = BBSquareIntensity( TDiskT[i], 
                                  orbit.minlambda[band],
                                  orbit.maxlambda[band])
        else:
            for i in range(1, disk.Ntiles):
                TDiskI[band][i] = BBFilterIntensity( 
                                  TDiskT[i], orbit.filter[band])

    if( control.diagnostics == "INSPECTDISKTILES"):
        InspectDiskTiles( targetarea, nringMain, maintiles, 
                            nringEdge, edgetiles)
        Quit("Quit after INSPECTDISKTILES.")

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
        if( a < Disk.amin ):
            height = 0.0
        elif ( a > ( (1.0 + 1.0e-7) * Disk.amax ) ):
            height = 0.0
        else:
            x = (a - Disk.amin) / (Disk.amax - Disk.amin)
            if( x == 0.0 ):
	           height = Disk.Hmax
            else:
                 height = Disk.Hmax * math.pow( x, Disk.Hpow )
        return height

    