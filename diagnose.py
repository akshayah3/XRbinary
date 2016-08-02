# -*- coding: utf-8 -*-
"""
This file contains functions used to diagnose and test the
line broadening program.
"""
import math 

def InspectInput():
    """
    This function writes out the various input files to files
    called "inputfilename.inspect
    """
    WritePars()
    WriteGDTable()
    WriteLDTable()
    WriteIperpTable()
    WriteIBBfilterTable()
    WriteZzetaTable()
    if( data.nbands > 0):
        for band in range(1, data.nbands):
            WriteData( band )

    return



def WritePars():
    """
    This function writes out the parameters read from parfile.dat
    into a file called parfile.inspect
    """
    filename  = "parfile.inspect"
    out = open(filename, "w")
    if out == None:
        Quit("Cannot open file parfile.inspect.")
    out.write( "DIAGNOSTICS=       %s  %5.3f  %s\n"% control.diagnostics,
                                                      control.diagnosephase,
                                                      control.diagnoseband)

    out.write( "\n")
    out.write( "STAR1=             %s\n"% control.star1)
    out.write( "STAR2=             %s\n"% control.star2)
    out.write( "DISK=              %s\n"% control.disk)
    out.write( "DISKRIM=           %s\n"% control.diskrim)
    out.write( "DISKTORUS=         %s\n"% control.disktorus)
    out.write( "INNERDISK=         %s\n"% control.innerdisk)
    out.write( "DISKSPOTS=         %s\n"% control.diskspots)
    out.write( "ADC=               %s\n"% control.adc)
    out.write( "THIRDLIGHT=        %s\n"% control.thirdlight)
    out.write( "IRRADIATION=       %s\n"% control.irradiation)

    out.write( "\n")
    out.write( "PHASES=            %6.3f  %6.3f  %6.3f\n"% 
	                   orbit.phasemin, orbit.phasemax, orbit.deltaphase)
    out.write( "PHASEOFFSET=      %7.4f\n"%  orbit.phaseoffset)
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"): 
            out.write( "BANDPASS=          SQUARE  %6.0f  %6.0f\n"%
                   orbit.minlambda[band], orbit.maxlambda[band])
        else:
            out.write( "BANDPASS=          FILTER   %s\n"% orbit.filter[band])
    if( orbit.normalize == "OFF"):
        out.write("NORMALIZE=         %s\n"% orbit.normalize)
    else:
        if( orbit.normalize == "MAXVALUE"):
            out.write( "NORMALIZE=   %s  %12.4e\n"% orbit.normalize,
	                                         orbit.normvalue)
        if( orbit.normalize == "FITDATA"):
            if( orbit.normfilter == "SQUARE"):
                out.write( "NORMALIZE=   %s  %s  %7.1f  %7.1f\n"%
		             orbit.normalize, orbit.normfilter, 
                             orbit.normMinlambda, orbit.normMaxlambda)
            else:
	          out.write( "NORMALIZE=   %s  %s\n"%
		   orbit.normalize, orbit.normfilter)

    out.write( "\n")
    out.write("PERIOD=            %10.8f\n"% syspars.p)
    out.write( "K2=                %6.2f\n"%  syspars.K2)
    out.write( "MASSRATIO=         %5.3f\n"%  syspars.q)
    out.write( "INCLINATION=       %5.2f\n"%  syspars.i)

    if( control.star1 == "ON"):
        out.write( "\n")
        out.write( "STAR1LUM=          %9.3e\n"% star1.L)
        out.write( "STAR1TEMP=         %9.3e\n"% star1.T)

    if( control.star2 == "ON"):
        out.write( "\n")
        out.write( "STAR2TILES=       %5ld\n"%  star2.targetNtiles)
        out.write( "STAR2TEMP=        %5.0f\n"% star2.meanT)
        out.write( "STAR2ALBEDO=      %5.2f\n"% star2.albedo)

    if( control.disk == "ON"):
        out.write( "\n")
        out.write( "DISKTILES=         %5ld\n"%   disk.targetNtiles)
        out.write( "DISKE=             %5.3f\n"%  disk.e)
        out.write( "DISKZETAZERO=     %5.1f\n"%  disk.zetazero)
        out.write( "DISKALBEDO=        %5.2f\n"% diskalbedo)

        out.write("\n")
        out.write( "MAINDISKA=         %5.3f  %5.3f\n"% maindisk.amin,
	                                              maindisk.amax)
        out.write( "MAINDISKH=         %5.3f  %4.1f\n"% maindisk.Hmax,
	                                               maindisk.Hpow)
        out.write( "MAINDISKT=        %5.0f  %4.1f\n"% maindisk.Tamax,
	                                              maindisk.Tpow)
        out.write( "\n")
        out.write( "DISKEDGET=        %5.0f  %5.0f  %5.1f  %5.1f\n"%
                                      diskedge.T, diskedge.Tspot, 
                                      diskedge.ZetaMid, diskedge.ZetaWidth)

    if( control.diskrim == "ON"):
        out.write( "\n")
        out.write( "DISKRIMAWIDTH=     %5.3f\n"%      diskrim.awidth)
        if( diskrim.type == "SINUSOID"):
            out.write( "DISKRIMPARS=       %s   %5.3f %5.3f %5.1f   %6.0f %6.0f %5.1f\n"%
	                  diskrim.type, 
                          diskrim.Hmax, diskrim.Hmin, diskrim.ZetaHmax,
                          diskrim.Tmax, diskrim.Tmin, diskrim.ZetaTmax)
        if( diskrim.type == "POINT"):
            for i in range(1, diskrim.points):
                out.write( "DISKRIMPARS=       %s  %5.1f  %5.3f  %7.1f\n"%
		     diskrim.type, diskrim.PointZeta[i],
                         diskrim.PointH[i], diskrim.PointT[i])

    if( control.disktorus == "ON"):
        out.write( "\n")
        out.write( "DISKTORUSAZERO=    %5.3f\n"%    disktorus.azero)
        out.write( "DISKTORUSAWIDTH=   %5.3f\n"%    disktorus.awidth)
        if( disktorus.type == "SINUSOID"):
            out.write( "DISKTORUSPARS=     %s   %5.3f %5.3f %5.1f   %6.0f %6.0f %5.1f\n"%
	                  disktorus.type, 
                          disktorus.Hmax, disktorus.Hmin, disktorus.ZetaHmax,
                          disktorus.Tmax, disktorus.Tmin, disktorus.ZetaTmax)
        if( disktorus.type == "POINT"):
            for i in range(1, disktorus.points):
                out.write( "DISKTORUSPARS=     %s  %5.1f  %5.3f  %7.1f\n"%
		     disktorus.type, disktorus.PointZeta[i],
                         disktorus.PointH[i], disktorus.PointT[i])

    if( control.diskspots == "ON"):
        out.write( "\n")
        for i in range(1, diskspot.nspots):
            out.write( "DISKSPOT=    %5.1f  %5.1f  %5.3f  %5.3f %6.3f\n"%
	             diskspot.zetamin[i], diskspot.zetamax[i],
                     diskspot.amin[i], diskspot.amax[i], diskspot.spotToverT[i])

    if( control.innerdisk == "ON"):
        out.write( "\n")
        out.write( "INNERDISKT=       %9.2e\n"% innerdisk.T)
        out.write( "INNERDISKL=       %9.2e\n"% innerdisk.L)

    if( control.adc == "ON"):
        out.write( "\n")
        out.write( "ADCL=               =  %10.3e\n"% adc.L)
        out.write( "ADCHEIGHT=          =  %11.4e\n"% adc.height)

    if( control.thirdlight == "ON"):
        out.write( "\n")
        out.write( "3rdLIGHTPHASE =   %6.3f\n"%  thirdlight.orbphase)
        for i in range(1, thirdlight.nbands):
            if( thirdlight.filter[i] == "SQUARE"):
                out.write( "3rdLIGHTFRACTION=  %s %6.0f %6.0f  %5.3f\n"%
		        thirdlight.filter[i], thirdlight.minlambda[i],
                        thirdlight.maxlambda[i], thirdlight.fraction[i] )
            else:
	          out.write( "3rdLIGHTFRACTION=  FILTER    %s           %5.3f\n"%
		        thirdlight.filter[i], thirdlight.fraction[i] )

    if( data.nbands > 0):
        out.write( "\n")
        for i in range(1, data.nbands):
            if( data.filter[i] == "SQUARE"):
                out.write( "READDATA=  %s %6.0f %6.0f   %s\n"%
		        data.filter[i], data.minlambda[i],
                        data.maxlambda[i], data.filename[i] )
            else:
	          out.write( "READDATA=  FILTER     %s           %s\n"%
		        data.filter[i], data.filename[i] )

    out.write( "\nEND\n")

    out.close()
    return
    
def WriteGDTable():
    """
    This function writes out the gravity darkening table read 
    from GDTable.dat into a file called GDTable.inspect.
    """
    filename = "GDTable.inspect"
    out = open(filename, "w")
    if out == None: 
        Quit("Cannot open file GDTable.inspect.")

    for i in range(0, maxGDindex):
        out.write("  %5.0f  %5.3f\n"% GDT[i], fourbeta[i])

    out.close()

    return


def WriteLDTable():
    """
    Write the limb darkening table to the file
    LDTable.inspect
    """

    outfile = "LDTable.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file LDTable.inspect.")

    delta = LDT[1] - LDT[0]
    out.write( "        %5.0f  %5.0f   %4.0f\n"% 
                             LDT[0], LDT[maxLDTindex], delta)
    delta = LDlogg[1] - LDlogg[0]
    out.write( "        %4.1f  %4.1f  %4.1f\n"%
	                     LDlogg[0], LDlogg[maxLDgindex], delta)
    nfilters = maxLDfilterindex + 1
    outputline = "     %2ld  "% (nfilters)
    for findex in range(0, maxLDfilterindex):
        outputline +=  " "
        outputline += LDfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, maxLDTindex):
        for gindex in range(0, maxLDgindex):
            for findex in range(0, maxLDfilterindex):
                out.write( " %5.1f  %5.0f  %s   %7.4f   %7.4f   %7.4f   %7.4f\n"%
		   LDlogg[gindex], LDT[Tindex], LDfilterName[findex],
			  LDtable[gindex][Tindex][findex][1],
			  LDtable[gindex][Tindex][findex][2],
			  LDtable[gindex][Tindex][findex][3],
			  LDtable[gindex][Tindex][findex][4])
    out.close()

    return


def WriteIperpTable():
    """
    Write the IperpTable to a file named IperpTable.inspect  
    """
    outfile = "IperpTable.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file IperpTable.inspect.")

    delta = IperpT[1] - IperpT[0]
    out.write( "     %5.0f  %5.0f   %4.0f\n"% 
                             IperpT[0], IperpT[maxIperpTindex], delta)
    delta = Iperplogg[1] - Iperplogg[0]
    out.write( "        %4.1f  %4.1f  %4.1f\n"%
	                     Iperplogg[0], Iperplogg[maxIperpgindex], delta)
    nfilters = maxIperpfilterindex + 1
    outputline = "     %2ld  "%(nfilters)
    for findex in range(0, maxIperpfliterindex):
        outputline +=  " "
        outputline += IperpfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, maxIperpTindex):
        for gindex in range(0, maxIperpgindex):
            outputline = "%5.0f %4.2f"% (IperpT[Tindex], Iperplogg[gindex])
            for findex in range(0, maxIperpfilterindex):
                dummy = "  %10.3e"% (Iperptable[gindex][Tindex][findex])
                outputline +=  dummy
            outputline += ("\n")
            out.write("%s"% outputline)
    out.close()

    return

def WriteIBBfilterTable():
    """
    Write the IBBfilterTable to a file named IBBfilter.inspect
    """
    outfile =  "IBBfilterTable.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file IBBfilter.inspect.")

    delta = IBBT[1] - IBBT[0]
    out.write( "     %5.0f  %5.0f   %4.0f\n"% 
                             IBBT[0], IBBT[maxIBBTindex], delta)
    nfilters = maxIBBfilterindex + 1
    outputline = "     %2ld  "% (nfilters)
    for findex in range(0, maxIBBfilterindex):
        outputline +=  " "
        outputline += IBBfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, maxIBBTindex):
        outputline =  "%7.0f"% (IBBT[Tindex])
        for findex in range(0, maxIBBfilterindex):
            dummy =  "  %10.3e"% (IBBtable[Tindex][findex])
            outputline +=  dummy
    outputline +=  "\n"
    out.write("%s"% outputline)
    out.close()

    return


def WriteZzetaTable():
    """
    Write the ZBBzeta to a file named ZzetaTable.inspect  
    """

    outfile = "ZzetaTable.inspect" 
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file ZzetaTable.inspect.")

    out.write( "  %6ld  %6.4f\n"% maxBBzetaindex, deltaBBzeta)
    for i in range(0, maxBBzetaindex):
        BBzeta = i * deltaBBzeta
        out.write( " %7.4f   %14.7e\n"% BBzeta, ZBBzeta[i])

    out.close()

    return 

def InspectStar2Tiles():
    """
    This function writes out the various properties of the tiles
    covering star 2.  For convenience of inspection, it writes out
    the properties into several different files.
    """
    outfile = "Star2TilesA.inspect" 
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Star2TilesA.inspect.")

    out.write( "\n")
    out.write( "  star2.Ntiles   =   %5ld\n"%   star2.Ntiles)
    out.write( "  star2.volume   =  %10.3e\n"% star2.volume)
    out.write( "  star2.meanr    =  %10.3e\n"% star2.meanr)
    out.write( "  star2.meang    =  %10.3e\n"% star2.meang)
    out.write( "  star2.logg     =  %6.3f\n"%  star2.logg)
    out.write( "  star2.meanT    =  %5.0f\n"%  star2.meanT)
    out.write( "  star2.beta     =   %5.3f\n"%  star2.beta)
    out.write( "  star2.albedo   =  %5.2f\n"%  star2.albedo)

    out.write( "\n")
    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile        r         theta     phi          x            y            z\n")
    for i in range(1, star2.Ntiles):
        theta = T2theta[i] * ( 360.0 / 2*math.pi )
        phi   = T2phi[i]   * ( 360.0 / 2*math.pi )
        out.write( "%6ld  %12.5e %8.3f %8.3f   %12.5e %12.5e %12.5e\n"%
            i, T2r[i], theta, phi, T2x[i], T2y[i], T2z[i])

    out.close()

    outfile =  "Star2TilesB.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Star2TilesB.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile    T2gradV.r   T2gradV.t   T2gradV.p\n")
    for i in range(1, star2.Ntiles):
        out.write( "%6ld   %10.3e  %10.3e  %10.3e\n"%
	       i, T2gradV[i].r, T2gradV[i].theta, T2gradV[i].phi)

    out.close()

    outfile = "Star2TilesC.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Star2TilesC.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile normSphere.r normSphere.t normSphere.p normCart.x normCart.y normCart.z\n")
    for i in range(1, star2.Ntiles):
        out.write( "%6ld  %9.6f     %8.5f     %8.5f   %8.5f   %8.5f   %8.5f\n"%
            i, T2normSphere[i].r, T2normSphere[i].theta, T2normSphere[i].phi,
               T2normCart[i].x, T2normCart[i].y, T2normCart[i].z)

    out.close()

    outfile = "Star2TilesD.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Star2TilesD.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile        g      log(g)       dS          T\n")
    for i in range(1, star2.Ntiles):
        out.write( "%6ld   %10.3e %6.3f   %11.4e   %5.0f\n"%
            i, T2g[i], T2logg[i], T2dS[i], T2T[i] )

    out.close()

    outfile = "Star2TilesE.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Star2TilesE.inspect.")

    out.write( "  Star 2 Tile Specific Intensities:\n\n")

    outputline = "     Tile "
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            outputline +=  "   "
            outputline += orbit.filter[band]
            outputline +=  "   "
        else:
            outputline +=  "      "
            outputline += orbit.filter[band]
            outputline += "     "
    outputline +=  "\n"
    out.write("%s"% outputline)

    for i in range(1, star2.Ntiles):
        outputline =  "  %6ld "% (i)
        for band in range(1, orbit.nbands):
            dummy =  " %10.3e "% T2I[band][i]
            outputline += dummy
        outputline += "\n"
        out.write("%s"% outputline)
    out.close()

    return


def InspectDiskTiles( targetarea, 
                       nringMain, maintiles, 
                       nringEdge, edgetiles ):
    """
    This function writes out the various properties of the tiles
    covering star 2.  For convenience of inspection, it writes out
    the properties into several different files.
    """
    outfile = "DiskTilesA.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file DiskTilesA.inspect.")

    out.write( "\n")
    out.write( " Disk Tiles:\n")

    out.write( "\n")
    out.write( " disk.targetNtiles  =   %5ld\n"%  disk.targetNtiles)
    out.write( " disk.Ntiles        =   %5ld\n"%  disk.Ntiles)
    out.write( " disk.e             =   %5.3f\n"% disk.e)
    x = disk.zetazero * (360.0/2*math.pi)
    out.write( " disk.zetazero      =   %5.1f\n"% x)
    out.write( " disk.albedo        =   %5.2f\n"% disk.albedo)

    out.write( "\n")
    out.write( " maindisk.amin      =  %11.4e\n"% maindisk.amin)
    out.write( " maindisk.amax      =  %11.4e\n"% maindisk.amax)
    out.write( " maindisk.Hmax      =  %11.4e\n"% maindisk.Hmax)
    out.write( " maindisk.Hpow      =  %5.2f\n"%  maindisk.Hpow)
    out.write( " maindisk.Tamax     =  %10.3e\n"% maindisk.Tamax)
    out.write( " maindisk.Tamin     =  %10.3e\n"% maindisk.Tamin)
    out.write( " maindisk.Tpow      =  %5.2f\n"%  maindisk.Tpow)

    out.write( "\n")
    out.write( " diskedge.T         =   %9.3e\n"% diskedge.T)
    out.write( " diskedge.Tspot     =   %9.3e\n"% diskedge.Tspot)
    x = diskedge.ZetaMid * (360.0/2*math.pi)
    out.write( " diskedge.ZetaMid   =   %5.1f\n"%  x)
    x = diskedge.ZetaWidth * (360.0/2*math.pi)
    out.write( " diskedge.ZetaWidth =   %5.1f\n"%  x)

    if( control.diskrim == "ON" ):
        out.write( "\n")
        out.write( " diskrim.awidth     =  %11.4e\n"% diskrim.awidth)
        out.write( " diskrim.type       =   %s\n"%     diskrim.type)
        out.write( " diskrim.Hmax       =  %11.4e\n"% diskrim.Hmax)
        out.write( " diskrim.Hmin       =  %11.4e\n"% diskrim.Hmin)
        x = diskrim.ZetaHmax * (360.0/2*math.pi)
        out.write( " diskrim.ZetaHmax   =   %5.1f\n"%  x)
        out.write( " diskrim.Tmax       =   %9.3e\n"% diskrim.Tmax)
        out.write( " diskrim.Tmin       =   %9.3e\n"% diskrim.Tmin)
        x = diskrim.ZetaTmax * (360.0/2*math.pi)
        out.write( " diskrim.ZetaTmax   =   %5.1f\n"%  x)
        if( diskrim.type == "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, diskrim.points):
                x = diskrim.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%
                 i, x, diskrim.PointH[i], diskrim.PointT[i])
    if( control.disktorus ==  "ON"):
        out.write( "\n")
        out.write( " disktorus.azero    =  %11.4e\n"% disktorus.azero)
        out.write( " disktorus.awidth   =  %11.4e\n"% disktorus.awidth)
        out.write( " disktorus.type     =   %s\n"%     diskrim.type)
        out.write( " disktorus.Hmax     =  %11.4e\n"% disktorus.Hmax)
        out.write( " disktorus.Hmin     =  %11.4e\n"% disktorus.Hmin)
        x = disktorus.ZetaHmax * (360.0/2*math.pi)
        out.write( " disktorus.ZetaHmax =   %5.1f\n"%  x)
        out.write( " disktorus.Tmax     =   %9.3e\n"% disktorus.Tmax)
        out.write( " disktorus.Tmin     =   %9.3e\n"% disktorus.Tmin)
        x = disktorus.ZetaTmax * (360.0/2*math.pi)
        out.write( " disktorus.ZetaTmax =   %5.1f\n"%  x)
        if disktorus.type == "POINT":
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, disktorus.points):
                x = disktorus.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%
                 i, x, disktorus.PointH[i], disktorus.PointT[i])

    if( diskspot.nspots >= 1 ):
        out.write( "\n")
        out.write( "diskspot.npoints = %2ld\n"% diskspot.nspots)
        out.write( "   i  ZetaMin[i]  ZetaMax[i]   aMin[i]     aMax[i]  spotToverT[i] \n")
        for i in range(1, diskspot.nspots):
            x = diskspot.zetamin[i] * (360.0 / 2*math.pi)
            y = diskspot.zetamax[i] * (360.0 / 2*math.pi)
            out.write( "  %2ld    %5.1f       %5.1f    %10.4e  %10.4e    %5.2f\n"%
                   i, x, y, diskspot.amin[i], diskspot.amax[i], 
                   diskspot.spotToverT[i])

    out.write( "\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile      a      zeta     rho         h          x          y          z\n")
    for i in range(1, diak.Ntiles):
        zeta = TDiskZeta[i] * ( 360.0 / 2*math.pi )
        out.write( "%6ld %10.3e %5.1f %10.3e %10.3e %10.3e %10.3e %10.3e\n"%
	       i, TDiska[i] ,zeta, TDiskRho[i], TDiskH[i],
                       TDiskx[i], TDisky[i], TDiskz[i])
    out.close()

    outfile = "DiskTilesB.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file DiskTilesB.inspect.");

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile  normCyl.rho  normCyl.zeta   normCyl.h  normCart.x normCart.y normCart.z\n")
    for i in range(1, disk.Ntiles):
        out.write( "%6ld   %9.6f     %8.5f     %8.5f    %8.5f   %8.5f   %8.5f\n"%
            i, TDisknormCyl[i].rho, TDisknormCyl[i].zeta, TDisknormCyl[i].h,
               TDisknormCart[i].x,  TDisknormCart[i].y,   TDisknormCart[i].z)
    out.close()

    outfile = "DiskTilesC.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file DiskTilesC.inspect.")

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile      DiskdS         DiskT\n")
    for i in range(1, disk.Ntiles): 
        out.write( "%6ld  %12.4e  %12.4e\n"%
	      i, TDiskdS[i], TDiskT[i] )
    out.close()

    outfile = "DiskTilesD.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file DiskTilesD.inspect.")

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n\n"% edgetiles)

    out.write( "  Disk Tile Specific Intensities:\n\n")

    outputline =  "     Tile "
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            outputline += "   "
            outputline += orbit.filter[band]
            outputline  += "   "
        else:
            outputline += "      "
            outputline += orbit.filter[band]
            outputline += "     "
    outputline += "\n"
    out.write("%s"% outputline)

    for i in range(1, disk.Ntiles):
        outputline =  "  %6ld "% i
        for band in range(1, orbit.nbands):
            dummy =  " %10.3e "% (TDiskI[band][i])
            outputline +=  dummy
        outputline += "\n"
        out.write("%s"% outputline)
    out.close()

    return


def InspectYlimits():
    """
    This function writes out Grid.Topy[][] and Grid.Bottomy[][].
    """
    outfile = "ylimits.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file ylimits.inspect.")

    out.write( "\n")
    out.write( "  Grid.Nxtiles  =  %4ld\n"%    Grid.Nxtiles)
    out.write( "  Grid.Nztiles  =  %4ld\n"%    Grid.Nztiles)
    out.write( "  Grid.deltax   =  %12.4e\n"% Grid.deltax )
    out.write( "  Grid.deltaz   =  %12.4e\n"% Grid.deltaz )
    out.write( "  Grid.deltal   =  %12.4e\n"% Grid.deltal )
    out.write( "  Grid.xmin     =  %12.4e\n"% Grid.xmin   )
    out.write( "  Grid.xmax     =  %12.4e\n"% Grid.xmax   )
    out.write( "  Grid.ymin     =  %12.4e\n"% Grid.ymin   )
    out.write( "  Grid.ymax     =  %12.4e\n"% Grid.ymax   )
    out.write( "  Grid.zmin     =  %12.4e\n"% Grid.zmin   )
    out.write( "  Grid.zmax     =  %12.4e\n"% Grid.zmax   )

    out.write("\n")
    out.write("   ix   iz         x             z           Topy         Bottomy\n")
    for ix in range(1, Grid.Nxtiles):
        x = Grid.xmin + (ix - 1) * Grid.deltax
        for iz in range(1, Grid.Nztiles):
            z = Grid.zmin + (iz - 1) * Grid.deltaz
            out.write("  %3ld  %3ld   %12.5e  %12.5e  %12.5e  %12.5e\n"%
	            ix, iz, x, z, Grid.Topy[ix][iz], Grid.Bottomy[ix][iz])
    out.close()

    return


def InspectEscape( sunvector, f,
                    Star1Emitted, Star1Escape,
                    T2mu, T2Emitted, T2Escape,
                    TDiskmu, TDiskEmitted,
                    TDiskEscape ):
    """
    This function writes out diagnostic information about the
    flux emitted from the tiles and whether or not the flux
    is seen from the Earth.
    """
    if( control.star1 == "ON" ):
        outfile = "escape1.inspect"
        out = open(outfile, "w")
        if out == None:
            Quit("Cannot open file escape1.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)

        out.write("\n")
        out.write( "  filter       =  %s\n"% f)
        out.write( "  Star1Emitted = %12.5e\n"% Star1Emitted)
        out.write( "  Star1Escape  =  %6.4f\n\n"% Star1Escape)
        out.close()
  
    if( control.star2 == "ON" ):
        outfile = "escape2.inspect"
        out = open(outfile, "w")
        if out == None:
            Quit("Cannot open file escape2.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)
        out.write( "  filter       =  %s\n"% f)

        out.write("\n")
        out.write( "                       Emitted    Fraction\n")
        out.write( "  itile       mu         Flux     Escaping\n")
        for itile in range(1, start2.Ntiles):
            out.write("  %5ld   %8.5f  %12.4e   %6.4f\n"%
	        itile, T2mu[itile], T2Emitted[itile], T2Escape[itile] )
        out.close()

    if( control.disk == "ON" ):
        outfile =  "escapeDisk.inspect" 
        out = open(outfile, "w")
        if out == None:
            Quit("Cannot open file escapeDisk.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)
        out.write( "  filter       =  %s\n"% f)

        out.write("\n")
        out.write( "                       Emitted     Fraction\n")
        out.write( "  itile       mu         Flux      Escaping\n")
        for itile in range(1, disk.Ntiles):
            out.write("  %5ld   %8.5f  %12.4e    %6.4f\n"%
	            itile, TDiskmu[itile], TDiskEmitted[itile], 
                    TDiskEscape[itile] )
        out.close()

    return


def InspectHeatDiskBy1(TDiskTold, muA1toD,
                        DeltaT41toD, transmit1toD):
    """
    This function writes out diagnostic information about the
    heating of the outer accretion disk by star 1.
    """
    outfile = "HeatDiskBy1.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file HeatDiskBy1.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT41toDisk  1->Disk   Temperature\n")
    for itile in range(1, disk.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muA1toD[itile],
                DeltaT41toD[itile], transmit1toD[itile], TDiskT[itile] )
    out.close()

    return

def InspectHeatDiskByID(TDiskTold, muAidtoD,
                        DeltaT4idtoD, transmitidtoD):
    """
    This function writes out diagnostic information about the
    heating of the outer disk by the inner disk.
    """
    outfile = "HeatDiskByID.inspect"
    out = open(outfile, "w")
    if out == None:      
        Quit("Cannot open file HeatDiskByID.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT4idtoDisk id->Disk   Temperature\n")
    for itile in range(1, disk.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muAidtoD[itile],
                DeltaT4idtoD[itile], transmitidtoD[itile], TDiskT[itile] )
    out.close()

    return


def InspectHeatDiskByADC( whatside,
                           TDiskTold, muAADCtoD,
                           DeltaT4ADCtoD, transmitADCtoD):
    """
    This function writes out diagnostic information about the
    heating of the outer disk by the ADC.
    """
    if( whatside == "TOP"):
        outfile = "HeatDiskByTopADC.inspect"
    if( whatside == "BOTTOM"):
        outfile = "HeatDiskByBotADC.inspect"
    else:
        Quit("Unrecognized side in function InspectHeatDiskByADC.")
    out = open(outfile, "w")
    if( out== NULL):
        Quit("Cannot open file HeatDiskByxxxADC.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT4ADC ADC->Disk   Temperature\n")
    for itile in range(1, diak.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muAADCtoD[itile],
                DeltaT4ADCtoD[itile], transmitADCtoD[itile], TDiskT[itile] )
    out.close()

    return


def InspectHeat2By1( T2Told, muA1to2, 
                      DeltaT41to2, transmit1to2 ):
    """
    This function writes out diagnostic information about the
    heating of star 2 by star 1.
    """
    outfile = "Heat2By1.inspect"
    out = open(outfile, "w")
    if (out == NULL):
        Quit("Cannot open file Heat2By1.inspect.")

    out.write( "           Original                          Fraction    Heated\n")
    out.write( "   Tile  Temperature     muA    DeltaT41to2    1->2    Temperature\n")
    for itile in range(1, star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e   %5.4f   %8.0f\n"%
	         itile, T2Told[itile], muA1to2[itile], 
                 DeltaT41to2[itile], transmit1to2[itile], T2T[itile] )

    out.close()

    return


def InspectHeat2ByID( T2Told, muAidto2, 
                     DeltaT4idto2, transmitidto2):
    """
    This function writes out diagnostic information about the
    heating of star 2 by the inner disk.
    """
    outfile = "Heat2ByInDisk.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Heat2ByInDisk.inspect.")
    out.write( "           Original                          Fraction    Heated\n")
    out.write( "   Tile  Temperature     muA    DeltaT4idto2   id->2    Temperature\n")
    for itile in range(1, star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e   %5.4f   %8.0f\n"%
	         itile, T2Told[itile], muAidto2[itile], 
                 DeltaT4idto2[itile], transmitidto2[itile], T2T[itile] )

    out.close()

    return


def InspectHeat2ByDisk( T2Told, DeltaT4Dto2, T2T):
    """
    This function writes out diagnostic information about the
    heating of star 2 by the disk.
    """
    outfile =  "Heat2ByDisk.inspect"
    out = open(outfile, "w")
    if out == None:
        Quit("Cannot open file Heat2ByDisk.inspect.")

    out.write( "           Original                      Heated\n")
    out.write( "   Tile  Temperature    DeltaT4Dto2    Temperature\n")
    for itile in range(1, star2.Ntiles):
        out.write("  %5ld  %8.0f      %11.4e    %8.0f\n"%
	         itile, T2Told[itile], DeltaT4Dto2[itile], T2T[itile] )

    out.close()

    return


def InspectHeat2ByADC( whatside,
                        T2Told, muAADCto2,
                        DeltaT4ADCto2, transmitADCto2):
    """
    This function writes out diagnostic information about the
    heating of star 2 by the ADC.
    """
    if( whatside == "TOP"):
        outfile = "Heat2ByTopADC.inspect"
    if( whatside == "BOTTOM"):
        outfile = "Heat2ByBotADC.inspect"
    else:
        Quit("Unrecognized side in function InspectHeat2ByADC.")
    out = open(outfile, "w")   
    if (out == NULL):
        Quit("Cannot open file Heat2ByxxxADC.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA     DeltaT4ADC    ADC->2   Temperature\n")
    for itile in range(1, star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, T2Told[itile], muAADCto2[itile],
                DeltaT4ADCto2[itile], transmitADCto2[itile], T2T[itile] )
    out.close()

    return


def WriteData(band):
    """
    This function writes out the observed light curves to files
    called band1data.inspect, band2data.inspect, etc.
    """
    bandnum =  "%2ld"% (band)
    filename = "band"
    filename += bandnum
    filename += "data.inspect"
   
    out = open(filename, "w")
    if out == None:
        print("Cannot open file %s\n"% filename)
        Quit("")

    out.write( "Band %2ld from file %s\n"% band, data.filename[band])
    if( data.filter[band] == "SQUARE"):
        out.write( "Filter = %s  %5.0f  %5.0f\n"% data.filter[band], 
                            data.minlambda[band], data.maxlambda[band])
    else:
        out.write( "Filter = %s\n"% data.filter[band])

    for i in range(1, data.npoints[band]):
        out.write( " %5.4f  %11.3e  %11.2e\n"% data.phase[band][i],
	      data.flux[band][i], data.standdev[band][i])

    out.close()
   
    return
