# -*- coding: utf-8 -*-
"""
This file contains functions used to diagnose and test the
line broadening program.
"""
import math
import sys
from .diskflux import maindisk
from .star1 import Star1
from .star2 import Star2
from .parmeter import filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, XYGrid, dataparams, globalvar


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
    if( dataparams.nbands > 0):
        for band in range(1, dataparams.nbands):
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
        sys.exit("Cannot open file parfile.inspect.")
    out.write( "DIAGNOSTICS=       %s  %5.3f  %s\n"% flowcontrol.diagnostics,
                                                      flowcontrol.diagnosephase,
                                                      flowcontrol.diagnoseband)

    out.write( "\n")
    out.write( "STAR1=             %s\n"% flowcontrol.star1)
    out.write( "STAR2=             %s\n"% flowcontrol.star2)
    out.write( "DISK=              %s\n"% flowcontrol.disk)
    out.write( "DISKRIM=           %s\n"% flowcontrol.diskrim)
    out.write( "DISKTORUS=         %s\n"% flowcontrol.disktorus)
    out.write( "INNERDISK=         %s\n"% flowcontrol.innerdisk)
    out.write( "DISKSPOTS=         %s\n"% flowcontrol.diskspots)
    out.write( "ADC=               %s\n"% flowcontrol.adc)
    out.write( "THIRDLIGHT=        %s\n"% flowcontrol.thirdlight)
    out.write( "IRRADIATION=       %s\n"% flowcontrol.irradiation)

    out.write( "\n")
    out.write( "PHASES=            %6.3f  %6.3f  %6.3f\n"% 
	                   orbitparams.phasemin, orbitparams.phasemax, orbitparams.deltaphase)
    out.write( "PHASEOFFSET=      %7.4f\n"%  orbitparams.phaseoffset)
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"): 
            out.write( "BANDPASS=          SQUARE  %6.0f  %6.0f\n"%
                   orbitparams.minlambda[band], orbitparams.maxlambda[band])
        else:
            out.write( "BANDPASS=          FILTER   %s\n"% orbitparams.filter[band])
    if( orbitparams.normalize == "OFF"):
        out.write("NORMALIZE=         %s\n"% orbitparams.normalize)
    else:
        if( orbitparams.normalize == "MAXVALUE"):
            out.write( "NORMALIZE=   %s  %12.4e\n"% orbitparams.normalize,
	                                         orbitparams.normvalue)
        if( orbitparams.normalize == "FITDATA"):
            if( orbitparams.normfilter == "SQUARE"):
                out.write( "NORMALIZE=   %s  %s  %7.1f  %7.1f\n"%
		             orbitparams.normalize, orbitparams.normfilter, 
                             orbitparams.normMinlambda, orbitparams.normMaxlambda)
            else:
	          out.write( "NORMALIZE=   %s  %s\n"%
		   orbitparams.normalize, orbitparams.normfilter)

    out.write( "\n")
    out.write("PERIOD=            %10.8f\n"% systemparams.p)
    out.write( "K2=                %6.2f\n"%  systemparams.K2)
    out.write( "MASSRATIO=         %5.3f\n"%  systemparams.q)
    out.write( "INCLINATION=       %5.2f\n"%  systemparams.i)

    if( flowcontrol.star1 == "ON"):
        out.write( "\n")
        out.write( "STAR1LUM=          %9.3e\n"% Star1.L)
        out.write( "STAR1TEMP=         %9.3e\n"% Star1.T)

    if( flowcontrol.star2 == "ON"):
        out.write( "\n")
        out.write( "STAR2TILES=       %5ld\n"%  Star2.targetNtiles)
        out.write( "STAR2TEMP=        %5.0f\n"% Star2.meanT)
        out.write( "STAR2ALBEDO=      %5.2f\n"% Star2.albedo)

    if( flowcontrol.disk == "ON"):
        out.write( "\n")
        out.write( "DISKTILES=         %5ld\n"%   wholediskpars.targetNtiles)
        out.write( "DISKE=             %5.3f\n"%  wholediskpars.e)
        out.write( "DISKZETAZERO=     %5.1f\n"%  wholediskpars.zetazero)
        out.write( "DISKALBEDO=        %5.2f\n"% wholediskpars.albedo)

        out.write("\n")
        out.write( "MAINDISKA=         %5.3f  %5.3f\n"% maindisk.amin,
	                                              maindisk.amax)
        out.write( "MAINDISKH=         %5.3f  %4.1f\n"% maindisk.Hmax,
	                                               maindisk.Hpow)
        out.write( "MAINDISKT=        %5.0f  %4.1f\n"% maindisk.Tamax,
	                                              maindisk.Tpow)
        out.write( "\n")
        out.write( "DISKEDGET=        %5.0f  %5.0f  %5.1f  %5.1f\n"%
                                      diskedgepars.T, diskedgepars.Tspot, 
                                      diskedgepars.ZetaMid, diskedgepars.ZetaWidth)

    if( flowcontrol.diskrim == "ON"):
        out.write( "\n")
        out.write( "DISKRIMAWIDTH=     %5.3f\n"%      diskrimpars.awidth)
        if( diskrimpars.type == "SINUSOID"):
            out.write( "DISKRIMPARS=       %s   %5.3f %5.3f %5.1f   %6.0f %6.0f %5.1f\n"%
	                  diskrimpars.type, 
                          diskrimpars.Hmax, diskrimpars.Hmin, diskrimpars.ZetaHmax,
                          diskrimpars.Tmax, diskrimpars.Tmin, diskrimpars.ZetaTmax)
        if( diskrimpars.type == "POINT"):
            for i in range(1, diskrimpars.points):
                out.write( "DISKRIMPARS=       %s  %5.1f  %5.3f  %7.1f\n"%
		     diskrimpars.type, diskrimpars.PointZeta[i],
                         diskrimpars.PointH[i], diskrimpars.PointT[i])

    if( flowcontrol.disktorus == "ON"):
        out.write( "\n")
        out.write( "DISKTORUSAZERO=    %5.3f\n"%    disktorusparams.azero)
        out.write( "DISKTORUSAWIDTH=   %5.3f\n"%    disktorusparams.awidth)
        if( disktorusparams.type == "SINUSOID"):
            out.write( "DISKTORUSPARS=     %s   %5.3f %5.3f %5.1f   %6.0f %6.0f %5.1f\n"%
	                  disktorusparams.type, 
                          disktorusparams.Hmax, disktorusparams.Hmin, disktorusparams.ZetaHmax,
                          disktorusparams.Tmax, disktorusparams.Tmin, disktorusparams.ZetaTmax)
        if( disktorusparams.type == "POINT"):
            for i in range(1, disktorusparams.points):
                out.write( "DISKTORUSPARS=     %s  %5.1f  %5.3f  %7.1f\n"%
		     disktorusparams.type, disktorusparams.PointZeta[i],
                         disktorusparams.PointH[i], disktorusparams.PointT[i])

    if( flowcontrol.diskspots == "ON"):
        out.write( "\n")
        for i in range(1, diskspotpars.nspots):
            out.write( "DISKSPOT=    %5.1f  %5.1f  %5.3f  %5.3f %6.3f\n"%
	             diskspotpars.zetamin[i], diskspotpars.zetamax[i],
                     diskspotpars.amin[i], diskspotpars.amax[i], diskspotpars.spotToverT[i])

    if( flowcontrol.innerdisk == "ON"):
        out.write( "\n")
        out.write( "INNERDISKT=       %9.2e\n"% innerdiskpars.T)
        out.write( "INNERDISKL=       %9.2e\n"% innerdiskpars.L)

    if( flowcontrol.adc == "ON"):
        out.write( "\n")
        out.write( "ADCL=               =  %10.3e\n"% adcpars.L)
        out.write( "ADCHEIGHT=          =  %11.4e\n"% adcpars.height)

    if( flowcontrol.thirdlight == "ON"):
        out.write( "\n")
        out.write( "3rdLIGHTPHASE =   %6.3f\n"%  thirdlightparams.orbphase)
        for i in range(1, thirdlightparams.nbands):
            if( thirdlightparams.filter[i] == "SQUARE"):
                out.write( "3rdLIGHTFRACTION=  %s %6.0f %6.0f  %5.3f\n"%
		        thirdlightparams.filter[i], thirdlightparams.minlambda[i],
                        thirdlightparams.maxlambda[i], thirdlightparams.fraction[i] )
            else:
	          out.write( "3rdLIGHTFRACTION=  FILTER    %s           %5.3f\n"%
		        thirdlightparams.filter[i], thirdlightparams.fraction[i] )

    if( dataparams.nbands > 0):
        out.write( "\n")
        for i in range(1, dataparams.nbands):
            if( dataparams.filter[i] == "SQUARE"):
                out.write( "READDATA=  %s %6.0f %6.0f   %s\n"%
		        dataparams.filter[i], dataparams.minlambda[i],
                        dataparams.maxlambda[i], dataparams.filename[i] )
            else:
	          out.write( "READDATA=  FILTER     %s           %s\n"%
		        dataparams.filter[i], dataparams.filename[i] )

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
        sys.exit("Cannot open file GDTable.inspect.")

    for i in range(0, globalvar.maxGDindex):
        out.write("  %5.0f  %5.3f\n"% globalvar.GDT[i], globalvar.fourbeta[i])

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
        sys.exit("Cannot open file LDTable.inspect.")

    delta = globalvar.LDT[1] - globalvar.LDT[0]
    out.write( "        %5.0f  %5.0f   %4.0f\n"% 
                             globalvar.LDT[0], globalvar.LDT[globalvar.maxLDTindex], delta)
    delta = globalvar.LDlogg[1] - globalvar.LDlogg[0]
    out.write( "        %4.1f  %4.1f  %4.1f\n"%
	                     globalvar.LDlogg[0], globalvar.LDlogg[globalvar.maxLDgindex], delta)
    nfilters = globalvar.maxLDfilterindex + 1
    outputline = "     %2ld  "% (nfilters)
    for findex in range(0, globalvar.maxLDfilterindex):
        outputline +=  " "
        outputline += LDfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, globalvar.maxLDTindex):
        for gindex in range(0, globalvar.maxLDgindex):
            for findex in range(0, globalvar.maxLDfilterindex):
                out.write( " %5.1f  %5.0f  %s   %7.4f   %7.4f   %7.4f   %7.4f\n"%
		   globalvar.LDlogg[gindex], globalvar.LDT[Tindex], LDfilterName[findex],
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
        sys.exit("Cannot open file IperpTable.inspect.")

    delta = globalvar.IperpT[1] - globalvar.IperpT[0]
    out.write( "     %5.0f  %5.0f   %4.0f\n"% 
                             globalvar.IperpT[0], globalvar.IperpT[globalvar.maxIperpTindex], delta)
    delta = globalvar.Iperplogg[1] - globalvar.Iperplogg[0]
    out.write( "        %4.1f  %4.1f  %4.1f\n"%
	                     globalvar.Iperplogg[0], globalvar.Iperplogg[globalvar.maxIperpgindex], delta)
    nfilters = globalvar.maxIperpfilterindex + 1
    outputline = "     %2ld  "%(nfilters)
    for findex in range(0, globalvar.maxIperpfliterindex):
        outputline +=  " "
        outputline += IperpfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, globalvar.maxIperpTindex):
        for gindex in range(0, globalvar.maxIperpgindex):
            outputline = "%5.0f %4.2f"% (globalvar.IperpT[Tindex], globalvar.Iperplogg[gindex])
            for findex in range(0, globalvar.maxIperpfilterindex):
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
        sys.exit("Cannot open file IBBfilter.inspect.")

    delta = globalvar.IBBT[1] - globalvar.IBBT[0]
    out.write( "     %5.0f  %5.0f   %4.0f\n"% 
                             globalvar.IBBT[0], globalvar.IBBT[globalvar.maxIBBTindex], delta)
    nfilters = globalvar.maxIBBfilterindex + 1
    outputline = "     %2ld  "% (nfilters)
    for findex in range(0, globalvar.maxIBBfilterindex):
        outputline +=  " "
        outputline += IBBfilterName[findex]
    outputline += "\n"
    out.write("%s"% outputline)

    for Tindex in range(0, globalvar.maxIBBTindex):
        outputline =  "%7.0f"% (globalvar.IBBT[Tindex])
        for findex in range(0, globalvar.maxIBBfilterindex):
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
        sys.exit("Cannot open file ZzetaTable.inspect.")

    out.write( "  %6ld  %6.4f\n"% globalvar.maxBBzetaindex, globalvar.deltaBBzeta)
    for i in range(0, globalvar.maxBBzetaindex):
        globalvar.BBzeta.append(i * globalvar.deltaBBzeta)
        out.write( " %7.4f   %14.7e\n"% globalvar.BBzeta, globalvar.ZBBzeta[i])

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
        sys.exit("Cannot open file Star2TilesA.inspect.")

    out.write( "\n")
    out.write( "  star2.Ntiles   =   %5ld\n"%   Star2.Ntiles)
    out.write( "  star2.volume   =  %10.3e\n"% Star2.volume)
    out.write( "  star2.meanr    =  %10.3e\n"% Star2.meanr)
    out.write( "  star2.meang    =  %10.3e\n"% Star2.meang)
    out.write( "  star2.logg     =  %6.3f\n"%  Star2.logg)
    out.write( "  star2.meanT    =  %5.0f\n"%  Star2.meanT)
    out.write( "  star2.beta     =   %5.3f\n"%  Star2.beta)
    out.write( "  star2.albedo   =  %5.2f\n"%  Star2.albedo)

    out.write( "\n")
    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile        r         theta     phi          x            y            z\n")
    for i in range(1, Star2.Ntiles):
        theta = globalvar.T2theta[i] * ( 360.0 / 2*math.pi )
        phi   = globalvar.T2phi[i]   * ( 360.0 / 2*math.pi )
        out.write( "%6ld  %12.5e %8.3f %8.3f   %12.5e %12.5e %12.5e\n"%
            i, globalvar.T2r[i], theta, phi, globalvar.T2x[i], globalvar.T2y[i], globalvar.T2z[i])

    out.close()

    outfile =  "Star2TilesB.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file Star2TilesB.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile    T2gradV.r   T2gradV.t   T2gradV.p\n")
    for i in range(1, Star2.Ntiles):
        out.write( "%6ld   %10.3e  %10.3e  %10.3e\n"%
	       i, T2gradV[i].r, T2gradV[i].theta, T2gradV[i].phi)

    out.close()

    outfile = "Star2TilesC.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file Star2TilesC.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile normSphere.r normSphere.t normSphere.p normCart.x normCart.y normCart.z\n")
    for i in range(1, Star2.Ntiles):
        out.write( "%6ld  %9.6f     %8.5f     %8.5f   %8.5f   %8.5f   %8.5f\n"%
            i, T2normSphere[i].r, T2normSphere[i].theta, T2normSphere[i].phi,
               T2normCart[i].x, T2normCart[i].y, T2normCart[i].z)

    out.close()

    outfile = "Star2TilesD.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file Star2TilesD.inspect.")

    out.write( "  Star 2 Tiles:\n\n")
    out.write( "  tile        g      log(g)       dS          T\n")
    for i in range(1, Star2.Ntiles):
        out.write( "%6ld   %10.3e %6.3f   %11.4e   %5.0f\n"%
            i, globalvar.T2g[i], globalvar.T2logg[i], globalvar.T2dS[i], globalvar.T2T[i] )

    out.close()

    outfile = "Star2TilesE.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file Star2TilesE.inspect.")

    out.write( "  Star 2 Tile Specific Intensities:\n\n")

    outputline = "     Tile "
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            outputline +=  "   "
            outputline += orbitparams.filter[band]
            outputline +=  "   "
        else:
            outputline +=  "      "
            outputline += orbitparams.filter[band]
            outputline += "     "
    outputline +=  "\n"
    out.write("%s"% outputline)

    for i in range(1, Star2.Ntiles):
        outputline =  "  %6ld "% (i)
        for band in range(1, orbitparams.nbands):
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
        sys.exit("Cannot open file DiskTilesA.inspect.")

    out.write( "\n")
    out.write( " Disk Tiles:\n")

    out.write( "\n")
    out.write( " disk.targetNtiles  =   %5ld\n"%  wholediskpars.targetNtiles)
    out.write( " disk.Ntiles        =   %5ld\n"%  wholediskpars.Ntiles)
    out.write( " disk.e             =   %5.3f\n"% wholediskpars.e)
    x = wholediskpars.zetazero * (360.0/2*math.pi)
    out.write( " disk.zetazero      =   %5.1f\n"% x)
    out.write( " disk.albedo        =   %5.2f\n"% wholediskpars.albedo)

    out.write( "\n")
    out.write( " maindisk.amin      =  %11.4e\n"% maindisk.amin)
    out.write( " maindisk.amax      =  %11.4e\n"% maindisk.amax)
    out.write( " maindisk.Hmax      =  %11.4e\n"% maindisk.Hmax)
    out.write( " maindisk.Hpow      =  %5.2f\n"%  maindisk.Hpow)
    out.write( " maindisk.Tamax     =  %10.3e\n"% maindisk.Tamax)
    out.write( " maindisk.Tamin     =  %10.3e\n"% maindisk.Tamin)
    out.write( " maindisk.Tpow      =  %5.2f\n"%  maindisk.Tpow)

    out.write( "\n")
    out.write( " diskedge.T         =   %9.3e\n"% diskedgepars.T)
    out.write( " diskedge.Tspot     =   %9.3e\n"% diskedgepars.Tspot)
    x = diskedgepars.ZetaMid * (360.0/2*math.pi)
    out.write( " diskedge.ZetaMid   =   %5.1f\n"%  x)
    x = diskedgepars.ZetaWidth * (360.0/2*math.pi)
    out.write( " diskedge.ZetaWidth =   %5.1f\n"%  x)

    if( flowcontrol.diskrim == "ON" ):
        out.write( "\n")
        out.write( " diskrim.awidth     =  %11.4e\n"% diskrimpars.awidth)
        out.write( " diskrim.type       =   %s\n"%     diskrimpars.type)
        out.write( " diskrim.Hmax       =  %11.4e\n"% diskrimpars.Hmax)
        out.write( " diskrim.Hmin       =  %11.4e\n"% diskrimpars.Hmin)
        x = diskrimpars.ZetaHmax * (360.0/2*math.pi)
        out.write( " diskrim.ZetaHmax   =   %5.1f\n"%  x)
        out.write( " diskrim.Tmax       =   %9.3e\n"% diskrimpars.Tmax)
        out.write( " diskrim.Tmin       =   %9.3e\n"% diskrimpars.Tmin)
        x = diskrimpars.ZetaTmax * (360.0/2*math.pi)
        out.write( " diskrim.ZetaTmax   =   %5.1f\n"%  x)
        if( diskrimpars.type == "POINT"):
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, diskrimpars.points):
                x = diskrimpars.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%
                 i, x, diskrimpars.PointH[i], diskrimpars.PointT[i])
    if( flowcontrol.disktorus ==  "ON"):
        out.write( "\n")
        out.write( " disktorus.azero    =  %11.4e\n"% disktorusparams.azero)
        out.write( " disktorus.awidth   =  %11.4e\n"% disktorusparams.awidth)
        out.write( " disktorus.type     =   %s\n"%     diskrimpars.type)
        out.write( " disktorus.Hmax     =  %11.4e\n"% disktorusparams.Hmax)
        out.write( " disktorus.Hmin     =  %11.4e\n"% disktorusparams.Hmin)
        x = disktorusparams.ZetaHmax * (360.0/2*math.pi)
        out.write( " disktorus.ZetaHmax =   %5.1f\n"%  x)
        out.write( " disktorus.Tmax     =   %9.3e\n"% disktorusparams.Tmax)
        out.write( " disktorus.Tmin     =   %9.3e\n"% disktorusparams.Tmin)
        x = disktorusparams.ZetaTmax * (360.0/2*math.pi)
        out.write( " disktorus.ZetaTmax =   %5.1f\n"%  x)
        if disktorusparams.type == "POINT":
            out.write( "   i     Zeta[i]       H[i]         T[i]\n")
            for i in range(1, disktorusparams.points):
                x = disktorusparams.PointZeta[i] * (360.0/2*math.pi)
                out.write( "  %2ld     %5.1f     %10.4e    %8.1f\n"%
                 i, x, disktorusparams.PointH[i], disktorusparams.PointT[i])

    if( diskspotpars.nspots >= 1 ):
        out.write( "\n")
        out.write( "diskspot.npoints = %2ld\n"% diskspotpars.nspots)
        out.write( "   i  ZetaMin[i]  ZetaMax[i]   aMin[i]     aMax[i]  spotToverT[i] \n")
        for i in range(1, diskspotpars.nspots):
            x = diskspotpars.zetamin[i] * (360.0 / 2*math.pi)
            y = diskspotpars.zetamax[i] * (360.0 / 2*math.pi)
            out.write( "  %2ld    %5.1f       %5.1f    %10.4e  %10.4e    %5.2f\n"%
                   i, x, y, diskspotpars.amin[i], diskspotpars.amax[i], 
                   diskspotpars.spotToverT[i])

    out.write( "\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile      a      zeta     rho         h          x          y          z\n")
    for i in range(1, wholediskpars.Ntiles):
        zeta = TDiskZeta[i] * ( 360.0 / 2*math.pi )
        out.write( "%6ld %10.3e %5.1f %10.3e %10.3e %10.3e %10.3e %10.3e\n"%
	       i, globalvar.TDiska[i] ,zeta, globalvar.TDiskRho[i], globalvar.TDiskH[i],
                       globalvar.TDiskx[i], globalvar.TDisky[i], globalvar.TDiskz[i])
    out.close()

    outfile = "DiskTilesB.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file DiskTilesB.inspect.");

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile  normCyl.rho  normCyl.zeta   normCyl.h  normCart.x normCart.y normCart.z\n")
    for i in range(1, wholediskpars.Ntiles):
        out.write( "%6ld   %9.6f     %8.5f     %8.5f    %8.5f   %8.5f   %8.5f\n"%
            i, TDisknormCyl[i].rho, TDisknormCyl[i].zeta, TDisknormCyl[i].h,
               TDisknormCart[i].x,  TDisknormCart[i].y,   TDisknormCart[i].z)
    out.close()

    outfile = "DiskTilesC.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file DiskTilesC.inspect.")

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n"% edgetiles)

    out.write( "\n")
    out.write( "  Tile      DiskdS         DiskT\n")
    for i in range(1, wholediskpars.Ntiles): 
        out.write( "%6ld  %12.4e  %12.4e\n"%
	      i, globalvar.TDiskdS[i], globalvar.TDiskT[i] )
    out.close()

    outfile = "DiskTilesD.inspect"
    out = open(outfile, "w")
    if out == None:
        sys.exit("Cannot open file DiskTilesD.inspect.")

    out.write( " Disk Tiles:\n\n")
    out.write( " targetarea = %12.4e\n"% targetarea)
    out.write( " nringMain  = %5ld\n"% nringMain)
    out.write( " maintiles  = %5ld\n"% maintiles)
    out.write( " nringEdge  = %5ld\n"% nringEdge)
    out.write( " edgetiles  = %5ld\n\n"% edgetiles)

    out.write( "  Disk Tile Specific Intensities:\n\n")

    outputline =  "     Tile "
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            outputline += "   "
            outputline += orbitparams.filter[band]
            outputline  += "   "
        else:
            outputline += "      "
            outputline += orbitparams.filter[band]
            outputline += "     "
    outputline += "\n"
    out.write("%s"% outputline)

    for i in range(1, wholediskpars.Ntiles):
        outputline =  "  %6ld "% i
        for band in range(1, orbitparams.nbands):
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
        sys.exit("Cannot open file ylimits.inspect.")

    out.write( "\n")
    out.write( "  Grid.Nxtiles  =  %4ld\n"%    XYGrid.Nxtiles)
    out.write( "  Grid.Nztiles  =  %4ld\n"%    XYGrid.Nztiles)
    out.write( "  Grid.deltax   =  %12.4e\n"% XYGrid.deltax )
    out.write( "  Grid.deltaz   =  %12.4e\n"% XYGrid.deltaz )
    out.write( "  Grid.deltal   =  %12.4e\n"% XYGrid.deltal )
    out.write( "  Grid.xmin     =  %12.4e\n"% XYGrid.xmin   )
    out.write( "  Grid.xmax     =  %12.4e\n"% XYGrid.xmax   )
    out.write( "  Grid.ymin     =  %12.4e\n"% XYGrid.ymin   )
    out.write( "  Grid.ymax     =  %12.4e\n"% XYGrid.ymax   )
    out.write( "  Grid.zmin     =  %12.4e\n"% XYGrid.zmin   )
    out.write( "  Grid.zmax     =  %12.4e\n"% XYGrid.zmax   )

    out.write("\n")
    out.write("   ix   iz         x             z           Topy         Bottomy\n")
    for ix in range(1, XYGrid.Nxtiles):
        x = XYGrid.xmin + (ix - 1) * XYGrid.deltax
        for iz in range(1, XYGrid.Nztiles):
            z = XYGrid.zmin + (iz - 1) * XYGrid.deltaz
            out.write("  %3ld  %3ld   %12.5e  %12.5e  %12.5e  %12.5e\n"%
	            ix, iz, x, z, XYGrid.Topy[ix][iz], XYGrid.Bottomy[ix][iz])
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
    if( flowcontrol.star1 == "ON" ):
        outfile = "escape1.inspect"
        out = open(outfile, "w")
        if out == None:
            sys.exit("Cannot open file escape1.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)

        out.write("\n")
        out.write( "  filter       =  %s\n"% f)
        out.write( "  Star1Emitted = %12.5e\n"% Star1Emitted)
        out.write( "  Star1Escape  =  %6.4f\n\n"% Star1Escape)
        out.close()
  
    if( flowcontrol.star2 == "ON" ):
        outfile = "escape2.inspect"
        out = open(outfile, "w")
        if out == None:
            sys.exit("Cannot open file escape2.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)
        out.write( "  filter       =  %s\n"% f)

        out.write("\n")
        out.write( "                       Emitted    Fraction\n")
        out.write( "  itile       mu         Flux     Escaping\n")
        for itile in range(1, Star2.Ntiles):
            out.write("  %5ld   %8.5f  %12.4e   %6.4f\n"%
	        itile, T2mu[itile], T2Emitted[itile], T2Escape[itile] )
        out.close()

    if( flowcontrol.disk == "ON" ):
        outfile =  "escapeDisk.inspect" 
        out = open(outfile, "w")
        if out == None:
            sys.exit("Cannot open file escapeDisk.inspect.")

        out.write("\n")
        out.write( "  sunvector.x  = %7.4f\n"% sunvector.x)
        out.write( "  sunvector.y  = %7.4f\n"% sunvector.y)
        out.write( "  sunvector.z  = %7.4f\n"% sunvector.z)
        out.write( "  filter       =  %s\n"% f)

        out.write("\n")
        out.write( "                       Emitted     Fraction\n")
        out.write( "  itile       mu         Flux      Escaping\n")
        for itile in range(1, wholediskpars.Ntiles):
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
        sys.exit("Cannot open file HeatDiskBy1.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT41toDisk  1->Disk   Temperature\n")
    for itile in range(1, wholediskpars.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muA1toD[itile],
                DeltaT41toD[itile], transmit1toD[itile], globalvar.TDiskT[itile] )
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
        sys.exit("Cannot open file HeatDiskByID.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT4idtoDisk id->Disk   Temperature\n")
    for itile in range(1, wholediskpars.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muAidtoD[itile],
                DeltaT4idtoD[itile], transmitidtoD[itile], globalvar.TDiskT[itile] )
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
        sys.exit("Unrecognized side in function InspectHeatDiskByADC.")
    out = open(outfile, "w")
    if( out== None):
        sys.exit("Cannot open file HeatDiskByxxxADC.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA   DeltaT4ADC ADC->Disk   Temperature\n")
    for itile in range(1, wholediskpars.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, TDiskTold[itile], muAADCtoD[itile],
                DeltaT4ADCtoD[itile], transmitADCtoD[itile], globalvar.TDiskT[itile] )
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
    if (out == None):
        sys.exit("Cannot open file Heat2By1.inspect.")

    out.write( "           Original                          Fraction    Heated\n")
    out.write( "   Tile  Temperature     muA    DeltaT41to2    1->2    Temperature\n")
    for itile in range(1, Star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e   %5.4f   %8.0f\n"%
	         itile, T2Told[itile], muA1to2[itile], 
                 DeltaT41to2[itile], transmit1to2[itile], globalvar.T2T[itile] )

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
        sys.exit("Cannot open file Heat2ByInDisk.inspect.")
    out.write( "           Original                          Fraction    Heated\n")
    out.write( "   Tile  Temperature     muA    DeltaT4idto2   id->2    Temperature\n")
    for itile in range(1, Star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e   %5.4f   %8.0f\n"%
	         itile, T2Told[itile], muAidto2[itile], 
                 DeltaT4idto2[itile], transmitidto2[itile], globalvar.T2T[itile] )

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
        sys.exit("Cannot open file Heat2ByDisk.inspect.")

    out.write( "           Original                      Heated\n")
    out.write( "   Tile  Temperature    DeltaT4Dto2    Temperature\n")
    for itile in range(1, Star2.Ntiles):
        out.write("  %5ld  %8.0f      %11.4e    %8.0f\n"%
	         itile, T2Told[itile], DeltaT4Dto2[itile], globalvar.T2T[itile] )

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
        sys.exit("Unrecognized side in function InspectHeat2ByADC.")
    out = open(outfile, "w")   
    if (out == None):
        sys.exit("Cannot open file Heat2ByxxxADC.inspect.")

    out.write( "           Original                           Fraction     Heated\n")
    out.write( "   Tile  Temperature     muA     DeltaT4ADC    ADC->2   Temperature\n")
    for itile in range(1, Star2.Ntiles):
        out.write("  %5ld  %8.0f     %7.4f   %11.4e    %5.4f    %8.0f\n"%
                itile, T2Told[itile], muAADCto2[itile],
                DeltaT4ADCto2[itile], transmitADCto2[itile], globalvar.T2T[itile] )
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
        sys.exit("")

    out.write( "Band %2ld from file %s\n"% band, dataparams.filename[band])
    if( dataparams.filter[band] == "SQUARE"):
        out.write( "Filter = %s  %5.0f  %5.0f\n"% dataparams.filter[band], 
                            dataparams.minlambda[band], dataparams.maxlambda[band])
    else:
        out.write( "Filter = %s\n"% dataparams.filter[band])

    for i in range(1, dataparams.npoints[band]):
        out.write( " %5.4f  %11.3e  %11.2e\n"% dataparams.phase[band][i],
	      dataparams.flux[band][i], dataparams.standdev[band][i])

    out.close()
   
    return
