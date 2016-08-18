# -*- coding: utf-8 -*-
"""
All the input functions are in this file.
"""
import sys
from .diskflux import maindisk
from .star1 import Star1
from .star2 import Star2
from .fitdata import ReadData 
from .diagnose import InspectInput
from .parmeter import CartVector, CylVector, SphereVector, filenames, flowcontrol, orbitparams, systemparams, star2spotparams, wholediskpars, diskedgepars
from .parmeter import diskrimpars, disktorusparams, diskspotpars, innerdiskpars, adcpars, thirdlightparams, dataparams, globalvar

global verbose 
MAXBANDPASSES  = 21
MAXZETAPOINTS  = 101
MAXPHASES     =  501
MAX2TILES     = 40506
MAXDISKTILES  = 40506
MAXFILTERS    =   12

def ReadInput():
    """
    This is the main function responsible for reading all input files, whether 
    data or parameters.
    """
    if( flowcontrol.diagnostics == "NOCHECKPARS"):
        ReadPars()
    else:
        ReadPars()
        CheckPars()
    ReadGDTable()
    ReadLDTable()
    ReadIperpTable()
    ReadIBBfilterTable()
    ReadZzetaTable()
    if( flowcontrol.diagnostics == "INSPECTINPUT"):
        InspectInput()
        sys.exit("Quit after INSPECTINPUT.")

    return

def ReadPars():
    """
    This function reads the input parameters from a file named 
    "parfile.dat".
    The file "input.txt" describes the parameter file.                    
    """
    out = open(filenames.parfile, "r")
    lines = out.readlines()
    if out == None:
        sys.exit("Cannot open the input parameter file.")

    while(lines != None):
        for line in lines:
            linesplit = line.split()
            nfields = len(linesplit)
            keyword = linesplit[0]
            if(nfields < 1):
                keyword = "" 
            elif( keyword == "END"):
                break
            elif( keyword == ""): 
                x = 1.0
            elif( keyword == "COMMENT="):
                x = 1.0

            elif( keyword == "VERBOSE="):
      	     globalvar.verbose =  linesplit[1]
            elif( keyword == "DIAGNOSTICS="):
                if( nfields < 4 ):
                    sys.exit("Too few parameters for keyword DIAGNOSTICS.")
                flowcontrol.diagnostics = "%s"%linesplit[1]
                flowcontrol.diagnosephase = "%lf"%linesplit[2]
                flowcontrol.diagnoseband =  "%s"%linesplit[3]

            elif( keyword == "STAR1="):
                flowcontrol.star1 = linesplit[1]
            elif( keyword == "STAR2="):
	           flowcontrol.star2 = linesplit[1]
            elif( keyword == "STAR2SPOTS="):
	           flowcontrol.star2spots = linesplit[1]
            elif( keyword == "DISK="):
	           flowcontrol.disk = linesplit[1]
            elif( keyword == "DISKRIM="):
	           flowcontrol.diskrim = linesplit[1]
            elif( keyword == "DISKTORUS="):
	           flowcontrol.disktorus = linesplit[1]
            elif( keyword == "INNERDISK="):
	           flowcontrol.innerdisk = linesplit[1]
            elif( keyword == "DISKSPOTS="):
	           flowcontrol.diskspots = linesplit[1]
            elif( keyword == "ADC="):
	           flowcontrol.adc = linesplit[1]
            elif( keyword == "THIRDLIGHT="):
	           flowcontrol.thirdlight = linesplit[1]
            elif( keyword == "IRRADIATION="):
	           flowcontrol.irradiation = linesplit[1]

            elif( keyword == "PHASES="):
                if( nfields < 4 ):
	               sys.exit("Too few parameters for keyword PHASES.")
                orbitparams.phasemin = "%lf"%linesplit[1]
                orbitparams.phasemax = "%lf"%linesplit[2]
                orbitparams.deltaphase = "%lf"%linesplit[3]
            elif( keyword == "PHASEOFFSET="):
                orbitparams.phaseoffset = "%lf"%linesplit[1]
            elif( keyword == "BANDPASS="):
                if( nfields < 3 ):
                    sys.exit("Too few parameters for keyword BANDPASS.")
                    orbitparams.nbands += 1
                if( orbitparams.nbands > (MAXBANDPASSES - 1) ):
	               sys.exit("Too many bandpasses.")
                if( linesplit[1] == "FILTER"):
                    orbitparams.filter = "%s"%linesplit[2]
                    orbitparams.minlambda = -1.0
                    orbitparams.maxlambda = -1.0
                elif( linesplit[1] == "SQUARE"):
                    if( nfields < 4 ):
                        sys.exit("Too few parameters for BANDPASS= SQUARE.")
                    orbitparams.filter = "SQUARE"
                    orbitparams.minlambda = "%lf"%linesplit[2]
                    orbitparams.maxlambda = "%lf"%linesplit[3]
                else:
	              sys.exit("BANDPASS: Unrecognized bandpass type.")
            elif( keyword == "NORMALIZE="):
                orbitparams.normalize = "%s"%linesplit[1]
                if( orbitparams.normalize == "MAXVALUE"):
                    if( nfields < 3 ):
                        sys.exit("Too few parameters for keyword NORMALIZE MAXVALUE.")
                orbitparams.normvalue = "%lf"%linesplit[2]
                if( orbitparams.normalize == "FITDATA"):
                    if( nfields < 3 ):
                        sys.exit()("Too few parameters for keyword NORMALIZE FITDATA.")
                    orbitparams.normfilter = "%s"%linesplit[2]
                if( orbitparams.normfilter == "SQUARE"):
                    if( nfields < 5 ):
                        sys.exit()("Too few parameters for keyword NORMALIZE FITDATA.");
                    orbitparams.normMinlambda = "%lf"%linesplit[3]
                    orbitparams.normMaxlambda = "%lf"%linesplit[4] 

            elif( keyword == "PERIOD="):
                systemparams.p = "%lf"%linesplit[1]
            elif( keyword == "K2="): 
                systemparams.K2 = "%lf"%linesplit[1]
            elif( keyword == "M1="):
                systemparams.M1 = "%lf"%linesplit[1]
            elif( keyword == "MASSRATIO="):
                systemparams.q = "%lf"%linesplit[1]
            elif( keyword == "INCLINATION="):
                systemparams.i = "%lf"%linesplit[1]


            elif( keyword == "STAR1LUM="):
	          Star1.L = "%lf"%linesplit[1]
            elif( keyword == "STAR1TEMP="):
	          Star1.T = "%lf"%linesplit[1]

            elif( keyword == "STAR2TILES="):
                Star2.targetNtiles = "%ld"%linesplit[1]
            elif( keyword == "STAR2TEMP="):
                Star2.meanT = "%lf"%linesplit[1]
            elif( keyword == "STAR2ALBEDO="):
                Star2.albedo = "%lf"%linesplit[1]

            elif( keyword == "STAR2SPOT="):
                if( nfields < 5 ):
	               sys.exit("Too few parameters in keyword STAR2SPOT.")
                star2spotparams.nspots += 1
                if( star2spotparams.nspots >= 20 ):
	               sys.exit("Too many star 2 spots.")
                star2spotparams.theta = "%lf"%linesplit[1]
                star2spotparams.phi   = "%lf"%linesplit[2]
                star2spotparams.radius = "%lf"%linesplit[3]
                star2spotparams.SpotToverStarT = "%lf"%linesplit[4]

            elif( keyword == "DISKTILES="):
                wholediskpars.targetNtiles = "%ld"%linesplit[1]
            elif( keyword == "DISKE="):
                wholediskpars.e = "%lf"%linesplit[1]
            elif( keyword == "DISKZETAZERO="):
                wholediskpars.zetazero = "%lf"%linesplit[1]
            elif( keyword == "DISKALBEDO="):
                wholediskpars.albedo = "%lf"%linesplit[1]

            elif( keyword == "MAINDISKA="):
                if( nfields < 3 ):
                    sys.exit("Too few parameters for keyword MAINDISKRHO.")
                maindisk.amin = "%lf"%linesplit[1] 
                maindisk.amax = "%lf"%linesplit[2]
            elif( keyword == "MAINDISKH="):
                if( nfields < 3 ):
                    sys.exit("Too few parameters for keyword MAINDISKH.")
                maindisk.Hmax = "%lf"%linesplit[1]
                maindisk.Hpow = "%lf"%linesplit[2]
            elif( keyword == "MAINDISKT="):
                if( nfields < 2 ):
                    sys.exit("Too few parameters for keyword MAINDISKT.");
                maindisk.Ttype = "%s"%linesplit[1]
                if( maindisk.Ttype == "POWERLAW"):
                    if( nfields < 4 ):
                        sys.exit("Too few parameters for keyword MAINDISKT.");
                    maindisk.Tpow = "%lf"%linesplit[1]
                    maindisk.maindiskL = "%lf"%linesplit[2]
                elif( maindisk.Ttype == "VISCOUS"):
                    if( nfields < 3 ):
                        sys.exit("Too few parameters for keyword MAINDISKT.")
                    maindisk.maindiskL = "%lf"%linesplit[1]
                else:
                    sys.exit("Unrecognized temperature distribution for MAINDISKT")

            elif( keyword == "DISKEDGET="):
                if( nfields < 5 ):
                    sys.exit("Too few parameters for keyword DISKRIMT.")
                diskedgepars.T = "%lf"%linesplit[1]
                diskedgepars.Tspot = "%lf"%linesplit[2]
                diskedgepars.ZetaMid = "%lf"%linesplit[3]
                diskedgepars.ZetaWidth = "%lf"%linesplit[4]

            elif( keyword == "INNERDISKT="):
                if( nfields < 2 ):
                    sys.exit("Too few parameters for keyword INNERDISKT.")
                innerdiskpars.T = "%lf"%linesplit[1]
            elif( keyword == "INNERDISKL="):
                if( nfields < 2 ):
                    sys.exit("Too few parameters for keyword INNERDISKL.")
                innerdiskpars.L = "%lf"%linesplit[1]

            elif( keyword == "DISKRIMAWIDTH="):
                diskrimpars.awidth = "%lf"%linesplit[1]
            elif( keyword == "DISKRIMPARS="):
                if( linesplit[1] == "SINUSOID"):
                    if( diskrimpars.type == "MISSING"):
                        diskrimpars.type = linesplit[1]
                    elif( diskrimpars.type != "SINUSOID"):
                        sys.exit("DISKRIMPARS: Inconsistent disk rim types.")
                    if( nfields < 8 ):
                        sys.exit("Too few parameters in keyword DISKRIMPARS.")
                    diskrimpars.Hmax = "%lf"%linesplit[2]
                    diskrimpars.Hmin = "%lf"%linesplit[3]
                    diskrimpars.ZetaHmax = "%lf"%linesplit[4]
                    diskrimpars.Tmax = "%lf"%linesplit[5]
                    diskrimpars.Tmin = "%lf"%linesplit[6]
                    diskrimpars.ZetaTmax = "%lf"%linesplit[7]
                elif( linesplit[1] == "POINT"):
                    if( diskrimpars.type == "MISSING"):
                        diskrimpars.type = linesplit[1]
                    elif( diskrimpars.type != "POINT"):
                        sys.exit("DISKRIMPARS: Inconsistent disk rim types.")
                    if ( nfields < 5 ):
                        sys.exit("Too few parameters in keyword DISKRIMPARS.")
                    diskrimpars.points += 1
                    if( diskrimpars.points > (MAXZETAPOINTS - 1) ):
                        sys.exit("DISKRIMPARS: Too many points in the POINT rim.")
                    diskrimpars.PointZeta = "%lf"%linesplit[2]
                    diskrimpars.PointH = "%lf"%linesplit[3]
                    diskrimpars.PointT = "%lf"%linesplit[4]
                else:
                    sys.exit("DISKRIMPARS: Unrecognized rim type.")

            elif( keyword == "DISKTORUSAZERO="):
                disktorusparams.azero = "%lf"%linesplit[1]
            elif( keyword == "DISKTORUSAWIDTH="):
                disktorusparams.awidth = "%lf"%linesplit[1]
            elif( keyword == "DISKTORUSPARS="):
                if( linesplit[1] == "SINUSOID"):
                    if( disktorusparams.type == "MISSING"):
                        disktorusparams.type = linesplit[1]
                    elif( disktorusparams.type != "SINUSOID"):
                        sys.exit("DISKTORUSPARS: Inconsistent disk torus types.")
                    if( nfields < 8 ):
                        sys.exit("Too few parameters in keyword DISKTORUSPARS.")
                    disktorusparams.Hmax = "%lf"%linesplit[2]
                    disktorusparams.Hmin = "%lf"%linesplit[3]
                    disktorusparams.ZetaHmax = "%lf"%linesplit[4]
                    disktorusparams.Tmax = "%lf"%linesplit[5]
                    disktorusparams.Tmin = "%lf"%linesplit[6]
                    disktorusparams.ZetaTmax = "%lf"%linesplit[7]
                elif( linesplit[1] == "POINT"):
                    if( disktorusparams.type == "MISSING"):
                        disktorusparams.type = linesplit[1]
                    elif( disktorusparams.type != "POINT"):
                        sys.exit("DISKTORUSPARS: Inconsistent disk torus types.")
                    if( nfields < 5 ):
                        sys.exit("Too few parameters in keyword DISKTORUSPARS.")
                    disktorusparams.points += 1
                    if( disktorusparams.points > (MAXZETAPOINTS - 1) ):
                        sys.exit("DISKTORUSPARS: Too many points in the POINT torus.")
                    disktorusparams.PointZeta = "%lf"%linesplit[2]
                    disktorusparams.PointH = "%lf"%linesplit[3]
                    disktorusparams.PointT = "%lf"%linesplit[4]
                else:
                    sys.exit("DISKTORUSPARS: Unrecognized torus type.")

            elif( keyword == "DISKSPOT="):
                if( nfields < 6 ):
                    sys.exit("Too few parameters in keyword DISKSPOT.");
                diskspotpars.nspots += 1
                if( diskspotpars.nspots >= 20 ):
                    sys.exit("Too many disk spots.")
                diskspotpars.zetamin = "%lf"%linesplit[1]
                diskspotpars.zetamax = "%lf"%linesplit[2]
                diskspotpars.amin = "%lf"%linesplit[3]
                diskspotpars.amax = "%lf"%linesplit[4]
                diskspotpars.spotToverT = "%lf"%linesplit[5]

            elif( keyword == "ADCLUM="):
                adcpars.L = "%lf"%linesplit[1]
            elif( keyword == "ADCHEIGHT="):
                adcpars.height = "%lf"%linesplit[1]

            elif( keyword == "3rdLIGHTPHASE="):
                thirdlightparams.orbphase = "%lf"%linesplit[1]
            elif( keyword == "3rdLIGHTFRACTION="):
                if( nfields < 4 ):
                    sys.exit("Too few parameters for keyword 3rdLIGHTFRACTION.")
                thirdlightparams.nbands += 1
                if( thirdlightparams.nbands > (MAXBANDPASSES - 1) ):
                    sys.exit("Too many 3rdLIGHT bandpasses.")
                if( linesplit[1] == "FILTER"):
                    thirdlightparams.filter = "%s"%linesplit[2]
                    thirdlightparams.minlambda = -1.0
                    thirdlightparams.maxlambda = -1.0
                    thirdlightparams.fraction = "%s"%linesplit[3]
                elif( linesplit[1] == "SQUARE"):
                    if( nfields < 5 ):
                        sys.exit("Too few parameters for BANDPASS= SQUARE.")
                    thirdlightparams.filter = "SQUARE"
                    thirdlightparams.minlambda = "%lf"%linesplit[2]
                    thirdlightparams.maxlambda = "%lf"%linesplit[3]
                    thirdlightparams.fraction  = "%lf"%linesplit[4] 
                else:
                    sys.exit("BANDPASS: Unrecognized bandpass type.")

            elif( keyword == "READDATA="):
                if( nfields < 4 ):
                    sys.exit("Too few parameters for keyword READDATA.")
                dataparams.nbands += 1
                if( dataparams.nbands > (MAXBANDPASSES - 1) ):
                    sys.exit("Too many data bandpasses.")
                if( linesplit[1] == "FILTER"):
                    dataparams.filter[dataparams.nbands] = "%s"%linesplit[2]
                    dataparams.minlambda = -1.0
                    dataparams.maxlambda = -1.0;
                    dataparams.filename = "%s"%linesplit[3]
                elif( linesplit[1] == "SQUARE"):
                    if( nfields < 5 ):
                        sys.exit("Too few parameters for READDATA= SQUARE.")
                    dataparams.filter[dataparams.nbands] == "SQUARE"
                    dataparams.minlambda[dataparams.nbands] = "%lf"%linesplit[2]
                    dataparams.maxlambda[dataparams.nbands] = "%lf"%linesplit[3]
                    dataparams.filename[dataparams.nbands] = "%lf"%linesplit[4] 
                else:
                    sys.exit("BANDPASS: Unrecognized bandpass type.")
                dataparams.npoints = 0
                ReadData( dataparams.nbands )

            else:
                print("Unrecognized keyword in get_data.\n")
                print("   keyword =%20s\n", keyword)
                sys.exit("")
            x = x
            out.close()
            return

def CheckPars():
    """
    This function checks the input parameters to insure that
    they are reasonable.
    """
    if( flowcontrol.diagnostics == "ON"):
        if( (flowcontrol.diagnosephase < -0.5) or (flowcontrol.diagnosephase > 1.0) ):
            sys.exit("DIAGNOSTICS: diagnosephase out of range.")
    if (( flowcontrol.star1 != "ON") and (flowcontrol.star1 != "OFF")):
        sys.exit("STAR1: neither ON nor OFF")
    if( (flowcontrol.star2 != "ON")) and (flowcontrol.star2 != "OFF"):
        sys.exit("STAR2: neither ON nor OFF")
    if (flowcontrol.star2spots != "ON") and (flowcontrol.star2 != "OFF"):
        sys.exit("STAR2SPOTS cannot be ON if STAR2 is OFF.")
    if( (flowcontrol.disk == "ON")) and (flowcontrol.disk == "OFF"):
        sys.exit("DISK: neither ON nor OFF")
    if( flowcontrol.diskrim != "ON") and (flowcontrol.diskrim != "OFF"):
        sys.exit("DISKRIM: neither ON nor OFF")
    if( flowcontrol.diskrim == "ON") and (flowcontrol.disk == "OFF"):
        sys.exit("DISKRIM cannot be ON if DISK is OFF.")
    if( flowcontrol.disktorus == "ON") and (flowcontrol.disk == "OFF"):
        sys.exit("DISKTORUS cannot be ON if DISK is OFF.")
    if( flowcontrol.disktorus !="ON") and (flowcontrol.disktorus != "OFF"):
        sys.exit("DISKTORUS: neither ON nor OFF")
    if( flowcontrol.innerdisk !="ON") and (flowcontrol.innerdisk != "OFF"):
        sys.exit("INNERDISK: neither ON nor OFF")
    if( flowcontrol.innerdisk == "ON") and (flowcontrol.disk == "OFF"):
        sys.exit("INNERDISK cannot be ON if DISK is OFF.")
    if( flowcontrol.diskspots != "ON") and (flowcontrol.diskspots != "OFF"):
        sys.exit("DISKSPOTS: neither ON nor OFF")
    if( flowcontrol.diskspots == "ON") and( flowcontrol.disk != "OFF"):
        sys.exit("DISKSPOTS cannot be ON if DISK is OFF.")
    if( flowcontrol.adc != "ON") and (flowcontrol.adc != "OFF"):
        sys.exit("ADC: neither ON nor OFF")
    if( flowcontrol.thirdlight != "ON") and ( flowcontrol.thirdlight != "OFF"):
        sys.exit("THIRDLIGHT: neither ON nor OFF")
    if( flowcontrol.irradiation != "ON") and (flowcontrol.irradiation != "OFF"):
        sys.exit("IRRADIATION: neither ON nor OFF")
    if( (orbitparams.phasemin < -0.5) or (orbitparams.phasemin > 1.0) ):
        sys.exit("PHASES: phasemin out of range.")
    if( (orbitparams.phasemax < -0.5) or (orbitparams.phasemax > 1.0) ):
        sys.exit("PHASES: phasemax out of range.")
    if( orbitparams.phasemin > orbitparams.phasemax ):
        sys.exit("PHASES: phasemax must be greater than or equal to phasemin.")
    if( (orbitparams.phaseoffset < -0.5) or (orbitparams.phaseoffset > 0.5) ):
        sys.exit("PHASEOFFSET: deltaphase must be ge -0.5 and le 0.5.")
    idummy = 1.0 + (orbitparams.phasemax - orbitparams.phasemin) / orbitparams.deltaphase
    if( idummy > MAXPHASES ):
        sys.exit("Number of orbital phases is greater than MAXPHASES.")
    if( orbitparams.nbands == 0):
        sys.exit("No bandpasses specified for the light curves.")
    for band in range(1, orbitparams.nbands):
        if( orbitparams.filter[band] == "SQUARE"):
            if( (orbitparams.minlambda[band] < 0.0) or (orbitparams.minlambda[band] > 30000.) ):
	           sys.exit("One of the BANDPASS= SQUARE minlambdas is out of range.")
            if( (orbitparams.maxlambda[band] < 0.0) or (orbitparams.maxlambda[band] > 30000.) ):
	           sys.exit("One of the BANDPASS= SQUARE maxlambdas is out of range.")
            if( orbitparams.maxlambda[band] <= orbitparams.minlambda[band] ):
	           sys.exit("BANDPASS: orbit.minlambda must be le than orbit.maxlambda.")
    if( orbitparams.normalize == "MISSING"):
        sys.exit("NORMALIZE keyword missing from parfile.dat.")
    if( orbitparams.normalize != "OFF"):
        if( orbitparams.normalize == "MAXVALUE"):
            if( orbitparams.normvalue <= 0.0 ):
	           sys.exit("NORMALIZE:  normalization value out of range.")
        if( orbitparams.normalize == "FITDATA"):
            found = 0
            for band in range(1, orbitparams.nbands):
                if( orbitparams.normfilter == orbitparams.filter[band]):
                    if( orbitparams.normfilter == "SQUARE"):
                        if( (orbitparams.normMinlambda == orbitparams.minlambda[band]) and (orbitparams.normMaxlambda == orbitparams.maxlambda[band]) ):
		                 found = 1
                    else:
		             found = 1
            if( found == 0 ):
	           sys.exit("No light curve calculated for FITDATA filter.")
            else:
                 sys.exit("Unrecognized normalization type for keyword NORMALIZE.")


    if( (systemparams.p < 0.001) or (systemparams.p > 365.0) ):
        sys.exit("Orbital period out of range.")
    if( (systemparams.K2 < 0.0) and (systemparams.M1 < 0.0) ):
        sys.exit("Either M1 or K2 but not both must be specified.")
    if( (systemparams.K2 > 0.0) and (systemparams.M1 > 0.0) ):
        sys.exit("Either M1 or K2 but not both must be specified.")
    if( systemparams.M1 < 0.0 ):
        if( (systemparams.K2 < 1.0) or (systemparams.K2 > 1000.0) ):
            sys.exit("K2 out of range.")

    if( systemparams.K2 < 0.0 ):
        if( (systemparams.M1 < 0.1) or (systemparams.M1 > 100.0) ):
            sys.exit("M1 out of range.")
    if( (systemparams.q < 0.01)  or  (systemparams.q > 1.0) ):
        sys.exit("Mass ratio must lie between 0.01 and 1.00.")
    if( (systemparams.i <= 0.0)  or  (systemparams.i > 90.0) ):
        sys.exit("Inclination must be gt 0.0 and le 90.0.")


    if( flowcontrol.star1 == "ON"):
        if( (Star1.L < 0.0) or (Star1.L > 1.0e40) ):
            sys.exit("Luminosity of star 1 out of range.")
        if( (Star1.T < 0.0) or (Star1.T > 3.0e7) ):
            sys.exit("Temperature of star 1 out of range.")


    if( (Star2.targetNtiles < 100) or (Star2.targetNtiles > 0.99 * MAX2TILES) ):
        sys.exit("STAR2TILES must be ge 100 and le 0.99*MAX2TILES.")
    if( (Star2.meanT < 3.0) or (Star2.meanT > 1.0e4) ):
        sys.exit("Temperature of star 2 out of range.")
    if( (Star2.albedo < 0.0) or (Star2.albedo > 1.0) ):
        sys.exit("STAR2ALBEDO out of range.")

    if( flowcontrol.star2spots == "ON"):
        if( star2spotparams.nspots <= 0): 
            sys.exit("STAR2SPOTS= ON, but no spots specified.")
        for i in range(1, star2spotparams.nspots):
            if( (star2spotparams.theta[i] < 0.0) or (star2spotparams.theta[i] > 180.0) ):
                sys.exit("STAR2SPOT: theta out of range.")
            if( (star2spotparams.phi[i] < 0.0) or (star2spotparams.phi[i] > 360.0) ):
                sys.exit("STAR2SPOT: phi out of range.")
            if( (star2spotparams.radius[i] <= 0.0) or (star2spotparams.radius[i] > 90.0) ):
                sys.exit("STAR2SPOT: radius out of range.")
            if( (star2spotparams.SpotToverStarT[i] <= 0.0) or (star2spotparams.SpotToverStarT[i] > 2.0) ):
                sys.exit("STAR2SPOT: SpotToverStarT out of range.")

    if( flowcontrol.disk == "ON"):
        if( (wholediskpars.targetNtiles < 100) or (wholediskpars.targetNtiles > 0.99 * MAXDISKTILES) ):
            sys.exit("DISKTILES must be ge 100 and le 0.99*MAXDISKTILES.")
        if( (wholediskpars.e < 0.0) or (wholediskpars.e >= 1.0) ):
            sys.exit("DISKE: disk.e out of range.")
        if( (wholediskpars.zetazero < 0.0) or (wholediskpars.zetazero >= 360.0) ):
            sys.exit("DISKZETAZERO: disk.zetazero out of range.")
        if( (wholediskpars.albedo < 0.0) or (wholediskpars.albedo > 1.0) ):
	      sys.exit("DISKALBEDO out of range.")

        if( (maindisk.amin <0.0) or (maindisk.amin > 0.6) ):
            sys.exit("MAINDISKA: maindisk.amin out of range.")
        if( (maindisk.amax <0.0) or (maindisk.amax > 0.6) ):
            sys.exit("MAINDISKA: maindisk.amax out of range.")
        if( maindisk.amax <= maindisk.amin ):
            sys.exit("MAINDISKA: maindisk.amax must be greater than maindisk.amin.");
        if( (maindisk.Hmax <= 0.0) or (maindisk.Hmax > maindisk.amax) ):
            sys.exit("MAINDISKH: maindisk.Hmax must be gt 0.0 and le maindisk.amax.")
        if( (maindisk.Hpow < 0.0) or (maindisk.Hpow > 2.0) ):
            sys.exit("MAINDISKH: maindisk.Hpow must be between 0.0 and 2.0.")
        if( (maindisk.maindiskL < 1.0e28) or (maindisk.maindiskL > 1.0e39)):
            sys.exit("MAINDISKT: maindisk luminosity must be between 1.0e28 and 1.0e39");
        if( maindisk.Ttype == "POWERLAW"):
            if( (maindisk.Tpow < -3.0) or (maindisk.Tpow > 3.0) ):
                sys.exit("MAINDISKT: maindisk.Tpow must be between -3.0 and 3.0.")
        if( (diskedgepars.T < 0.0) or (diskedgepars.T > 1.0e6) ):
            sys.exit("DISKEDGET: Edge T out of range.")
        if( (diskedgepars.Tspot < 0.0) or (diskedgepars.Tspot > 1.0e6) ):
            sys.exit("DISKEDGET: Tspot out of range.")
        if( (diskedgepars.ZetaMid < 0.0) or (diskedgepars.ZetaMid >= 360.0) ):
            sys.exit("DISKEDGET: ZetaMid must be ge 0.0 and lt 360.0.")
        if( (diskedgepars.ZetaWidth < 0.0) or (diskedgepars.ZetaWidth >= 360.0) ):
            sys.exit("DISKEDGET: ZetaWidth must be ge 0.0 and lt 360.0.")

    if( flowcontrol.innerdisk == "ON"):
        if( (innerdiskpars.T < 0.0) or (innerdiskpars.T > 1.0e7) ):
	      sys.exit("Inner disk temperature out of range.")
        if( (innerdiskpars.L < 500.0) or (innerdiskpars.L > 1.0e39) ):
	      sys.exit("Inner disk luminosity out of range.")

    if( flowcontrol.diskrim == "ON"):
        if( diskrimpars.awidth <= 0.0 ):
            sys.exit("DISKRIMAWIDTH is out of range.")
        if( diskrimpars.awidth > (maindisk.amax - maindisk.amin) ):
            sys.exit("DISKRIMAWIDTH is greater than the disk width.")
        if( diskrimpars.type == "SINUSOID"):
            if( (diskrimpars.Hmax < 0.0) or (diskrimpars.Hmax > maindisk.amax) ):
                sys.exit("DISKRIMPARS: Hmax must be ge 0.0 and le maindisk.amax.")
            if( (diskrimpars.Hmin < 0.0) or (diskrimpars.Hmin > maindisk.amax) ):
                sys.exit("DISKRIMPARS: Hmin must be ge 0.0 and le maindisk.amax.")
            if( diskrimpars.Hmax < diskrimpars.Hmin ):
                sys.exit("DISKRIMPARS: Hmax must be ge Hmin.")
            if( (diskrimpars.ZetaHmax < 0.0) or (diskrimpars.ZetaHmax >= 360.0) ):
                sys.exit("DISKRIMPARS: ZetaHmax must be ge 0.0 and lt 360.0.")
            if( (diskrimpars.Tmax < 0.0) or (diskrimpars.Tmax > 1.0e6) ):
                sys.exit("DISKRIMPARS: Tmax out of range.")
            if( (diskrimpars.Tmin < 0.0) or (diskrimpars.Tmin > 1.0e6) ):
                sys.exit("DISKRIMPARS: Tmin out of range.")
            if( diskrimpars.Tmax < diskrimpars.Tmin ):
                sys.exit("DISKRIMPARS: Tmax must be ge than Tmin.")
            if( (diskrimpars.ZetaTmax < 0.0) or (diskrimpars.ZetaTmax >= 360.0) ): 
	           sys.exit("DISKRIMPARS: ZetaTmax must be ge 0.0 and lt 360.0.")
        if( diskrimpars.type == "POINT"):
            for i in range(1, diskrimpars.points):
                if( (diskrimpars.PointZeta[i] < 0.0) or (diskrimpars.PointZeta[i] >= 360.0) ):
	               sys.exit("DISKRIMPARS: The PointZetas must be ge 0.0 and lt 360.0.")
                if( (diskrimpars.PointH[i] < 0.0) or (diskrimpars.PointH[i] > maindisk.amax) ):
 	               sys.exit("DISKRIMPARS: Rim H must be ge 0 and le maindisk.amax.")
                if( (diskrimpars.PointT[i] < 0.0) or (diskrimpars.PointT[i] > 1.0e6) ):
                    sys.exit("DISKRIMPARS: At least one rim T is out of range.")
        else:
            sys.exit("DISKRIMPARS: Unrecognized type.")

    if( flowcontrol.disktorus == "ON"):
        if( (disktorusparams.azero >= maindisk.amax) or (disktorusparams.azero <= maindisk.amin) ):
            sys.exit("DISKTORUSAZERO is outside the disk.")
        if( (disktorusparams.azero - 0.5 * disktorusparams.awidth) < maindisk.amin):
            sys.exit("DISKTORUSAWIDTH: torus extends past the disk inner edge.");
        if( (disktorusparams.azero + 0.5 * disktorusparams.awidth) > maindisk.amax):
            sys.exit("DISKTORUSAWIDTH: torus extends past the disk outer edge.")
        if( flowcontrol.diskrim == "ON"):
            if( (maindisk.amax - diskrimpars.awidth) < (disktorusparams.azero + 0.5 * disktorusparams.awidth) ):
                sys.exit("Disk rim and disk torus overlap.")
        if( disktorusparams.type == "SINUSOID"):
            if( (disktorusparams.Hmax < 0.0) or (disktorusparams.Hmax > maindisk.amax) ):
                sys.exit("DISKTORUSPARS: Hmax must be ge 0.0 and le maindisk.amax.")
            if( (disktorusparams.Hmin < 0.0) or (disktorusparams.Hmin > maindisk.amax) ):
                sys.exit("DISKTORUSPARS: Hmin must be ge 0.0 and le maindisk.amax.")
            if( disktorusparams.Hmax < disktorusparams.Hmin ):
                sys.exit("DISKTORUSPARS: Hmax must be ge Hmin.")
            if( (disktorusparams.ZetaHmax < 0.0) or (disktorusparams.ZetaHmax >= 360.0) ):
                sys.exit("DISKTORUSPARS: ZetaHmax must be ge 0.0 and lt 360.0.")
            if( (disktorusparams.Tmax < 0.0) or (disktorusparams.Tmax > 1.0e6) ):
                sys.exit("DISKTORUSPARS: Tmax out of range.")
            if( (disktorusparams.Tmin < 0.0) or (disktorusparams.Tmin > 1.0e6) ):
                sys.exit("DISKTORUSPARS: Tmin out of range.")
            if( disktorusparams.Tmax < disktorusparams.Tmin ):
                sys.exit("DISKTORUSPARS: Tmax must be ge than Tmin.")
            if( (disktorusparams.ZetaTmax < 0.0) or (disktorusparams.ZetaTmax >= 360.0) ): 
	           sys.exit("DISKTORUSPARS: ZetaTmax must be ge 0.0 and lt 360.0.")
        if( disktorusparams.type == "POINT"):
            for i in range(1, disktorusparams.points):
                if( (disktorusparams.PointZeta[i] < 0.0) or (disktorusparams.PointZeta[i] >= 360.0) ):
                    sys.exit("DISKTORUSPARS: The PointZetas must be ge 0.0 and lt 360.0.")
                if( (disktorusparams.PointH[i] < 0.0) or (disktorusparams.PointH[i] > maindisk.amax) ):
 	               sys.exit("DISKTORUSPARS: Torus H must be ge 0 and le maindisk.amax.")
                if( (disktorusparams.PointT[i] < 0.0) or (disktorusparams.PointT[i] > 1.0e6) ):
                     sys.exit("DISKTORUSPARS: At least one torus T is out of range.")
        else:
            sys.exit("DISKTORUSPARS: Unrecognized type.")

    if( flowcontrol.diskspots == "ON"):
        if( diskspotpars.nspots <= 0): 
            sys.exit("DISKSPOTS= ON, but no spots specified.")
        for i in range(1, diskspotpars.nspots):
            if( (diskspotpars.zetamin[i] < 0.0) or (diskspotpars.zetamin[i] > 360.0) ):
                sys.exit("diskspot.zetamin out of range.")
            if( (diskspotpars.zetamax[i] < 0.0) or (diskspotpars.zetamin[i] > 360.0) ):
                sys.exit("diskspot.zetamax out of range.")
            if( diskspotpars.amin[i] >= diskspotpars.amax[i] ):
                sys.exit("diskspot.amax must be greater than diskspot.amin.")
            if( (diskspotpars.amin[i] < 0.0) or (diskspotpars.amin[i] > 0.6) ):
                sys.exit("diskspot.amin out of range.")
            if( (diskspotpars.amax[i] < 0.0) or (diskspotpars.amax[i] > 0.6) ):
                sys.exit("diskspot.amin out of range.")
            if( (diskspotpars.spotToverT[i] < 0.0) or (diskspotpars.spotToverT[i] > 100.0) ):
                sys.exit("diskspot.spotToverT out of range.")

    if( flowcontrol.adc == "ON"):
        if( (adcpars.L < 0.0) or (adcpars.L > 1.0e40) ):
            sys.exit("Luminosity of the ADC is out of range.")
        if( (adcpars.height <= 0.0) or (adcpars.height > 0.5) ):
            sys.exit("Height of the point-approx ADC is out of range.")

    if( flowcontrol.thirdlight == "ON" ):
        if( (thirdlightparams.orbphase < -0.5) or (thirdlightparams.orbphase >= 1.0) ):
            sys.exit("3rdLIGHTPHASE out of range.")
        if( thirdlightparams.nbands <= 0 ):
            sys.exit("Third light is on but no 3rdLIGHTFRACTIONs specified.")
        for i in range(1, thirdlightparams.nbands):
            if( orbitparams.filter[band] == "SQUARE"):
                if( (thirdlightparams.minlambda[i] < 0.0) or (thirdlightparams.minlambda[i] > 30000.) ):
                    sys.exit("One of the 3rd light SQUARE minlambdas is out of range.")
                if( (thirdlightparams.maxlambda[i] < 0.0) or (thirdlightparams.maxlambda[i] > 30000.) ):
	               sys.exit("One of the 3rd light SQUARE maxlambdas is out of range.")
                if( thirdlightparams.maxlambda[i] <= orbitparams.minlambda[i] ):
	               sys.exit("The 3rd light orbit.minlambda must be le than orbit.maxlambda.")
            if( (thirdlightparams.fraction[i] < 0.0) or (thirdlightparams.fraction[i] >= 1.0) ):
                sys.exit("3rd light fraction must be ge 0.0 and lt 1.0.")

    if( globalvar.verbose != "ON") and (globalvar.verbose !="OFF"):
        sys.exit("VERBOSE must be ON or OFF.")
    return


def ReadGDTable():
    """
    This function read the table containing the gravity
    darkening coefficients.
    Note that this function assumes that the
    coefficients in the table are 4*beta, where 
       Teff = <Teff> * ( g / <g> )^beta
 
    The gravity darkening data are stored in the global variables
      long maxGDindex
      double GDT[], fourbeta[]
    """

    filename = "GDTable.dat"
    out = open(filename, "r")
    if out == None:
        sys.exit("Cannot open file GDTable.dat.")

    i = -1
    for inputline in out:
        if( inputline[0] != '*' ):
            i = i +1
            a = inputline.split()
            globalvar.GDT[i], globalvar.fourbeta[i] = ("%lf")%(a[0]), ("%lf")%(a[1])
            #sscanf(inputline, "%lf %lf", &GDT[i], &fourbeta[i] )
    globalvar.maxGDindex = i
    out.close()

    return


def ReadLDTable():
    """
    This function reads the limb darkening table LDTable.dat.   
    The limb darkening law is the four-parameter law advocated by
         A. Claret (2000, A&A, 363, 1081):
         I(mu) / I(1) =  1 - a1*(1 - mu**(1/2))
                           - a2*(1 - mu**(1)  )
                           - a3*(1 - mu**(3/2))
                           - a4*(1 - mu**(2)  )
 
    The file must have the format:
        
        Tmin    Tmax    deltaT    =  The minimum and maximum temperature
                                     and the temperature spacing.
        gmin    gmax    deltag    =  The minimum and maximum LOG g
                                     and the spacing in log g.
          N filtername1 filtername2 filtername3 ... filternameN
         T   log(g)   lambda    a1    a2    a3    a4
         .     .        .       .     .     .     .
         .     .        .       .     .     .     .
         .     .        .       .     .     .     .
 
 
    The limb darkening table is stored in the global variables
       long maxgindex, maxTindex, maxlindex
       double LDlogg[gindex], LDT[Tindex], LDlambda[lindex]
       double LDtable[gindex][Tindex][lindex][aindex]
   
    The maximum values of the indices in these arrays is set
    in header.h.  Be careful.
    """
    filename = "LDTable.dat"
    out = open(filename, "r")
    lines = out.readlines()    
    if(out == None):
        sys.exit("Cannot open file LDTable.dat.")

    while(True):
        #for i in out: #fgets(inputline, 80, out)
        if( lines[0][0] != '*' ):
            e = lines[0].split() #sscanf( inputline, "%lf %lf %lf", &Tmin, &Tmax, &deltaT);
            Tmin, Tmax, deltaT = float(e[0]), float(e[1]), float(e[2])
            maxLDTindex = ( Tmax - Tmin + 0.1 ) / deltaT
            for i in range(1, maxLDTindex):
                globalvar.LDT[i] = Tmin + i * deltaT
            b = lines[1].split()
            gmin, gmax, deltag = float(b[0]), float(b[1]), float(b[2])
            maxLDgindex = ( gmax - gmin + 0.001) / deltag
            for i in range(0, maxLDgindex):
                globalvar.LDlogg[i] = gmin + i * deltag
            c = lines[2].split()
            nfilters = int(c[0])
            LDfilterName = [0 for i in range(0,8)]
            LDfilterName[0] = float(c[1])
            LDfilterName[1] = float(c[2])
            LDfilterName[2] = float(c[3]) 
            LDfilterName[3] = float(c[4])
            LDfilterName[4] = float(c[5])
            LDfilterName[5] = float(c[6])
            LDfilterName[6] = float(c[8])
            LDfilterName[7] = float(c[9])
            maxLDfilterindex = nfilters - 1
            if( maxLDfilterindex > (MAXFILTERS - 1) ):
	          sys.exit("Too many filters in the LD table.")
            break

    d = lines[3].split()
    logg = float(d[0])
    T =    float(d[1])
    filtername = (d[2])
    a = [0 for i in range(0,5)]
    a[1] = float(d[3]) 
    a[2] = float(d[4])
    a[3] = float(d[5])
    a[4] = float(d[6])
    gindex = (logg - gmin + 0.1) / deltag
    if( (gindex < 0) or (gindex > maxLDgindex) ):
         sys.exit("gindex out of range in ReadLDTable.")
    Tindex = ( T - Tmin + 1.0) / deltaT
    if( (Tindex < 0) or (Tindex > maxLDTindex) ):
         sys.exit("Tindex out of range in ReadLDTable.")
    findex = -1
    for i in range(0, maxLDfilterindex):
        if( filtername == globalvar.LDfilterName[i]): 
            findex = i
    if( (findex < 0) or (findex > maxLDfilterindex) ):
        sys.exit("Unrecognized filter name in ReadLDTable.")
    for i in range(1,4):
        globalvar.LDtable[gindex][Tindex][findex][i] = a[i]
    out.close()
    return


def ReadIperpTable():
    """
    This function reads IperpTable.dat.   The file containing the Iperp
    data must have the format:
        Tmin    Tmax    deltaT    =  The minimum and maximum temperature
                                     and the temperature spacing.
        gmin    gmax    deltag    =  The minimum and maximum LOG g
                                     and the spacing in log g.
          n    filter1  filter2 filter3 filter4 ........
         t   logg    Iz(lambda1)  Iz(lambda2)  Iz(lambda3)  Iz(lambda4)......
         .     .           .            .            .            .
         .     .           .            .            .            .
         .     .           .            .            .            .
 
     n is the number of filters.  "filter1" "filter2", etc. are the
       names of the filters, eg U B V R I
 
    The Iperp data are stored in the global variables
       long maxIperpgindex, maxIperpTindex, maxIperplindex
       double Iperplogg[gindex], IperpT[Tindex], Iperplambda[lindex]
       double Iperptable[gindex][Tindex][lindex]
    """
    filename = "IperpTable.dat"
    out = open(filename, "r")
    lines = out.readlines()
    if out == None:
        sys.exit("Cannot open file IperpTable.dat.")

    b = lines[0].split()#fscanf( in, "%lf  %lf  %lf", &Tmin, &Tmax, &deltaT);
    Tmin ,Tmax, deltaT = float(b[0]), float(b[1]), float(b[2])
    maxIperpTindex = ( Tmax - Tmin + 1.0 ) / deltaT
    for Tindex in range(0, maxIperpTindex):
        globalvar.IperpT[Tindex] = Tmin + Tindex * deltaT
    c = lines[1].split()#fscanf( in, "%lf  %lf  %lf", &gmin, &gmax, &deltag);
    gmin, gmax, deltag = float(c[0]), float(c[1]), float(c[2])
    maxIperpgindex = ( gmax - gmin + 0.01) / deltag
    for gindex in range(0, maxIperpgindex):
        globalvar.Iperplogg[gindex] = gmin + gindex * deltag
    d = lines[2].split()
    nfilters = int(d[0])
    IperpfilterName = [0 for i in range(0,8)]
    IperpfilterName[0] = float(d[1])
    IperpfilterName[1] = float(d[2])
    IperpfilterName[2] = float(d[3]) 
    IperpfilterName[3] = float(d[4])
    IperpfilterName[4] = float(d[5])
    IperpfilterName[5] = float(d[6])
    IperpfilterName[6] = float(d[8])
    IperpfilterName[7] = float(d[9])
    maxIperpfilterindex = nfilters - 1
    if( globalvar.maxIperpfilterindex > (MAXFILTERS - 1) ):
        sys.exit("Too many filters in the Iperp table.")

    while(True):
        e = lines[3].split()
        xT = float(d[0])
        xlogg =    float(e[1])
        xIperp = [0 for i in range(0, 8)]
        xIperp[0] = float(e[2])
        xIperp[1] = float(e[3]) 
        xIperp[2] = float(e[4])
        xIperp[3] = float(e[5])
        xIperp[4] = float(e[6])
        xIperp[5] = float(e[7])
        xIperp[6] = float(e[8])
        xIperp[7] = float(e[9])
        gindex = (xlogg - gmin + 0.01) / deltag
        if( (gindex < 0) or (gindex > maxIperpgindex) ):
            sys.exit("gindex out of range in ReadIperpTable.")
        Tindex = ( xT - Tmin + 1.0) / deltaT
        if( (Tindex < 0) or (Tindex > maxIperpTindex) ):
            sys.exit("Tindex out of range in ReadIperpTable.")
        for findex in range(0, maxIperpfilterindex):
            globalvar.Iperptable[gindex][Tindex][findex] = xIperp[findex]
    out.close()

    return

def ReadIBBfilterTable():
    """
    This function reads the intensities of a black body
    observed through a filter.
 
    The file name must be "IBBfilter.dat" and must have the format:
        
        Tmin    Tmax    deltaT    =  The minimum and maximum temperature
                                     and the temperature spacing.
          N filtername1 filtername2 filtername3 ... filternameN
         T   Ifilter1  Ifilter2  Ifilter3 ... IfilterN
         .      .         .         .            .
         .      .         .         .            .
         .      .         .         .            .
 
    The filter can have comment lines beginning with a "*" at the
    beginning.
    """
    filename = "IBBfilterTable.dat"
    out = open(filename, "r")
    lines = out.readlines()
    if out == None:
        sys.exit("Cannot open file IBBfilter.dat.");

    while(True):
        if( lines[0][0] != '*' ):
            e = lines[0].split() #sscanf( inputline, "%lf %lf %lf", &Tmin, &Tmax, &deltaT);
            IBBTmin, IBBTmax, IBBdeltaT = float(e[0]), float(e[1]), float(e[2])
            globalvar.maxLDTindex = ( IBBTmax - IBBTmin + 0.1 ) / IBBdeltaT
            for i in range(1, globalvar.maxIBBindex):
                globalvar.IBBT[i] = IBBTmin + i * IBBdeltaT
            c = lines[1].split()
            nfilters = int(c[0])
            IBBfilterName = [0 for i in range(0, 8)]
            IBBfilterName[0] = float(c[1])
            IBBfilterName[1] = float(c[2])
            IBBfilterName[2] = float(c[3]) 
            IBBfilterName[3] = float(c[4])
            IBBfilterName[4] = float(c[5])
            IBBfilterName[5] = float(c[6])
            IBBfilterName[6] = float(c[8])
            IBBfilterName[7] = float(c[9])
            globalvar.maxLDfilterindex = nfilters - 1
            if( globalvar.maxIBBfilterindex > (MAXFILTERS - 1) ):
	          sys.exit("Too many filters in the LD table.")
            break

    d = lines[2].split()
    xT = float(d[0])
    xIBB = [0 for i in range(0, 8)]
    xIBB[0] = float(d[1])
    xIBB[1]=  float(d[2])
    xIBB[2] = float(d[3]) 
    xIBB[3] = float(d[4])
    xIBB[4] = float(d[5])
    xIBB[5] = float(d[6])
    xIBB[6] = float(d[7])
    xIBB[7] = float(d[8]) 
    Tindex = ( xT - IBBTmin + 1.0) / IBBdeltaT
    if( (Tindex < 0) or (Tindex > globalvar.maxIBBTindex) ):
        sys.exit("Tindex out of range in ReadIBBTable.")
    for findex in range(0,globalvar.maxIBBfilterindex):
         globalvar.IBBtable[Tindex][findex] = xIBB[findex]
    out.close()

    return

def ReadZzetaTable():
    """ 
    Read the ZBBzeta table

    """
    filename = "ZzetaTable.dat"
    out = open(filename, "r")
    lines = out.readlines()
    if out == None:
        sys.exit("Cannot open file ZzetaTable.dat.")
   
    b = lines[0].split()   
    maxBBzetaindex , deltaBBzeta= float(b[0]), float(b[1])   
    for i in range(1, maxBBzetaindex):
        dummy, globalvar.ZBBzeta[i] = float(lines[i].split()[0]), float(lines[i].split()[1])
    globalvar.BBzetamax = maxBBzetaindex * deltaBBzeta

    out.close()

    return