# -*- coding: utf-8 -*-
"""
All the input functions are in this file.
"""

def ReadInput():
    """
    This is the main function responsible for reading all input files, whether 
    data or parameters.
    """
    if( control.diagnostics == "NOCHECKPARS"):
        ReadPars()
    else:
        ReadPars()
        CheckPars()
    ReadGDTable()
    ReadLDTable()
    ReadIperpTable()
    ReadIBBfilterTable()
    ReadZzetaTable()
    if( control.diagnostics == "INSPECTINPUT"):
        InspectInput()
        Quit("Quit after INSPECTINPUT.")

    return

def CheckPars():
    """
    This function checks the input parameters to insure that
    they are reasonable.
    """
    if( control.diagnostics == "ON"):
        if( (control.diagnosephase < -0.5) or (control.diagnosephase > 1.0) ):
            Quit("DIAGNOSTICS: diagnosephase out of range.")
    if (( control.star1 != "ON") and (control.star1 != "OFF")):
        Quit("STAR1: neither ON nor OFF")
    if( (control.star2 != "ON")) and (control.star2 != "OFF"):
        Quit("STAR2: neither ON nor OFF")
    if (control.star2spots != "ON") and (control.star2 != "OFF"):
        Quit("STAR2SPOTS cannot be ON if STAR2 is OFF.")
    if( (control.disk == "ON")) and (control.disk == "OFF"):
        Quit("DISK: neither ON nor OFF")
    if( control.diskrim != "ON") and (control.diskrim != "OFF"):
        Quit("DISKRIM: neither ON nor OFF")
    if( control.diskrim == "ON") and (control.disk == "OFF"):
        Quit("DISKRIM cannot be ON if DISK is OFF.")
    if( control.disktorus == "ON") and (control.disk == "OFF"):
        Quit("DISKTORUS cannot be ON if DISK is OFF.")
    if( control.disktorus !="ON") and (control.disktorus != "OFF"):
        Quit("DISKTORUS: neither ON nor OFF")
    if( control.innerdisk !="ON") and (control.innerdisk != "OFF"):
        Quit("INNERDISK: neither ON nor OFF")
    if( control.innerdisk == "ON") and (control.disk == "OFF"):
        Quit("INNERDISK cannot be ON if DISK is OFF.")
    if( control.diskspots != "ON") and (control.diskspots != "OFF"):
        Quit("DISKSPOTS: neither ON nor OFF")
    if( control.diskspots == "ON") and( control.disk != "OFF"):
        Quit("DISKSPOTS cannot be ON if DISK is OFF.")
    if( control.adc != "ON") and (control.adc != "OFF"):
        Quit("ADC: neither ON nor OFF")
    if( control.thirdlight != "ON") and ( control.thirdlight != "OFF"):
        Quit("THIRDLIGHT: neither ON nor OFF")
    if( control.irradiation != "ON") and (control.irradiation != "OFF"):
        Quit("IRRADIATION: neither ON nor OFF")
    if( (orbit.phasemin < -0.5) or (orbit.phasemin > 1.0) ):
        Quit("PHASES: phasemin out of range.")
    if( (orbit.phasemax < -0.5) or (orbit.phasemax > 1.0) ):
        Quit("PHASES: phasemax out of range.")
    if( orbit.phasemin > orbit.phasemax ):
        Quit("PHASES: phasemax must be greater than or equal to phasemin.")
    if( (orbit.phaseoffset < -0.5) or (orbit.phaseoffset > 0.5) ):
        Quit("PHASEOFFSET: deltaphase must be ge -0.5 and le 0.5.")
    idummy = 1.0 + (orbit.phasemax - orbit.phasemin) / orbit.deltaphase
    if( idummy > MAXPHASES ):
        Quit("Number of orbital phases is greater than MAXPHASES.")
    if( orbit.nbands == 0):
        Quit("No bandpasses specified for the light curves.")
    for band in range(1, orbit.nbands):
        if( orbit.filter[band] == "SQUARE"):
            if( (orbit.minlambda[band] < 0.0) or (orbit.minlambda[band] > 30000.) ):
	           Quit("One of the BANDPASS= SQUARE minlambdas is out of range.")
            if( (orbit.maxlambda[band] < 0.0) or (orbit.maxlambda[band] > 30000.) ):
	           Quit("One of the BANDPASS= SQUARE maxlambdas is out of range.")
            if( orbit.maxlambda[band] <= orbit.minlambda[band] ):
	           Quit("BANDPASS: orbit.minlambda must be le than orbit.maxlambda.")
    if( orbit.normalize == "MISSING"):
        Quit("NORMALIZE keyword missing from parfile.dat.")
    if( orbit.normalize != "OFF"):
        if( orbit.normalize == "MAXVALUE"):
            if( orbit.normvalue <= 0.0 ):
	           Quit("NORMALIZE:  normalization value out of range.")
        if( orbit.normalize == "FITDATA"):
            found = 0
            for band in range(1, orbit.nbands):
                if( orbit.normfilter == orbit.filter[band]):
                    if( orbit.normfilter == "SQUARE"):
                        if( (orbit.normMinlambda == orbit.minlambda[band]) and (orbit.normMaxlambda == orbit.maxlambda[band]) ):
		                 found = 1
                    else:
		             found = 1
            if( found == 0 ):
	           Quit("No light curve calculated for FITDATA filter.")
            else:
                 Quit("Unrecognized normalization type for keyword NORMALIZE.")


    if( (syspars.p < 0.001) or (syspars.p > 365.0) ):
        Quit("Orbital period out of range.")
    if( (syspars.K2 < 0.0) and (syspars.M1 < 0.0) ):
        Quit("Either M1 or K2 but not both must be specified.")
    if( (syspars.K2 > 0.0) and (syspars.M1 > 0.0) ):
        Quit("Either M1 or K2 but not both must be specified.")
    if( syspars.M1 < 0.0 ):
        if( (syspars.K2 < 1.0) or (syspars.K2 > 1000.0) ):
            Quit("K2 out of range.")

    if( syspars.K2 < 0.0 ):
        if( (syspars.M1 < 0.1) or (syspars.M1 > 100.0) ):
            Quit("M1 out of range.")
    if( (syspars.q < 0.01)  or  (syspars.q > 1.0) ):
        Quit("Mass ratio must lie between 0.01 and 1.00.")
    if( (syspars.i <= 0.0)  or  (syspars.i > 90.0) ):
        Quit("Inclination must be gt 0.0 and le 90.0.")


    if( control.star1 == "ON"):
        if( (star1.L < 0.0) or (star1.L > 1.0e40) ):
            Quit("Luminosity of star 1 out of range.")
        if( (star1.T < 0.0) or (star1.T > 3.0e7) ):
            Quit("Temperature of star 1 out of range.")


    if( (star2.targetNtiles < 100) or (star2.targetNtiles > 0.99 * MAX2TILES) ):
        Quit("STAR2TILES must be ge 100 and le 0.99*MAX2TILES.")
    if( (star2.meanT < 3.0) or (star2.meanT > 1.0e4) ):
        Quit("Temperature of star 2 out of range.")
    if( (star2.albedo < 0.0) or (star2.albedo > 1.0) ):
        Quit("STAR2ALBEDO out of range.")

    if( control.star2spots == "ON"):
        if( star2spot.nspots <= 0): 
            Quit("STAR2SPOTS= ON, but no spots specified.")
        for i in range(1, star2spot.nspots):
            if( (star2spot.theta[i] < 0.0) or (star2spot.theta[i] > 180.0) ):
                Quit("STAR2SPOT: theta out of range.")
            if( (star2spot.phi[i] < 0.0) or (star2spot.phi[i] > 360.0) ):
                Quit("STAR2SPOT: phi out of range.")
            if( (star2spot.radius[i] <= 0.0) or (star2spot.radius[i] > 90.0) ):
                Quit("STAR2SPOT: radius out of range.")
            if( (star2spot.SpotToverStarT[i] <= 0.0) or (star2spot.SpotToverStarT[i] > 2.0) ):
                Quit("STAR2SPOT: SpotToverStarT out of range.")

    if( control.disk == "ON"):
        if( (disk.targetNtiles < 100) or (disk.targetNtiles > 0.99 * MAXDISKTILES) ):
            Quit("DISKTILES must be ge 100 and le 0.99*MAXDISKTILES.")
        if( (disk.e < 0.0) or (disk.e >= 1.0) ):
            Quit("DISKE: disk.e out of range.")
        if( (disk.zetazero < 0.0) or (disk.zetazero >= 360.0) ):
            Quit("DISKZETAZERO: disk.zetazero out of range.")
        if( (disk.albedo < 0.0) or (disk.albedo > 1.0) ):
	      Quit("DISKALBEDO out of range.")

        if( (maindisk.amin <0.0) or (maindisk.amin > 0.6) ):
            Quit("MAINDISKA: maindisk.amin out of range.")
        if( (maindisk.amax <0.0) or (maindisk.amax > 0.6) ):
            Quit("MAINDISKA: maindisk.amax out of range.")
        if( maindisk.amax <= maindisk.amin ):
            Quit("MAINDISKA: maindisk.amax must be greater than maindisk.amin.");
        if( (maindisk.Hmax <= 0.0) or (maindisk.Hmax > maindisk.amax) ):
            Quit("MAINDISKH: maindisk.Hmax must be gt 0.0 and le maindisk.amax.")
        if( (maindisk.Hpow < 0.0) or (maindisk.Hpow > 2.0) ):
            Quit("MAINDISKH: maindisk.Hpow must be between 0.0 and 2.0.")
        if( (maindisk.maindiskL < 1.0e28) or (maindisk.maindiskL > 1.0e39)):
            Quit("MAINDISKT: maindisk luminosity must be between 1.0e28 and 1.0e39");
        if( maindisk.Ttype == "POWERLAW"):
            if( (maindisk.Tpow < -3.0) or (maindisk.Tpow > 3.0) ):
                Quit("MAINDISKT: maindisk.Tpow must be between -3.0 and 3.0.")
        if( (diskedge.T < 0.0) or (diskedge.T > 1.0e6) ):
            Quit("DISKEDGET: Edge T out of range.")
        if( (diskedge.Tspot < 0.0) or (diskedge.Tspot > 1.0e6) ):
            Quit("DISKEDGET: Tspot out of range.")
        if( (diskedge.ZetaMid < 0.0) or (diskedge.ZetaMid >= 360.0) ):
            Quit("DISKEDGET: ZetaMid must be ge 0.0 and lt 360.0.")
        if( (diskedge.ZetaWidth < 0.0) or (diskedge.ZetaWidth >= 360.0) ):
            Quit("DISKEDGET: ZetaWidth must be ge 0.0 and lt 360.0.")

    if( control.innerdisk == "ON"):
        if( (innerdisk.T < 0.0) or (innerdisk.T > 1.0e7) ):
	      Quit("Inner disk temperature out of range.")
        if( (innerdisk.L < 500.0) or (innerdisk.L > 1.0e39) ):
	      Quit("Inner disk luminosity out of range.")

    if( control.diskrim == "ON"):
        if( diskrim.awidth <= 0.0 ):
            Quit("DISKRIMAWIDTH is out of range.")
        if( diskrim.awidth > (maindisk.amax - maindisk.amin) ):
            Quit("DISKRIMAWIDTH is greater than the disk width.")
        if( diskrim.type == "SINUSOID"):
            if( (diskrim.Hmax < 0.0) or (diskrim.Hmax > maindisk.amax) ):
                Quit("DISKRIMPARS: Hmax must be ge 0.0 and le maindisk.amax.")
            if( (diskrim.Hmin < 0.0) or (diskrim.Hmin > maindisk.amax) ):
                Quit("DISKRIMPARS: Hmin must be ge 0.0 and le maindisk.amax.")
            if( diskrim.Hmax < diskrim.Hmin ):
                Quit("DISKRIMPARS: Hmax must be ge Hmin.")
            if( (diskrim.ZetaHmax < 0.0) or (diskrim.ZetaHmax >= 360.0) ):
                Quit("DISKRIMPARS: ZetaHmax must be ge 0.0 and lt 360.0.")
            if( (diskrim.Tmax < 0.0) or (diskrim.Tmax > 1.0e6) ):
                Quit("DISKRIMPARS: Tmax out of range.")
            if( (diskrim.Tmin < 0.0) or (diskrim.Tmin > 1.0e6) ):
                Quit("DISKRIMPARS: Tmin out of range.")
            if( diskrim.Tmax < diskrim.Tmin ):
                Quit("DISKRIMPARS: Tmax must be ge than Tmin.")
            if( (diskrim.ZetaTmax < 0.0) or (diskrim.ZetaTmax >= 360.0) ): 
	           Quit("DISKRIMPARS: ZetaTmax must be ge 0.0 and lt 360.0.")
        if( diskrim.type == "POINT"):
            for i in range(1, diskrim.points):
                if( (diskrim.PointZeta[i] < 0.0) or (diskrim.PointZeta[i] >= 360.0) ):
	               Quit("DISKRIMPARS: The PointZetas must be ge 0.0 and lt 360.0.")
                if( (diskrim.PointH[i] < 0.0) or (diskrim.PointH[i] > maindisk.amax) ):
 	               Quit("DISKRIMPARS: Rim H must be ge 0 and le maindisk.amax.")
                if( (diskrim.PointT[i] < 0.0) or (diskrim.PointT[i] > 1.0e6) ):
                    Quit("DISKRIMPARS: At least one rim T is out of range.")
        else:
            Quit("DISKRIMPARS: Unrecognized type.")

    if( control.disktorus == "ON"):
        if( (disktorus.azero >= maindisk.amax) or (disktorus.azero <= maindisk.amin) ):
            Quit("DISKTORUSAZERO is outside the disk.")
        if( (disktorus.azero - 0.5 * disktorus.awidth) < maindisk.amin):
            Quit("DISKTORUSAWIDTH: torus extends past the disk inner edge.");
        if( (disktorus.azero + 0.5 * disktorus.awidth) > maindisk.amax):
            Quit("DISKTORUSAWIDTH: torus extends past the disk outer edge.")
        if( control.diskrim == "ON"):
            if( (maindisk.amax - diskrim.awidth) < (disktorus.azero + 0.5 * disktorus.awidth) ):
                Quit("Disk rim and disk torus overlap.")
        if( disktorus.type == "SINUSOID"):
            if( (disktorus.Hmax < 0.0) or (disktorus.Hmax > maindisk.amax) ):
                Quit("DISKTORUSPARS: Hmax must be ge 0.0 and le maindisk.amax.")
            if( (disktorus.Hmin < 0.0) or (disktorus.Hmin > maindisk.amax) ):
                Quit("DISKTORUSPARS: Hmin must be ge 0.0 and le maindisk.amax.")
            if( disktorus.Hmax < disktorus.Hmin ):
                Quit("DISKTORUSPARS: Hmax must be ge Hmin.")
            if( (disktorus.ZetaHmax < 0.0) or (disktorus.ZetaHmax >= 360.0) ):
                Quit("DISKTORUSPARS: ZetaHmax must be ge 0.0 and lt 360.0.")
            if( (disktorus.Tmax < 0.0) or (disktorus.Tmax > 1.0e6) ):
                Quit("DISKTORUSPARS: Tmax out of range.")
            if( (disktorus.Tmin < 0.0) or (disktorus.Tmin > 1.0e6) ):
                Quit("DISKTORUSPARS: Tmin out of range.")
            if( disktorus.Tmax < disktorus.Tmin ):
                Quit("DISKTORUSPARS: Tmax must be ge than Tmin.")
            if( (disktorus.ZetaTmax < 0.0) or (disktorus.ZetaTmax >= 360.0) ): 
	           Quit("DISKTORUSPARS: ZetaTmax must be ge 0.0 and lt 360.0.")
        if( disktorus.type == "POINT"):
            for i in range(1, disktorus.points):
                if( (disktorus.PointZeta[i] < 0.0) or (disktorus.PointZeta[i] >= 360.0) ):
                    Quit("DISKTORUSPARS: The PointZetas must be ge 0.0 and lt 360.0.")
                if( (disktorus.PointH[i] < 0.0) or (disktorus.PointH[i] > maindisk.amax) ):
 	               Quit("DISKTORUSPARS: Torus H must be ge 0 and le maindisk.amax.")
                if( (disktorus.PointT[i] < 0.0) or (disktorus.PointT[i] > 1.0e6) ):
                     Quit("DISKTORUSPARS: At least one torus T is out of range.")
        else:
            Quit("DISKTORUSPARS: Unrecognized type.")

    if( control.diskspots == "ON"):
        if( diskspot.nspots <= 0): 
            Quit("DISKSPOTS= ON, but no spots specified.")
        for i in range(1, diskspot.nspots):
            if( (diskspot.zetamin[i] < 0.0) or (diskspot.zetamin[i] > 360.0) ):
                Quit("diskspot.zetamin out of range.")
            if( (diskspot.zetamax[i] < 0.0) or (diskspot.zetamin[i] > 360.0) ):
                Quit("diskspot.zetamax out of range.")
            if( diskspot.amin[i] >= diskspot.amax[i] ):
                Quit("diskspot.amax must be greater than diskspot.amin.")
            if( (diskspot.amin[i] < 0.0) or (diskspot.amin[i] > 0.6) ):
                Quit("diskspot.amin out of range.")
            if( (diskspot.amax[i] < 0.0) or (diskspot.amax[i] > 0.6) ):
                Quit("diskspot.amin out of range.")
            if( (diskspot.spotToverT[i] < 0.0) or (diskspot.spotToverT[i] > 100.0) ):
                Quit("diskspot.spotToverT out of range.")

    if( control.adc == "ON"):
        if( (adc.L < 0.0) or (adc.L > 1.0e40) ):
            Quit("Luminosity of the ADC is out of range.")
        if( (adc.height <= 0.0) or (adc.height > 0.5) ):
            Quit("Height of the point-approx ADC is out of range.")

    if( control.thirdlight == "ON" ):
        if( (thirdlight.orbphase < -0.5) or (thirdlight.orbphase >= 1.0) ):
            Quit("3rdLIGHTPHASE out of range.")
        if( thirdlight.nbands <= 0 ):
            Quit("Third light is on but no 3rdLIGHTFRACTIONs specified.")
        for i in range(1, thirdlight.nbands):
            if( orbit.filter[band] == "SQUARE"):
                if( (thirdlight.minlambda[i] < 0.0) or (thirdlight.minlambda[i] > 30000.) ):
                    Quit("One of the 3rd light SQUARE minlambdas is out of range.")
                if( (thirdlight.maxlambda[i] < 0.0) or (thirdlight.maxlambda[i] > 30000.) ):
	               Quit("One of the 3rd light SQUARE maxlambdas is out of range.")
                if( thirdlight.maxlambda[i] <= orbit.minlambda[i] ):
	               Quit("The 3rd light orbit.minlambda must be le than orbit.maxlambda.")
            if( (thirdlight.fraction[i] < 0.0) or (thirdlight.fraction[i] >= 1.0) ):
                Quit("3rd light fraction must be ge 0.0 and lt 1.0.")

    if( verbose != "ON") and (verbose !="OFF"):
        Quit("VERBOSE must be ON or OFF.")
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
        Quit("Cannot open file GDTable.dat.")

    i = -1
    while(fgets(inputline, 80, in) != NULL):
        if( inputline[0] != '*' ):
            i = i +1
            sscanf(inputline, "%lf %lf", &GDT[i], &fourbeta[i] )
    maxGDindex = i
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
    if((out == None):
        Quit("Cannot open file LDTable.dat.")

   for(;;) {
      fgets(inputline, 80, out)
      if( inputline[0] != '*' ) {
            sscanf( inputline, "%lf %lf %lf", &Tmin, &Tmax, &deltaT);
            maxLDTindex = ( Tmax - Tmin + 0.1 ) / deltaT;
            for( i = 0; i <= maxLDTindex; i++) {
	       LDT[i] = Tmin + i * deltaT;
            }
         fgets(inputline, 80, in);
            sscanf( inputline, "%lf %lf %lf", &gmin, &gmax, &deltag);
            maxLDgindex = ( gmax - gmin + 0.001) / deltag;
            for( i = 0; i <= maxLDgindex; i++) {
	       LDlogg[i] = gmin + i * deltag;
            }
         fgets(inputline, 80, in);
            sscanf( inputline, "%ld %s %s %s %s %s %s %s %s %s %s %s %s", 
		      &nfilters,
                      LDfilterName[0], LDfilterName[1],  LDfilterName[2],
                      LDfilterName[3], LDfilterName[4],  LDfilterName[5],
                      LDfilterName[6], LDfilterName[7],  LDfilterName[8],
		      LDfilterName[9], LDfilterName[10], LDfilterName[11]);
            maxLDfilterindex = nfilters - 1;
            if( maxLDfilterindex > (MAXFILTERS - 1) )
	      Quit("Too many filters in the LD table.");
         break;
      }
   }

   while( fscanf( in , "%lf %lf %s %lf %lf %lf %lf",
	  &logg, &T, filtername, &a[1], &a[2], &a[3], &a[4] ) != EOF ) {
      gindex = (logg - gmin + 0.1) / deltag;
      if( (gindex < 0) || (gindex > maxLDgindex) )
         Quit("gindex out of range in ReadLDTable.");
      Tindex = ( T - Tmin + 1.0) / deltaT;
      if( (Tindex < 0) || (Tindex > maxLDTindex) )
         Quit("Tindex out of range in ReadLDTable.");
      findex = -1;
      for( i = 0; i <= maxLDfilterindex; i++ ) {
         if( strcmp( filtername, LDfilterName[i]) == 0 ) 
            findex = i;
      }
      if( (findex < 0) || (findex > maxLDfilterindex) )
         Quit("Unrecognized filter name in ReadLDTable.");
      for( i = 1; i <= 4; i++ )
            LDtable[gindex][Tindex][findex][i] = a[i];
   }
   fclose ( in );
   return;
}


void ReadIperpTable()
/*****************************************************************
*
*   This function reads IperpTable.dat.   The file containing the Iperp
*   data must have the format:
*       Tmin    Tmax    deltaT    =  The minimum and maximum temperature
*                                    and the temperature spacing.
*       gmin    gmax    deltag    =  The minimum and maximum LOG g
*                                    and the spacing in log g.
*         n    filter1  filter2 filter3 filter4 ........
*        t   logg    Iz(lambda1)  Iz(lambda2)  Iz(lambda3)  Iz(lambda4)......
*        .     .           .            .            .            .
*        .     .           .            .            .            .
*        .     .           .            .            .            .
*
*    n is the number of filters.  "filter1" "filter2", etc. are the
*      names of the filters, eg U B V R I
*
*   The Iperp data are stored in the global variables
*      long maxIperpgindex, maxIperpTindex, maxIperplindex
*      double Iperplogg[gindex], IperpT[Tindex], Iperplambda[lindex]
*      double Iperptable[gindex][Tindex][lindex]
*  
*******************************************************************/
{
   FILE *in;
   char filename[40], inputline[201], Iperpfilter[15][20];
   long Tindex, gindex, findex, nfilters;
   double Tmin, Tmax, deltaT, gmin, gmax, deltag, xlogg, xT, xIperp[20];

   strcpy(filename, "IperpTable.dat");
   if((in = fopen(filename, "r")) == NULL)
      Quit("Cannot open file IperpTable.dat.");

   fscanf( in, "%lf  %lf  %lf", &Tmin, &Tmax, &deltaT);
      maxIperpTindex = ( Tmax - Tmin + 1.0 ) / deltaT;
      for( Tindex = 0; Tindex <= maxIperpTindex; Tindex++) {
	 IperpT[Tindex] = Tmin + Tindex * deltaT;
      }
   fscanf( in, "%lf  %lf  %lf", &gmin, &gmax, &deltag);
      maxIperpgindex = ( gmax - gmin + 0.01) / deltag;
      for( gindex = 0; gindex <= maxIperpgindex; gindex++) {
         Iperplogg[gindex] = gmin + gindex * deltag;
      }
   fgets(inputline, 120, in);
   fgets(inputline, 120, in);
   sscanf( inputline, "%ld %s %s %s %s %s %s %s %s %s %s %s %s",
            &nfilters,
            IperpfilterName[0],  IperpfilterName[1],  IperpfilterName[2], 
            IperpfilterName[3],  IperpfilterName[4],  IperpfilterName[5],
	    IperpfilterName[6],  IperpfilterName[7],  IperpfilterName[8],
	    IperpfilterName[9], IperpfilterName[10], IperpfilterName[11] );
   maxIperpfilterindex = nfilters - 1;
   if( maxIperpfilterindex > (MAXFILTERS - 1) )
     Quit("Too many filters in the Iperp table.");

   while( fgets( inputline, 200, in) != NULL ) {
      sscanf( inputline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	        &xT, &xlogg, 
                &xIperp[0],  &xIperp[1],  &xIperp[2], 
                &xIperp[3],  &xIperp[4],  &xIperp[5],
	        &xIperp[6],  &xIperp[7],  &xIperp[8],
	        &xIperp[9],  &xIperp[10], &xIperp[11]);
      gindex = (xlogg - gmin + 0.01) / deltag;
      if( (gindex < 0) || (gindex > maxIperpgindex) )
         Quit("gindex out of range in ReadIperpTable.");
      Tindex = ( xT - Tmin + 1.0) / deltaT;
      if( (Tindex < 0) || (Tindex > maxIperpTindex) )
         Quit("Tindex out of range in ReadIperpTable.");
      for( findex = 0; findex <= maxIperpfilterindex; findex++) {
         Iperptable[gindex][Tindex][findex] = xIperp[findex];
      }
   }
   fclose ( in );

   return;
}


void ReadIBBfilterTable( void )
/*****************************************************************
*
*   This function reads the intensities of a black body
*   observed through a filter.
*
*   The file name must be "IBBfilter.dat" and must have the format:
*       
*       Tmin    Tmax    deltaT    =  The minimum and maximum temperature
*                                    and the temperature spacing.
*         N filtername1 filtername2 filtername3 ... filternameN
*        T   Ifilter1  Ifilter2  Ifilter3 ... IfilterN
*        .      .         .         .            .
*        .      .         .         .            .
*        .      .         .         .            .
*
*   The filter can have comment lines beginning with a "*" at the
*   beginning.
*
*******************************************************************/
{
   FILE *in;
   char filename[40], inputline[201];
   long i, nfilters, Tindex, findex;
   double xT, xIBB[20];

   strcpy(filename, "IBBfilterTable.dat");
   if((in = fopen(filename, "r")) == NULL)
      Quit("Cannot open file IBBfilter.dat.");

   /******************************************************
   *
   *   Skip over any initial comment lines (lines beginning
   *   with a '*') then read two lines containing the
   *   information about the contents of the table.
   *
   ********************************************************/
   for(;;) {
      fgets(inputline, 80, in);
      if( inputline[0] != '*' ) {
            sscanf( inputline, "%lf %lf %lf", 
                                &IBBTmin, &IBBTmax, &IBBdeltaT);   
            maxIBBTindex = ( IBBTmax - IBBTmin + 0.1 ) / IBBdeltaT;
            for( i = 0; i <= maxIBBTindex; i++) {
               IBBT[i] = IBBTmin + i * IBBdeltaT;
            }
         fgets(inputline, 80, in);
            sscanf( inputline, "%ld %s %s %s %s %s %s %s %s %s %s %s %s", 
	           &nfilters,
                   IBBfilterName[0], IBBfilterName[1],  IBBfilterName[2],
                   IBBfilterName[3], IBBfilterName[4],  IBBfilterName[5],
                   IBBfilterName[6], IBBfilterName[7],  IBBfilterName[8],
	           IBBfilterName[9], IBBfilterName[10], IBBfilterName[11]);
            maxIBBfilterindex = nfilters - 1;
            if( maxIBBfilterindex > (MAXFILTERS - 1) )
	      Quit("Too many filters in the IBBfilter table.");
         break;
      }
   }
   while( fgets( inputline, 200, in) != NULL ) {
      sscanf( inputline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	        &xT,
                &xIBB[0], &xIBB[1], &xIBB[2], &xIBB[3], &xIBB[4],  &xIBB[5],
	        &xIBB[6], &xIBB[7], &xIBB[8], &xIBB[9], &xIBB[10], &xIBB[11]);
      Tindex = ( xT - IBBTmin + 1.0) / IBBdeltaT;
      if( (Tindex < 0) || (Tindex > maxIBBTindex) )
         Quit("Tindex out of range in ReadIBBTable.");
      for( findex = 0; findex <= maxIBBfilterindex; findex++) {
         IBBtable[Tindex][findex] = xIBB[findex];
      }
   }
   fclose ( in );

   return;
}


void ReadZzetaTable( void )
/*******************************************************
*
*   Read the ZBBzeta table
*
***********************************************************/
{
   FILE *in;
   char filename[40];
   long i;
   double dummy;

   strcpy( filename, "ZzetaTable.dat");
   if((in = fopen(filename, "r")) == NULL)
      Quit("Cannot open file ZzetaTable.dat.");

   fscanf( in, "%ld %lf", &maxBBzetaindex, &deltaBBzeta);
   for( i = 0; i <= maxBBzetaindex; i++) {
     fscanf( in, "%lf %lf", &dummy, &ZBBzeta[i]);
   }
   BBzetamax = maxBBzetaindex * deltaBBzeta;

   fclose( in );

   return;
}

