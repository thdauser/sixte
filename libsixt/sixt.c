#include "sixt.h"
#include "headas_rand.h"


#ifdef USE_RCL
// Use the Remeis random number server.
#include <rcl.h>
#endif


/** MJDREF used in the FITS header of eROSITA event files
    [d]. Corresponds to 2000-01-01T00:00:00. */
const double eromjdref=51543.875;

/** MJDREF used in the FITS header of XMM event files [d]. Corresponds
    to 1998-01-01T00:00:00.00. */
const double xmmmjdref=50814.;


unsigned int getSeed(int seed)
{
  if (seed>=0) {
    return((unsigned int)seed);
  } else {
    // Determine the seed from the system clock.
    return((unsigned int)time(NULL));
  }
}


double sixt_get_random_number(int* const status)
{
  // Return a random value out of the interval [0,1).

#ifdef USE_RCL

  // Use the Remeis random number server running on leo.
  // The library for accessing the server is maintained
  // by Fritz-Walter Schwarm.
  int rcl_status=0;
  double rand=rcl_rand_net_get_double(NULL, NULL, &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

  return(rand);  
#else

  // Use the HEAdas random number generator.
  return(HDmtDrand());

  // Status variable is not needed.
  (void)(*status);
#endif
}


void sixt_init_rng(const unsigned int seed, int* const status)
{
  // Initialize HEAdasS random number generator.
  // Note that this has to be done in any case, even
  // if the RCL random number server is used, because
  // the HEAdas routines (like heasp) rely on HDmtDrand().
  HDmtInit(seed);

#ifdef USE_RCL

  // Call the RCL random number generator specifying the
  // server and method.
  int rcl_status=0;
  rcl_rand_net_get_double("draco", "rand", &rcl_status);

  if(RCL_RANDOM_SUCCESS!=rcl_status) {
    SIXT_ERROR("failed getting random number from RCL");
    *status=EXIT_FAILURE;
  }

#else 

  // The status variable is not used.
  (void)(*status);

#endif
}


void sixt_destroy_rng()
{
  // Release HEADAS random number generator:
  HDmtFree();
}


void sixt_get_gauss_random_numbers(double* const x, 
				   double* const y, 
				   int* const status)
{
  double sqrt_2rho=sqrt(-log(sixt_get_random_number(status))*2.);
  CHECK_STATUS_VOID(*status);
  double phi=sixt_get_random_number(status)*2.*M_PI;
  CHECK_STATUS_VOID(*status);

  *x=sqrt_2rho * cos(phi);
  *y=sqrt_2rho * sin(phi);
}


double rndexp(const double avgdist, int* const status)
{
  assert(avgdist>0.);

  double rand=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0.);
  if (rand<1.e-15) {
    rand=1.e-15;
  }

  return(-log(rand)*avgdist);
}


void strtoupper(char* const string) 
{
  int count=0;
  while (string[count]!='\0') {
    string[count]=toupper(string[count]);
    count++;
  }
}


void sixt_error(const char* const func, const char* const msg)
{
  // Use the HEADAS error output routine.
  char output[MAXMSG];
  sprintf(output, "Error in %s: %s!\n", func, msg);
  HD_ERROR_THROW(output, EXIT_FAILURE);
}


void sixt_warning(const char* const msg)
{
  // Print the formatted output message.
  headas_chat(1, "### Warning: %s!\n", msg);
}


void sixt_get_XMLFile(char* const filename,
		      const char* const xmlfile,
		      const char* const mission,
		      const char* const instrument,
		      const char* const mode,
		      int* const status)
{
  // Convert the user input to capital letters.
  char Mission[MAXMSG];
  char Instrument[MAXMSG];
  char Mode[MAXMSG];
  strcpy(Mission, mission);
  strcpy(Instrument, instrument);
  strcpy(Mode, mode);
  strtoupper(Mission);
  strtoupper(Instrument);
  strtoupper(Mode);

  // Check the available missions, instruments, and modes.
  char XMLFile[MAXFILENAME];
  strcpy(XMLFile, xmlfile);
  strtoupper(XMLFile);
  if (0==strcmp(XMLFile, "NONE")) {
    // Determine the base directory containing the XML
    // definition files.
    strcpy(filename, SIXT_DATA_PATH);
    strcat(filename, "/instruments");

    // Determine the XML filename according to the selected
    // mission, instrument, and mode.
    if (0==strcmp(Mission, "SRG")) {
      strcat(filename, "/srg");
      if (0==strcmp(Instrument, "EROSITA")) {
	strcat(filename, "/erosita.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "IXO")) {
      strcat(filename, "/ixo");
      if (0==strcmp(Instrument, "WFI")) {
	strcat(filename, "/wfi");
	if (0==strcmp(Mode, "FULLFRAME")) {
	  strcat(filename, "/fullframe.xml");
	} else {
	  *status=EXIT_FAILURE;
	  SIXT_ERROR("selected mode is not supported");
	  return;
	}
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "ATHENA")) {
      strcat(filename, "/athena");
      if (0==strcmp(Instrument, "WFI")) {
	strcat(filename, "/wfi");
	if (0==strcmp(Mode, "FULLFRAME")) {
	  strcat(filename, "/fullframe.xml");
	} else {
	  *status=EXIT_FAILURE;
	  SIXT_ERROR("selected mode is not supported");
	  return;
	}
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else if (0==strcmp(Mission, "GRAVITAS")) {
      strcat(filename, "/gravitas");
      if (0==strcmp(Instrument, "HIFI")) {
	strcat(filename, "/hifi.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }
      
    } else if (0==strcmp(Mission, "NUSTAR")) {
      strcat(filename, "/nustar");
      if (0==strcmp(Instrument, "NUSTAR")) {
	strcat(filename, "/nustar.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected instrument is not supported");
	return;
      }

    } else {
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "selected mission ('%s') is not supported", Mission);
      SIXT_ERROR(msg);
      return;
    }
    
  } else {
    // The XML filename has been given explicitly.
    strcpy(filename, xmlfile);
  }
}


void sixt_get_LADXMLFile(char* const filename,
			 const char* const xmlfile)
{
  // Check the available missions, instruments, and modes.
  char XMLFile[MAXFILENAME];
  strcpy(XMLFile, xmlfile);
  strtoupper(XMLFile);
  if (0==strcmp(XMLFile, "NONE")) {
    // Set default LAD XML file.
    strcpy(filename, SIXT_DATA_PATH);
    strcat(filename, "/instruments/loft/lad.xml");
    
  } else {
    // The XML filename has been given explicitly.
    strcpy(filename, xmlfile);
  }
}


void sixt_get_date_time(const double mjdref,
			const double t,
			char* const datestr,
			char* const timestr,
			int* const status)
{
  int int_day=(int)(mjdref-40587.0+t/86400.0);
  int int_sec=(int)((mjdref-40587.0-int_day)*86400.0+t);
  struct tm time_utc;
  time_utc.tm_sec=int_sec;
  time_utc.tm_min=0;
  time_utc.tm_hour=0;
  time_utc.tm_mday=1+int_day;
  time_utc.tm_mon=0;
  time_utc.tm_year=70;
  // avoid warning by valgrind
  time_utc.tm_isdst=-1;
  time_t timet=mktime(&time_utc);
  // Note that we have to use 'localtime' here, although we
  // want to determine the UTC, because this is the inverse
  // of 'mktime'.
  struct tm* time_utcn=localtime(&timet);
  if (NULL==time_utcn) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not determine UTC time");
    return;
  }
  if (10!=strftime(datestr, MAXMSG, "%Y-%m-%d", time_utcn)) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("failed formatting date string");
    return;
  }
  if (8!=strftime(timestr, MAXMSG, "%H:%M:%S", time_utcn)) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("failed formatting time string");
    return;
  }
}


void sixt_add_fits_stdkeywords(fitsfile* const fptr,
			       const int hdunum,
			       char* const telescop,
			       char* const instrume,
			       char* const filter,
			       char* const ancrfile,
			       char* const respfile,
			       double mjdref,
			       double timezero,
			       double tstart,
			       double tstop,
			       int* const status)
{
  // Determine the current HDU.
  int prev_hdunum=0;
  fits_get_hdu_num(fptr, &prev_hdunum);

  // Move to the desired HDU.
  if (prev_hdunum!=hdunum) {
    int hdutype=0;
    fits_movabs_hdu(fptr, hdunum, &hdutype, status);
    CHECK_STATUS_VOID(*status);
  }

  // Set the mission keywords.
  fits_update_key(fptr, TSTRING, "TELESCOP", telescop, "", status);
  fits_update_key(fptr, TSTRING, "INSTRUME", instrume, "", status);
  fits_update_key(fptr, TSTRING, "FILTER", filter, "", status);
  CHECK_STATUS_VOID(*status);

  // Set the names of the response files used for the simulation.
  fits_update_key(fptr, TSTRING, "ANCRFILE", ancrfile, 
		  "ancillary response file", status);
  fits_update_key(fptr, TSTRING, "RESPFILE", respfile, 
		  "response file", status);
  CHECK_STATUS_VOID(*status);

  // Set the timing keywords.
  // Determine the current date and time (Stevens, "Advanced Programming 
  // in the UNIX environment", p. 155 ff.).
  time_t current_time;
  if (0==time(&current_time)) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not determine current time");
    return;
  }
  struct tm* current_time_utc=gmtime(&current_time);
  if (NULL==current_time_utc) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not determine current UTC time");
    return;
  }   
  char current_datetimestr[MAXMSG];
  if (19!=strftime(current_datetimestr, MAXMSG, "%Y-%m-%dT%H:%M:%S", 
		   current_time_utc)) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("failed formatting date-time string");
    return;
  }
  fits_update_key(fptr, TSTRING, "DATE", current_datetimestr, 
		  "file creation date", status);
  CHECK_STATUS_VOID(*status);

  // Determine the start date and time.
  char datestr[MAXMSG], timestr[MAXMSG];
  sixt_get_date_time(mjdref, tstart, datestr, timestr, status);
  CHECK_STATUS_VOID(*status);
  fits_update_key(fptr, TSTRING, "DATE-OBS", datestr, 
		  "UT date of observation start", status);
  fits_update_key(fptr, TSTRING, "TIME-OBS", timestr, 
		  "UT time of observation start", status);
  CHECK_STATUS_VOID(*status);

  // Determine the stop date and time.
  sixt_get_date_time(mjdref, tstop, datestr, timestr, status);
  CHECK_STATUS_VOID(*status);
  fits_update_key(fptr, TSTRING, "DATE-END", datestr, 
		  "UT date of observation end", status);
  fits_update_key(fptr, TSTRING, "TIME-END", timestr, 
		  "UT time of observation end", status);
  CHECK_STATUS_VOID(*status);

  // MJDREF, TSTART, TSTOP.
  fits_update_key(fptr, TDOUBLE, "MJDREF", &mjdref, 
		  "reference MJD", status);
  fits_update_key(fptr, TDOUBLE, "TIMEZERO", &timezero,
		  "time offset", status);
  fits_update_key(fptr, TDOUBLE, "TSTART", &tstart,
		  "start time", status);
  fits_update_key(fptr, TDOUBLE, "TSTOP", &tstop,
		  "stop time", status);
  CHECK_STATUS_VOID(*status);


  // Add header information about program parameters if
  // this is the primary extension.
  if (1==hdunum) {
    HDpar_stamp(fptr, hdunum, status);
    CHECK_STATUS_VOID(*status);
  }

  // Move back to the original HDU.
  if (prev_hdunum!=hdunum) {
    int hdutype=0;
    fits_movabs_hdu(fptr, prev_hdunum, &hdutype, status);
    CHECK_STATUS_VOID(*status);
  }
}


void sixt_read_fits_stdkeywords(fitsfile* const ifptr,
			       char* const telescop,
			       char* const instrume,
			       char* const filter,
			       char* const ancrfile,
			       char* const respfile,
			       double *mjdref,
			       double *timezero,
			       double *tstart,
			       double *tstop, 
				int* const status)
{
  
  char comment[MAXMSG];
  
  fits_read_key(ifptr, TSTRING, "TELESCOP", telescop, comment, status);
  fits_read_key(ifptr, TSTRING, "INSTRUME", instrume, comment, status);
  fits_read_key(ifptr, TSTRING, "FILTER", filter, comment, status);
  CHECK_STATUS_VOID(*status);
  fits_read_key(ifptr, TSTRING, "ANCRFILE", ancrfile, comment, status);
  fits_read_key(ifptr, TSTRING, "RESPFILE", respfile, comment, status);
  CHECK_STATUS_VOID(*status);
  // MJDREF, TSTART, TSTOP.
  fits_read_key(ifptr, TDOUBLE, "MJDREF", mjdref, comment, status);
  fits_read_key(ifptr, TDOUBLE, "TIMEZERO", timezero, comment, status);
  fits_read_key(ifptr, TDOUBLE, "TSTART", tstart, comment, status);
  fits_read_key(ifptr, TDOUBLE, "TSTOP", tstop, comment, status);
  CHECK_STATUS_VOID(*status);
}


void sixt_add_fits_erostdkeywords(fitsfile* const fptr,
				  const int hdunum,
				  char* const creation_date,
				  char* const date_obs,
				  char* const time_obs,
				  char* const date_end,
				  char* const time_end,
				  double tstart,
				  double tstop,
				  double mjdref,
				  double timezero,
				  int ccdnr,
				  int* const status)
{
  // Determine the current HDU.
  int prev_hdunum=0;
  fits_get_hdu_num(fptr, &prev_hdunum);

  // Move to the desired HDU.
  if (prev_hdunum!=hdunum) {
    int hdutype=0;
    fits_movabs_hdu(fptr, hdunum, &hdutype, status);
    CHECK_STATUS_VOID(*status);
  }

  char origin[MAXMSG]="ECAP";
  fits_update_key(fptr, TSTRING, "ORIGIN", origin, "Origin of FITS file", status);
  char creator[MAXMSG]="SIXTE";
  fits_update_key(fptr, TSTRING, "CREATOR", creator, "Program that created this FITS file", status);
  char mission[MAXMSG]="SRG";
  fits_update_key(fptr, TSTRING, "MISSION", mission, "", status);
  char telescop[MAXMSG]="eROSITA";
  fits_update_key(fptr, TSTRING, "TELESCOP", telescop, "", status);
  char instrume[MAXMSG];
  sprintf(instrume, "FM%d", ccdnr);
  fits_update_key(fptr, TSTRING, "INSTRUME", instrume, "", status);
  fits_update_key(fptr, TSTRING, "INSTRUM1", instrume, "", status);
  int ninst=1;
  fits_update_key(fptr, TINT, "NINST", &ninst, "", status);

  char obsmode[MAXMSG]="SURVEY";
  fits_update_key(fptr, TSTRING, "OBS_MODE", obsmode, "", status);
  char datamode[MAXMSG]="";
  fits_update_key(fptr, TSTRING, "DATAMODE", datamode, "", status);
  
  float frametim=50.0;
  fits_update_key(fptr, TFLOAT, "FRAMETIM", &frametim, "[ms] nominal frame time", status);
  
  char filter[MAXMSG]="OPEN";
  fits_update_key(fptr, TSTRING, "FILTER", filter, "", status);

  long obs_id=0;
  fits_update_key(fptr, TLONG, "OBS_ID", &obs_id, "", status);
  long exp_id=0;
  fits_update_key(fptr, TLONG, "EXP_ID", &exp_id, "", status);

  char observer[MAXMSG]="";
  fits_update_key(fptr, TSTRING, "OBSERVER", observer, "", status);
  char object[MAXMSG]="";
  fits_update_key(fptr, TSTRING, "OBJECT", object, "", status);

  double ra_obj=0.0;
  fits_update_key(fptr, TDOUBLE, "RA_OBJ", &ra_obj, "[deg] J2000", status);
  double de_obj=0.0;
  fits_update_key(fptr, TDOUBLE, "DEC_OBJ", &de_obj, "[deg] J2000", status);

  fits_update_key(fptr, TSTRING, "DATE", creation_date, "File creation date", status);
  char date_obs_time[MAXMSG]="";
  sprintf(date_obs_time, "%sT%s", date_obs, time_obs);
  fits_update_key(fptr, TSTRING, "DATE-OBS", date_obs_time, "UT date of observation start", status);
  char date_end_time[MAXMSG]="";
  sprintf(date_end_time, "%sT%s", date_end, time_end);
  fits_update_key(fptr, TSTRING, "DATE-END", date_end_time, "UT date of observation end", status);

  fits_update_key(fptr, TDOUBLE, "TSTART", &tstart, "Start time of exposure in units of TIME column", status);
  fits_update_key(fptr, TDOUBLE, "TSTOP", &tstop, "Stop time of exposure in units of TIME column", status);

  fits_update_key(fptr, TDOUBLE, "MJDREF", &mjdref, "[d]", status);
  fits_update_key(fptr, TDOUBLE, "TIMEZERO", &timezero, "Time offset", status);
  
  char timeunit[MAXMSG]="s";
  fits_update_key(fptr, TSTRING, "TIMEUNIT", timeunit, "Time unit", status);
  char timesys[MAXMSG]="TT";
  fits_update_key(fptr, TSTRING, "TIMESYS", timesys, "Time system (Terrestial Time)", status);

  double ra_pnt=0.0;
  fits_update_key(fptr, TDOUBLE, "RA_PNT", &ra_pnt, "[deg] actual pointing RA J2000", status);
  double dec_pnt=0.0;
  fits_update_key(fptr, TDOUBLE, "DEC_PNT", &dec_pnt, "[deg] actual pointing DEC J2000", status);
  double pa_pnt=0.0;
  fits_update_key(fptr, TDOUBLE, "PA_PNT", &pa_pnt, "[deg] mean/median position angle of pointing", status);
  
  char radecsys[MAXMSG]="FK5";
  fits_update_key(fptr, TSTRING, "RADECSYS", radecsys, "Stellar reference frame", status);
  double equinox=2000.0;
  fits_update_key(fptr, TDOUBLE, "EQUINOX", &equinox, "Coordinate system equinox", status);

  char longstr[MAXMSG]="OGIP 1.0";
  fits_update_key(fptr, TSTRING, "LONGSTR", longstr, "", status);

  int ibuffer=384;
  fits_update_key(fptr, TINT, "NXDIM", &ibuffer, "", status);
  fits_update_key(fptr, TINT, "NYDIM", &ibuffer, "", status);
  float fbuffer=75.0;
  fits_update_key(fptr, TFLOAT, "PIXLEN_X", &fbuffer, "[micron]", status);
  fits_update_key(fptr, TFLOAT, "PIXLEN_Y", &fbuffer, "[micron]", status);
  
  CHECK_STATUS_VOID(*status);

  // Move back to the original HDU.
  if (prev_hdunum!=hdunum) {
    int hdutype=0;
    fits_movabs_hdu(fptr, prev_hdunum, &hdutype, status);
    CHECK_STATUS_VOID(*status);
  }
}


void verifyMJDREF(const double refmjdref,
		  const double mjdref,
		  const char* const description,
		  int* const status)
{
  if (fabs(mjdref-refmjdref)>1.e-6) {
    *status=EXIT_FAILURE;
    char insertmsg[MAXMSG];
    if (NULL==description) {
      strcpy(insertmsg, "");
    } else {
      strcpy(insertmsg, description);
      strcat(insertmsg, " ");
    }
    char msg[MAXMSG];
    sprintf(msg, "MJDREF %sdoes not match required value of '%.1lf'", 
	    insertmsg, refmjdref);
    SIXT_ERROR(msg);
  }
}


void verifyTIMEZERO(const double timezero,
		    int* const status)
{
  if (0.0!=timezero) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("implementation requires that TIMEZERO=0.0 in input file");
  }
}


float getEBOUNDSEnergy(const long channel,
		       const struct RMF* const rmf, 
		       int* const status)
{
  float lo, hi;
  getEBOUNDSEnergyLoHi(channel, rmf, &lo, &hi, status);
  CHECK_STATUS_RET(*status, 0.0);
  double r=sixt_get_random_number(status);
  CHECK_STATUS_RET(*status, 0.0);
  return(r*lo + (1.0-r)*hi);
}
