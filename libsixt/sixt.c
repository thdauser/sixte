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

/** FITS error message for error handling macro, which retrieves error status */
char _fits_err_msg[80];

unsigned long microtime(){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

unsigned int getSeed(int seed)
{
  if (seed>=0) {
    return((unsigned int)seed);
  } else {
    // Determine the seed from the system clock in MicroSecond precision.
    return((unsigned int)microtime());
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

void sixt_deprecated(const char* const fnc, const char* const alt)
{
  // Print a warning that this function is deprecated and propose the
  // alternative supplied in alt (if not NULL).
  headas_chat(1, "### Warning: The function %s is deprecated!\n", fnc);
  if(alt != NULL && strlen(alt) > 0) {
    headas_chat(1, "### Please consider using the function %s instead.\n", alt);
  }
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
  strncpy(Mission, mission,MAXMSG-1);
  Mission[MAXMSG-1]='\0';
  strncpy(Instrument, instrument,MAXMSG-1);
  Instrument[MAXMSG-1]='\0';
  strncpy(Mode, mode,MAXMSG-1);
  Mode[MAXMSG-1]='\0';
  strtoupper(Mission);
  strtoupper(Instrument);
  strtoupper(Mode);

  // Check the available missions, instruments, and modes.
  char XMLFile[MAXFILENAME];
  strncpy(XMLFile, xmlfile,MAXFILENAME-1);
  XMLFile[MAXFILENAME-1]='\0';
  strtoupper(XMLFile);

  // check whether XML filename has been given explicitly
  if (0!=strcmp(XMLFile,"NONE")) {
    // ... yes: just copy it and be happy
    strcpy(filename,xmlfile);
    return;
  }
    

  // Determine the base directory containing the XML
  // definition files.
  strcpy(filename, SIXT_DATA_PATH);
  strcat(filename, "/instruments");

  // Determine the XML filename according to the selected
  // mission, instrument, and mode.
  if (0==strcmp(Mission, "SRG")) {
    strcat(filename, "/srg");
    if (0==strcmp(Instrument, "EROSITA")) {
      strcat(filename, "/erosita_dummy.xml");
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected instrument is not supported");
    }
    return;
  }
 
  if (0==strcmp(Mission, "IXO")) {
    strcat(filename, "/ixo");
    if (0==strcmp(Instrument, "WFI")) {
      strcat(filename, "/wfi");
      if (0==strcmp(Mode, "FULLFRAME")) {
	strcat(filename, "/fullframe.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected mode is not supported");
      }
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected instrument is not supported");
    }
    return;
  }

  
  if (0==strcmp(Mission, "ATHENA")) {
    strcat(filename, "/athena");
    if (0==strcmp(Instrument, "WFI")) {
      strcat(filename, "/wfi");
      if (0==strcmp(Mode, "FULLFRAME")) {
	strcat(filename, "/fullframe.xml");
      } else {
	*status=EXIT_FAILURE;
	SIXT_ERROR("selected mode is not supported");
      }
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected instrument is not supported");
    }
    return;
  }

  if (0==strcmp(Mission, "GRAVITAS")) {
    strcat(filename, "/gravitas");
    if (0==strcmp(Instrument, "HIFI")) {
      strcat(filename, "/hifi.xml");
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected instrument is not supported");
    }
    return;
  }
    
  if (0==strcmp(Mission, "NUSTAR")) {
    strcat(filename, "/nustar");
    if (0==strcmp(Instrument, "NUSTAR")) {
      strcat(filename, "/nustar.xml");
    } else {
      *status=EXIT_FAILURE;
      SIXT_ERROR("selected instrument is not supported");
    }
    return;
  }
  
  *status=EXIT_FAILURE;
  char msg[MAXMSG];
  sprintf(msg, "selected mission ('%s') is not supported", Mission);
  SIXT_ERROR(msg);
  return;
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void sixt_get_eroXMLFile(char *filename,
			const int telescop_index,
			int* const status){
  
  sprintf(filename, "%s/instruments/srg/erosita_%d.xml", SIXT_DATA_PATH, telescop_index+1);
  
}
#pragma GCC diagnostic pop


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


void sixt_add_fits_stdkeywords_obsolete(fitsfile* const fptr,
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


void sixt_read_fits_stdkeywords_obsolete(fitsfile* const ifptr,
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
				  char* const filter,
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
  sprintf(instrume, "TM%d", ccdnr);
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
    sprintf(msg, "MJDREF %sdoes not match required value of '%.4lf'", 
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

/** Add standard FITS header keywords to the specified file using info
 *  contained in a SixtStdKeywords structure. */
void sixt_add_fits_stdkeywords(fitsfile* const fptr,
			       const int hdunum,
			       SixtStdKeywords * keyword_struct,
			       int* const status){
  // Run usual function
  sixt_add_fits_stdkeywords_obsolete(fptr,hdunum,keyword_struct->telescop,keyword_struct->instrume,keyword_struct->filter,
				     keyword_struct->ancrfile,keyword_struct->respfile,keyword_struct->mjdref,keyword_struct->timezero,
				     keyword_struct->tstart,keyword_struct->tstop,status);
  CHECK_STATUS_VOID(*status);
  
  // Add extname keyword if we are:
  //   * not on the primary header,
  //   * if the extname is non-NULL, and 
  //   * if no EXTNAME has been set
  if((hdunum>1) && (keyword_struct->extname!=NULL)){
    char keyval[FLEN_VALUE];
    fits_read_key(fptr,TSTRING,"EXTNAME",keyval,NULL,status);
    if (*status==VALUE_UNDEFINED) {
      *status=0;
      fits_update_key(fptr, TSTRING, "EXTNAME", keyword_struct->extname,
		      "Name of this extension", status);
      CHECK_STATUS_VOID(*status);
    }
  }
}


/** Reads standard header keywords from a FITS file using a
 *  SixtStdKeywords structure. Does at the same time the
 *  malloc of the different char arrays. */
void sixt_read_fits_stdkeywords(fitsfile* const ifptr,
		SixtStdKeywords* keyword_struct,
		int* const status){

	char comment[FLEN_COMMENT];
	char keyword[FLEN_KEYWORD];
	//Read string keywords
	fits_read_key(ifptr, TSTRING, "TELESCOP", keyword, comment, status);
	if(NULL!=keyword_struct->telescop){
		free(keyword_struct->telescop);
	}
	keyword_struct->telescop = strdup(keyword);
	fits_read_key(ifptr, TSTRING, "INSTRUME", keyword, comment, status);

	if(NULL!=keyword_struct->instrume){
		free(keyword_struct->instrume);
	}
	keyword_struct->instrume = strdup(keyword);

	fits_read_key(ifptr, TSTRING, "FILTER", keyword, comment, status);
	if(NULL!=keyword_struct->filter){
		free(keyword_struct->filter);
	}
	keyword_struct->filter = strdup(keyword);

	fits_read_key(ifptr, TSTRING, "ANCRFILE", keyword, comment, status);
	if(NULL!=keyword_struct->ancrfile){
		free(keyword_struct->ancrfile);
	}
	keyword_struct->ancrfile = strdup(keyword);

	fits_read_key(ifptr, TSTRING, "RESPFILE", keyword, comment, status);
	if(NULL!=keyword_struct->respfile){
		free(keyword_struct->respfile);
	}
	keyword_struct->respfile = strdup(keyword);
	CHECK_STATUS_VOID(*status);

	fits_read_key(ifptr, TDOUBLE, "MJDREF", &keyword_struct->mjdref, comment, status);
	fits_read_key(ifptr, TDOUBLE, "TIMEZERO", &keyword_struct->timezero, comment, status);
	fits_read_key(ifptr, TDOUBLE, "TSTART", &keyword_struct->tstart, comment, status);
	fits_read_key(ifptr, TDOUBLE, "TSTOP", &keyword_struct->tstop, comment, status);
	CHECK_STATUS_VOID(*status);

	// Determine the current HDU.
	int hdunum=0;
	fits_get_hdu_num(ifptr, &hdunum);
	CHECK_STATUS_VOID(*status);

	// If we are not on the primary HDU, read extname
	if(NULL!=keyword_struct->extname){
	  free(keyword_struct->extname);
	}
	if(hdunum>1){
		fits_read_key(ifptr,TSTRING, "EXTNAME", keyword, comment, status);
		keyword_struct->extname = strdup(keyword);
		CHECK_STATUS_VOID(*status);
	} else {
	  keyword_struct->extname=NULL;
	}
}

/** Constructor of the SixtStdKeywords structure: returns a pointer to an empty structure of this type */
SixtStdKeywords* newSixtStdKeywords(int* const status){
	SixtStdKeywords* keywords=malloc(sizeof(*keywords));
	if (NULL==keywords) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for SixtStdKeywords failed");
		return(keywords);
	}
	keywords->telescop=NULL;
	keywords->instrume=NULL;
	keywords->filter=NULL;
	keywords->ancrfile=NULL;
	keywords->respfile=NULL;
	keywords->extname=NULL;

	keywords->mjdref = 0;
	keywords->timezero = 0;
	keywords->tstart = 0;
	keywords->tstop = 0;

	return(keywords);
}


/** Builds a SixtStdKeywords struct from the individual keywords.
 * 	Does at the same time the malloc of the different char arrays. */
SixtStdKeywords* buildSixtStdKeywords(char* const telescop,
	       char* const instrume,
	       char* const filter,
	       char* const ancrfile,
	       char* const respfile,
	       char* const extname,
	       double mjdref,
	       double timezero,
	       double tstart,
	       double tstop,
	       int* const status){

	SixtStdKeywords* keywords = newSixtStdKeywords(status);
	CHECK_STATUS_RET(*status,keywords);

	keywords->telescop = strdup(telescop);
	keywords->instrume = strdup(instrume);
	keywords->filter = strdup(filter);
	keywords->ancrfile = strdup(ancrfile);
	keywords->respfile = strdup(respfile);
	if (extname!=NULL) {
	  keywords->extname = strdup(extname);
	} else {
	  keywords->extname=NULL;
	}

	keywords->mjdref = mjdref;
	keywords->timezero = timezero;
	keywords->tstart = tstart;
	keywords->tstop = tstop;

	return(keywords);
}

/** Duplicate a SixtStdKeywords structure **/
// we set the extname to NULL by default
SixtStdKeywords* duplicateSixtStdKeywords(const SixtStdKeywords *key,int* const status) {
  return buildSixtStdKeywords(key->telescop,
			      key->instrume,
			      key->filter,
			      key->ancrfile,
			      key->respfile,
			      NULL,
			      key->mjdref,
			      key->timezero,
			      key->tstart,
			      key->tstop,
			      status);
}


/** Destructor of the SixtStdKeywordsStructure */
void freeSixtStdKeywords(SixtStdKeywords* keyword_struct){
  if (keyword_struct==NULL) {
    return;
  }
  free(keyword_struct->telescop);
  free(keyword_struct->instrume);
  free(keyword_struct->filter);
  free(keyword_struct->ancrfile);
  free(keyword_struct->respfile);
  if (keyword_struct->extname!=NULL) {
    free(keyword_struct->extname);
  }
  free(keyword_struct);
  keyword_struct=NULL;
}

/** convenience function to create a FITS-file or to error out depending on the value of **/
/** clobber **/
int fits_create_file_clobber(fitsfile **fptr, char *filename, int clobber, int *status) {
  CHECK_STATUS_RET(*status,0);

  int exists;
  fits_file_exists(filename, &exists, status);
  if (0!=exists) {
    if (clobber) {
      // Delete the file.
      remove(filename);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      snprintf(msg,MAXMSG,"file '%s' already exists", filename);
      msg[MAXMSG-1]='\0';
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      CHECK_STATUS_RET(*status,0);
    }
  }
  int ret=fits_create_file(fptr, filename, status);
  return ret;
}

void fits_close_file_chksum(fitsfile *fptr,int *status) {
  // close a file after updating the checksum in the current HDU
  // and in the primary extension
  int iomode;
  fits_file_mode(fptr,&iomode,status);
  if (iomode==READWRITE){
    fits_write_chksum(fptr,status);
    fits_movabs_hdu(fptr,1,NULL,status);
    fits_write_chksum(fptr,status);
  }
  fits_close_file(fptr,status);
}


void sixt_check_obsolete_keyword(int* status){

	int num_obs_keys = 4;

	char* old_names[] = {"EventList","EventFile","PatternList","PatternFile","TesEventFile","eroEventList"};
	char* new_names[] = {"RawData","RawData","EvtFile","EvtFile","EvtFile","eroEvtFile"};

	int retval;

	for (int ii=0; ii<num_obs_keys; ii++){

		char* sbuffer;
		retval=ape_trad_query_string(old_names[ii], &sbuffer);
		if (retval==EXIT_SUCCESS){
			strtoupper(sbuffer);
			if (0!=strcmp(sbuffer,"NONE")) {
				printf("*** error reading parameters:  old name '%s' is obsolete. It is now called '%s'!  \n",
						old_names[ii],new_names[ii]);
				*status=EXIT_FAILURE;
			}
		}

	}

}


/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
int binary_search(double val, double* arr, int n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

int binary_search_long(double val, double* arr, long n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	long high=n-1;
	long low=0;
	long mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

int binary_search_float(float val, float* arr, int n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

int binary_search_float_long(float val, float* arr, long n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	long high=n-1;
	long low=0;
	long mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

