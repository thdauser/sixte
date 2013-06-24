#include "ero_rawevents.h"


void add_fits_erostdkeywords(fitsfile* const fptr, 
			     const int hdunum, 
			     char* const creation_date, 
			     char* const date_obs, 
			     char* const time_obs,
			     char* const date_end, 
			     char* const time_end, 
			     double tstart, 
			     double tstop, 
			     double timezero, 
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
  char instrume[MAXMSG]="INSTRUME";
  fits_update_key(fptr, TSTRING, "INSTRUME", instrume, "", status);

  char obsmode[MAXMSG]="";
  fits_update_key(fptr, TSTRING, "OBSMODE", obsmode, "", status);
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
  fits_update_key(fptr, TDOUBLE, "DE_OBJ", &de_obj, "[deg] J2000", status);

  fits_update_key(fptr, TSTRING, "DATE", creation_date, "File creation date", status);
  fits_update_key(fptr, TSTRING, "DATE-OBS", date_obs, "UT date of observation start", status);
  fits_update_key(fptr, TSTRING, "TIME-OBS", time_obs, "UT time of observation start", status);
  fits_update_key(fptr, TSTRING, "DATE-END", date_end, "UT date of observation end", status);
  fits_update_key(fptr, TSTRING, "TIME-END", time_end, "UT time of observation end", status);

  fits_update_key(fptr, TDOUBLE, "TSTART", &tstart, "Start time of exposure in units of TIME column", status);
  fits_update_key(fptr, TDOUBLE, "TSTOP", &tstop, "Stop time of exposure in units of TIME column", status);
  fits_update_key(fptr, TDOUBLE, "TEND", &tstop, "End time of exposure in units of TIME column", status);

  double mjdref=54101.0;
  fits_update_key(fptr, TDOUBLE, "MJDREF", &mjdref, "[d] 2007-01-01T00:00:00", status);
  
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


int ero_rawevents_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input event file.
  EventListFile* elf=NULL;

  // File pointer to the output eROSITA event file. 
  fitsfile* fptr=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("ero_rawevents");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop.

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;


    // Open the input event file.
    elf=openEventListFile(par.EventList, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read keywords from the input file.
    char comment[MAXMSG];
    float timezero=0.0;
    fits_read_key(elf->fptr, TFLOAT, "TIMEZERO", &timezero, comment, &status);
    CHECK_STATUS_BREAK(status);

    char date_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-OBS", date_obs, comment, &status);
    CHECK_STATUS_BREAK(status);

    char time_obs[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-OBS", time_obs, comment, &status);
    CHECK_STATUS_BREAK(status);

    char date_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "DATE-END", date_end, comment, &status);
    CHECK_STATUS_BREAK(status);

    char time_end[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TIME-END", time_end, comment, &status);
    CHECK_STATUS_BREAK(status);

    double tstart=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTART", &tstart, comment, &status);
    CHECK_STATUS_BREAK(status);

    double tstop=0.0;
    fits_read_key(elf->fptr, TDOUBLE, "TSTOP", &tstop, comment, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the creation date in the header.
    char creation_date[MAXMSG];
    int timeref;
    fits_get_system_time(creation_date, &timeref, &status);
    CHECK_STATUS_BREAK(status);

    // Check if the output file already exists.
    int exists;
    fits_file_exists(par.eroEventList, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.eroEventList);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.eroEventList);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }
    
    // Create and open a new FITS file using the template.
    headas_chat(3, "create new eROSITA event list file '%s' ...\n",
		par.eroEventList);
    fits_create_file(&fptr, par.eroEventList, &status);
    CHECK_STATUS_BREAK(status);

    // Create the event table.
    char *ttype[]={"TIME", "FRAME", "RAWX", "RAWY", "PHA"};
    char *tunit[]={"", "", "", "", "adu"};
    char *tform[]={"D", "J", "I", "I", "I"}; 
    fits_create_tbl(fptr, BINARY_TBL, 0, 5, ttype, tform, tunit, 
		    "EVENTS", &status);
    if (EXIT_SUCCESS!=status) {
      char msg[MAXMSG];
      sprintf(msg, "could not create binary table for events "
	      "in file '%s'", par.eroEventList);
      SIXT_ERROR(msg);
      break;
    }

    // Insert header keywords.
    char hduclass[MAXMSG]="OGIP";
    fits_update_key(fptr, TSTRING, "HDUCLASS", hduclass, "", &status);
    char hduclas1[MAXMSG]="EVENTS";
    fits_update_key(fptr, TSTRING, "HDUCLAS1", hduclas1, "", &status);
    CHECK_STATUS_BREAK(status);

    // Insert the standard eROSITA header keywords.
    add_fits_erostdkeywords(fptr, 1, creation_date, date_obs, time_obs,
			    date_end, time_end, tstart, tstop, 
			    timezero, &status);
    CHECK_STATUS_BREAK(status);
    add_fits_erostdkeywords(fptr, 2, creation_date, date_obs, time_obs,
			    date_end, time_end, tstart, tstop, 
			    timezero, &status);
    CHECK_STATUS_BREAK(status);

    // Determine the column numbers.
    int ctime, cframe, crawx, crawy, cpha;
    fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &status);
    fits_get_colnum(fptr, CASEINSEN, "FRAME", &cframe, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWX", &crawx, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWY", &crawy, &status);
    fits_get_colnum(fptr, CASEINSEN, "PHA", &cpha, &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---

    
    // --- Beginning of Copy Process ---

    headas_chat(3, "start copy process ...\n");

    // Loop over all events in the input file. 
    long row;
    for (row=0; row<elf->nrows; row++) {
      
      // Read the next event from the input file.
      Event event;
      getEventFromFile(elf, row+1, &event, &status);
      CHECK_STATUS_BREAK(status);

      // Store it in the output file.
      fits_insert_rows(fptr, row, 1, &status);
      CHECK_STATUS_BREAK(status);

      fits_write_col(fptr, TDOUBLE, ctime, row+1, 1, 1, 
		     &event.time, &status);
      fits_write_col(fptr, TLONG, cframe, row+1, 1, 1, 
		     &event.frame, &status);
      fits_write_col(fptr, TLONG, cpha, row+1, 1, 1, 
		     &event.pi, &status);
      int rawx=event.rawx+1;
      fits_write_col(fptr, TINT, crawx, row+1, 1, 1, &rawx, &status);
      int rawy=event.rawy+1;
      fits_write_col(fptr, TINT, crawy, row+1, 1, 1, &rawy, &status);

      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all events in the input file.

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  freeEventListFile(&elf, &status);

  // Append the check sum to the FITS header.
  if (NULL!=fptr) {
    int mode;
    fits_file_mode(fptr, &mode, &status);
    if (READWRITE==mode) {
      int ii;
      for (ii=0; ii<2; ii++) {
	int hdutype;
	fits_movabs_hdu(fptr, ii+1, &hdutype, &status);
	fits_write_chksum(fptr, &status);
      }
    }
    fits_close_file(fptr, &status);
  }
  
  if (EXIT_SUCCESS==status) headas_chat(3, "finished successfully\n\n");
  return(status);
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the input file");
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("eroEventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the output file");
    return(status);
  } 
  strcpy(par->eroEventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  return(status);
}

