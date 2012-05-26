#include "ero_events.h"


int ero_events_main() 
{
  // Containing all programm parameters read by PIL
  struct Parameters par; 

  // Input pattern file.
  PatternFile* plf=NULL;

  // File pointer to the output eROSITA event file. 
  fitsfile* fptr=NULL;

  // WCS data structure used for projection.
  struct wcsprm wcs = { .flag=-1 };
  // String buffer for FITS header.
  char* headerstr=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 


  // Register HEATOOL:
  set_toolname("ero_events");
  set_toolversion("0.05");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Check whether an appropriate WCS projection has been selected.
    if (strlen(par.Projection)!=3) {
      SIXT_ERROR("invalid WCS projection type");
      status=EXIT_FAILURE;
      break;
    }

    // Open the input pattern file.
    plf=openPatternFile(par.PatternList, READONLY, &status);
    if (EXIT_SUCCESS!=status) break;

    // Create and open the output eROSITA event file.
    // Remove old file, if it exists.
    remove(par.eroEventList);

    // Filename of the template file.
    char template[MAXMSG];
    strcpy(template, SIXT_DATA_PATH);
    strcat(template, "/templates/eroeventlist.tpl");
    
    // Create and open a new FITS file using the template.
    char buffer[MAXMSG];
    sprintf(buffer, "%s(%s)", par.eroEventList, template);
    headas_chat(4, "create new eROSITA event list file '%s' "
		"from template '%s' ...\n", 
		par.eroEventList, template);
    fits_create_file(&fptr, buffer, &status);
    if (EXIT_SUCCESS!=status) break;

    // Move to the binary table HDU.
    int hdutype;
    fits_movabs_hdu(fptr, 2, &hdutype, &status);
    if (EXIT_SUCCESS!=status) break;

    // Set the time-keyword in the Event List Header.
    // See also: Stevens, "Advanced Programming in the UNIX environment",
    // p. 155 ff.
    time_t current_time;
    if (0 != time(&current_time)) {
      struct tm* current_time_utc = gmtime(&current_time);
      if (NULL != current_time_utc) {
	char current_time_str[MAXMSG];
	if (strftime(current_time_str, MAXMSG, "%Y-%m-%dT%H:%M:%S", 
		     current_time_utc) > 0) {
	  // Return value should be == 19 !
	  fits_update_key(fptr, TSTRING, "DATE-OBS", current_time_str, 
			  "Start Time (UTC) of exposure", &status);
	  CHECK_STATUS_BREAK(status);
	}
      }
    } 
    // END of writing time information to Event File FITS header.

    // Determine the column numbers.
    int ctime, crawx, crawy, cframe, cpha, cenergy, cra, cdec, cx, cy, cccdnr;
    fits_get_colnum(fptr, CASEINSEN, "TIME", &ctime, &status);
    fits_get_colnum(fptr, CASEINSEN, "FRAME", &cframe, &status);
    fits_get_colnum(fptr, CASEINSEN, "PHA", &cpha, &status);
    fits_get_colnum(fptr, CASEINSEN, "ENERGY", &cenergy, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWX", &crawx, &status);
    fits_get_colnum(fptr, CASEINSEN, "RAWY", &crawy, &status);
    fits_get_colnum(fptr, CASEINSEN, "RA", &cra, &status);
    fits_get_colnum(fptr, CASEINSEN, "DEC", &cdec, &status);
    fits_get_colnum(fptr, CASEINSEN, "X", &cx, &status);
    fits_get_colnum(fptr, CASEINSEN, "Y", &cy, &status);
    fits_get_colnum(fptr, CASEINSEN, "CCDNR", &cccdnr, &status);
    CHECK_STATUS_BREAK(status);

    // Timing keywords.
    float frametime=50.0;
    fits_update_key(fptr, TFLOAT, "FRAMETIM", &frametime,
		    "[ms] nominal frame time", &status);
    float timezero=0.0;
    fits_update_key(fptr, TFLOAT, "TIMEZERO", &timezero,
		    "clock correction", &status);
    fits_update_key(fptr, TSTRING, "TIMEUNIT", "s",
		    "Time unit", &status);
    fits_update_key(fptr, TSTRING, "TIMESYS", "TT",
		    "Time system (Terrestial Time)", &status);
    CHECK_STATUS_BREAK(status);

    // Set up the WCS data structure.
    if (0!=wcsini(1, 2, &wcs)) {
      SIXT_ERROR("initalization of WCS data structure failed");
      status=EXIT_FAILURE;
      break;
    }
    wcs.naxis = 2;
    wcs.crpix[0] = 0.;
    wcs.crpix[1] = 0.;
    wcs.crval[0] = par.RefRA;
    wcs.crval[1] = par.RefDec;    
    wcs.cdelt[0] = 0.05/3600.;
    wcs.cdelt[1] = 0.05/3600.;
    strcpy(wcs.cunit[0], "deg");
    strcpy(wcs.cunit[1], "deg");
    strcpy(wcs.ctype[0], "RA---");
    strcat(wcs.ctype[0], par.Projection);
    strcpy(wcs.ctype[1], "DEC--");
    strcat(wcs.ctype[1], par.Projection);

    // Update the WCS keywords in the output event file.
    char keyword[MAXMSG];
    sprintf(keyword, "TCTYP%d", cx);
    fits_update_key(fptr, TSTRING, keyword, wcs.ctype[0], 
		    "projection type", &status);
    sprintf(keyword, "TCTYP%d", cy);
    fits_update_key(fptr, TSTRING, keyword, wcs.ctype[1], 
		    "projection type", &status);
    sprintf(keyword, "TCRVL%d", cx);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[0], 
		    "reference value", &status);
    sprintf(keyword, "TCRVL%d", cy);
    fits_update_key(fptr, TDOUBLE, keyword, &wcs.crval[1], 
		    "reference value", &status);
    fits_update_key(fptr, TSTRING, "RADECSYS", "FK5", "", &status);
    float equinox=2000.0;
    fits_update_key(fptr, TFLOAT, "EQUINOX", &equinox, "", &status);
    fits_update_key(fptr, TSTRING, "LONGSTR", "OGIP 1.0",
		    "to support multi-line COMMENT oder HISTORY records",
		    &status);
    CHECK_STATUS_BREAK(status);

    // --- END of Initialization ---

    
    // --- Beginning of Copy Process ---

    headas_chat(3, "start copy process ...\n");

    // Values for TLMIN and TLMAX header keywords in the output file.
    float tlmin_energy=0.;
    float tlmax_energy=0.;
    long tlmin_x=0;
    long tlmax_x=0;
    long tlmin_y=0;
    long tlmax_y=0;

    // Loop over all patterns in the FITS file. 
    long row;
    for (row=0; row<plf->nrows; row++) {
      
      // Read the next pattern from the input file.
      Pattern pattern;
      getPatternFromFile(plf, row+1, &pattern, &status);
      CHECK_STATUS_BREAK(status);

      // Store the pattern in the output file.
      fits_insert_rows(fptr, row, 1, &status);
      CHECK_STATUS_BREAK(status);

      fits_write_col(fptr, TDOUBLE, ctime, row+1, 1, 1, &pattern.time, &status);
      fits_write_col(fptr, TLONG, cframe, row+1, 1, 1, &pattern.frame, &status);
      fits_write_col(fptr, TLONG, cpha, row+1, 1, 1, &pattern.pha, &status);

      float energy = pattern.signal * 1000.; // [eV]
      fits_write_col(fptr, TFLOAT, cenergy, row+1, 1, 1, &energy, &status);
      if ((energy < tlmin_energy) || (0==row)) {
	tlmin_energy = energy;
      }
      if (energy > tlmax_energy) {
	tlmax_energy = energy;
      }

      int rawx = pattern.rawx+1;
      fits_write_col(fptr, TINT, crawx, row+1, 1, 1, &rawx, &status);
      int rawy = pattern.rawy+1;
      fits_write_col(fptr, TINT, crawy, row+1, 1, 1, &rawy, &status);

      long ra = (long)(pattern.ra*180./M_PI/1.e-6);
      if (pattern.ra < 0.) ra--;
      fits_write_col(fptr, TLONG, cra, row+1, 1, 1, &ra, &status);

      long dec = (long)(pattern.dec*180./M_PI/1.e-6);
      if (pattern.dec < 0.) dec--;
      fits_write_col(fptr, TLONG, cdec, row+1, 1, 1, &dec, &status);
      CHECK_STATUS_BREAK(status);
      
      // Convert world coordinates to image coordinates X and Y.
      double world[2] = { pattern.ra*180./M_PI, pattern.dec*180./M_PI };
      double imgcrd[2], pixcrd[2];
      double phi, theta;
      wcss2p(&wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &status);
      CHECK_STATUS_BREAK(status);
      long x = (long)pixcrd[0]; 
      if (pixcrd[0] < 0.) x--;
      long y = (long)pixcrd[1]; 
      if (pixcrd[1] < 0.) y--;
      fits_write_col(fptr, TLONG, cx, row+1, 1, 1, &x, &status);
      fits_write_col(fptr, TLONG, cy, row+1, 1, 1, &y, &status);
      if ((x < tlmin_x) || (0==row)) {
	tlmin_x = x;
      }
      if (x > tlmax_x) {
	tlmax_x = x;
      }
      if ((y < tlmin_y) || (0==row)) {
	tlmin_y = y;
      }
      if (y > tlmax_y) {
	tlmax_y = y;
      }

      fits_write_col(fptr, TINT, cccdnr, row+1, 1, 1, &par.CCDNr, &status);

      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);
    // END of loop over all patterns in the FITS file.


    // Update TLMIN and TLMAX header keywords.
    sprintf(keyword, "TLMIN%d", cenergy);
    fits_update_key(fptr, TFLOAT, keyword, &tlmin_energy, "", &status);
    sprintf(keyword, "TLMAX%d", cenergy);
    fits_update_key(fptr, TFLOAT, keyword, &tlmax_energy, "", &status);

    sprintf(keyword, "TLMIN%d", cx);
    fits_update_key(fptr, TLONG, keyword, &tlmin_x, 
		    "", &status);
    sprintf(keyword, "TLMAX%d", cx);
    fits_update_key(fptr, TLONG, keyword, &tlmax_x, 
		    "", &status);

    sprintf(keyword, "TLMIN%d", cy);
    fits_update_key(fptr, TLONG, keyword, &tlmin_y, 
		    "", &status);
    sprintf(keyword, "TLMAX%d", cy);
    fits_update_key(fptr, TLONG, keyword, &tlmax_y, 
		    "", &status);

    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Close the files.
  destroyPatternFile(&plf, &status);
  if (NULL!=fptr) {
    // If the file was opened in READWRITE mode, calculate
    // the check sum an append it to the FITS header.
    int mode;
    fits_file_mode(fptr, &mode, &status);
    if (READWRITE==mode) {
      fits_write_chksum(fptr, &status);
    }
    fits_close_file(fptr, &status);
  }
  
  // Release memory.
  wcsfree(&wcs);
  if (NULL!=headerstr) free(headerstr);

  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the input pattern list!\n", status);
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("eroEventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the output event list!\n", status);
    return(status);
  } 
  strcpy(par->eroEventList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("CCDNr", &par->CCDNr);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the CCDNr parameter!\n", status);
    return(status);
  }

  status=ape_trad_query_string("Projection", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the projection type!\n", status);
    return(status);
  } 
  strcpy(par->Projection, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("RefRA", &par->RefRA);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading RefRA!\n", status);
    return(status);
  } 

  status=ape_trad_query_float("RefDec", &par->RefDec);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading RefDEC!\n", status);
    return(status);
  } 

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


