#include "gti.h"


GTI* newGTI(int* const status)
{
  GTI* gti=(GTI*)malloc(sizeof(GTI));
  CHECK_NULL_RET(gti, *status, 
		 "memory allocation for GTI data structure failed", gti);

  // Initialize pointers with NULL.
  gti->start=NULL;
  gti->stop =NULL;

  // Initialize values.
  gti->nentries=0;

  return(gti);
}


void freeGTI(GTI** const gti)
{
  if (NULL!=*gti) {
    if (NULL!=(*gti)->start) {
      free((*gti)->start);
    }
    if (NULL!=(*gti)->stop) {
      free((*gti)->stop);
    }
    free(*gti);
    *gti=NULL;
  }
}


GTI* loadGTI(const char* const filename, int* const status)
{
  GTI* gti=NULL;
  fitsfile* fptr=NULL;

  do { // Beginning of Error handling loop.

    // Open the FITS file for reading:
    fits_open_file(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);

    // Move to the second HDU (GTI table).
    int hdutype;
    fits_movabs_hdu(fptr, 2, &hdutype, status);
    CHECK_STATUS_BREAK(*status);
    
    // Image HDU results in an error message.
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "no table extension available in file '%s'", 
	      filename);
      SIXT_ERROR(msg);
      break;
    }

    // Determine the number of rows in the table.
    long nrows;
    fits_get_num_rows(fptr, &nrows, status);
    CHECK_STATUS_BREAK(*status);

    // Determine the individual column numbers.
    int cstart, cstop;
    fits_get_colnum(fptr, CASEINSEN, "START", &cstart, status);
    fits_get_colnum(fptr, CASEINSEN,  "STOP", &cstop , status);
    CHECK_STATUS_BREAK(*status);

    // Allocate memory.
    gti=newGTI(status);
    CHECK_NULL_BREAK(gti, *status, 
		     "memory allocation for GTI data structure failed");
    gti->start=(double*)malloc(nrows*sizeof(double));
    CHECK_NULL_BREAK(gti->start, *status, 
		     "memory allocation for GTI data structure failed");
    gti->stop =(double*)malloc(nrows*sizeof(double));
    CHECK_NULL_BREAK(gti->stop , *status, 
		     "memory allocation for GTI data structure failed");
    
    // Read the data from the table.
    int anynul = 0;
    fits_read_col(fptr, TDOUBLE, cstart, 1, 1, nrows, 
		  NULL, gti->start, &anynul, status);
    fits_read_col(fptr, TDOUBLE, cstop , 1, 1, nrows, 
		  NULL, gti->stop , &anynul, status);
    CHECK_STATUS_BREAK(*status);

    gti->nentries=nrows;

  } while(0); // END of error handling loop.

  if (NULL!=fptr) fits_close_file(fptr, status);

  return(gti);
}


void saveGTI(GTI* const gti,
	     const char* const filename,
	     const char clobber,
	     int* const status)
{
  fitsfile* fptr=NULL;

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_VOID(*status);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(filename);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", filename);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      return;
    }
  }

  // Create a new FITS file.
  fits_create_file(&fptr, filename, status);
  CHECK_STATUS_VOID(*status);

  // Set the time-keyword in the event list header.
  char datestr[MAXMSG];
  int timeref;
  fits_get_system_time(datestr, &timeref, status);
  CHECK_STATUS_VOID(*status);
  fits_update_key(fptr, TSTRING, "DATE", datestr, 
		  "File creation date", status);
  CHECK_STATUS_VOID(*status);

  // Store the GTI extension.
  saveGTIExt(fptr, "GTI", gti, status);
  CHECK_STATUS_VOID(*status);
  
  // Close the FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
}


void saveGTIExt(fitsfile* const fptr,
		char* const extname,
		GTI* const gti,
		int* const status)
{
  // Create the GTI table.
  char *ttype[]={"START", "STOP"};
  char *tform[]={"D", "D"};
  char *tunit[]={"", ""};
  fits_create_tbl(fptr, BINARY_TBL, 0, 2, ttype, tform, tunit, 
		  extname, status);
  if (EXIT_SUCCESS!=*status) {
    SIXT_ERROR("could not create binary table for GTI extension");
    return;
  }

  // Insert header keywords.
  fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", "", status);
  fits_update_key(fptr, TSTRING, "HDUCLAS1", "GTI", "", status);
  fits_update_key(fptr, TSTRING, "HDUCLAS2", "STANDARD", "", status);
  CHECK_STATUS_VOID(*status);

  // Determine the individual column numbers.
  int cstart, cstop;
  fits_get_colnum(fptr, CASEINSEN, "START", &cstart, status);
  fits_get_colnum(fptr, CASEINSEN,  "STOP", &cstop , status);
  CHECK_STATUS_VOID(*status);
  
  // Store the data in the table.
  fits_write_col(fptr, TDOUBLE, cstart, 
		 1, 1, gti->nentries, gti->start, status);
  fits_write_col(fptr, TDOUBLE, cstop, 
		 1, 1, gti->nentries, gti->stop, status);
  CHECK_STATUS_VOID(*status);
}


void appendGTI(GTI* const gti, 
	       const double start, 
	       const double stop, 
	       int* const status)
{
  // Allocate memory.
  gti->start=realloc(gti->start, (gti->nentries+1)*sizeof(double));
  CHECK_NULL_VOID(gti, *status, 
		  "memory allocation for GTI data structure failed");
  gti->stop =realloc(gti->stop , (gti->nentries+1)*sizeof(double));
  CHECK_NULL_VOID(gti, *status, 
		  "memory allocation for GTI data structure failed");

  // Store the start and the stop time.
  gti->start[gti->nentries]=start;
  gti->stop[gti->nentries] =stop;
  gti->nentries++;
}

