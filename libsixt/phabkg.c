#include "phabkg.h"

#include <math.h>


////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////


PHABkg* newPHABkg(const char* const filename, int* const status) 
{
  // Allocate memory.
  PHABkg* phabkg=(PHABkg*)malloc(sizeof(PHABkg));
  if (NULL==phabkg) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for PHABkg failed");
    return(phabkg);
  }

  // Initialize all pointers with NULL.
  phabkg->channel     =NULL;
  phabkg->distribution=NULL;

  // Initialize values.
  phabkg->nbins=0;

  // Load the specified PHA file.
  fitsfile *fptr=NULL;

  do { // Beginning of ERROR handling loop.

    headas_chat(5, "load PHA background spectrum from file '%s' ...\n", filename);

    // Open the file:
    fits_open_table(&fptr, filename, READONLY, status);
    CHECK_STATUS_BREAK(*status);

    // Get the HDU type.
    int hdutype;
    fits_get_hdu_type(fptr, &hdutype, status);
    CHECK_STATUS_BREAK(*status);

    // Image HDU results in an error message.
    if (IMAGE_HDU==hdutype) {
      *status=EXIT_FAILURE;
      char msg[MAXMSG];
      sprintf(msg, "no table extension available in FITS file '%s'", 
	      filename);
      SIXT_ERROR(msg);
      break;
    }

    // Determine the number of rows in the file.
    fits_get_num_rows(fptr, &phabkg->nbins, status);
    CHECK_STATUS_BREAK(*status);

    // Allocate memory.
    phabkg->channel=(long*)malloc(phabkg->nbins*sizeof(long));
    CHECK_NULL_BREAK(phabkg->channel, *status, 
		     "memory allocation for PHA background spectrum failed");
    phabkg->distribution=(float*)malloc(phabkg->nbins*sizeof(float));
    CHECK_NULL_BREAK(phabkg->distribution, *status, 
		     "memory allocation for PHA background spectrum failed");

    // Determine the column numbers.
    int cchannel, crate;
    fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &cchannel, status);
    fits_get_colnum(fptr, CASEINSEN,    "RATE",    &crate, status);
    CHECK_STATUS_BREAK(*status);

    // Read the data.
    int anynul=0;
    fits_read_col(fptr, TLONG, cchannel, 1, 1, phabkg->nbins, 
		  0, phabkg->channel, &anynul, status);
    fits_read_col(fptr, TFLOAT, crate, 1, 1, phabkg->nbins, 
		  0, phabkg->distribution, &anynul, status);
    CHECK_STATUS_BREAK(*status);

    // Sum up the rate values in order to obtain an accumulative distribution.
    long ii;
    for (ii=1; ii<phabkg->nbins; ii++) {
      phabkg->distribution[ii]+=phabkg->distribution[ii-1];
    }

  } while (0); // END of error handling loop  

  // Clean up:
  if (fptr) fits_close_file(fptr, status);

  return(phabkg);
}


void destroyPHABkg(PHABkg** const phabkg)
{
  if (NULL!=*phabkg) {
    if (NULL!=(*phabkg)->channel) {
      free((*phabkg)->channel);
    }
    if (NULL!=(*phabkg)->distribution) {
      free((*phabkg)->distribution);
    }
    free(*phabkg);
    *phabkg=NULL;
  }
}


int getPHABkgEvent(const PHABkg* const phabkg,
		   const float scaling,
		   const double tstart,
		   const double tstop,
		   double* const t,
		   long* const pha,
		   int* const status)
{
  static double tnext=0.0;

  // Check if everything is set up properly.
  if (NULL==phabkg) {
    SIXT_ERROR("no PHA background model specified");
    *status=EXIT_FAILURE;
    return(0);
  }

  if (NULL==phabkg->distribution) {
    SIXT_ERROR("no PHA background spectrum loaded");
    *status=EXIT_FAILURE;
    return(0);
  }

  // Determine the average event rate.
  double rate=phabkg->distribution[phabkg->nbins-1]*scaling;

  // Update to the start time, if time of next background event lies
  // before that.
  if (tnext<tstart) {
    tnext=tstart+rndexp(1./rate, status);
    CHECK_STATUS_RET(*status, 0);
  }

  // Check if the time of the next background event lies before the
  // specified upper limit. In this case no event is produced.
  if (tnext>tstop) {
    return(0);
  }

  // Set the time of the background event.
  *t=tnext;

  // Determine the PHA value.
  double r=sixt_get_random_number(status)*phabkg->distribution[phabkg->nbins-1];
  CHECK_STATUS_RET(*status, 0);

  // Perform a binary search to obtain the corresponding channel.
  long min=0;
  long max=phabkg->nbins-1;
  long mid;
  while (max>min) {
    mid=(min+max)/2;
    if (phabkg->distribution[mid]<r) {
      min=mid+1;
    } else {
      max=mid;
    }
  }
  *pha=phabkg->channel[min];

  // Determine the time of the next background event.
  tnext+=rndexp(1./rate, status);
  CHECK_STATUS_RET(*status, 0);

  return(1);
}
