#include "phabkg.h"

#include <math.h>
#include <sys/timeb.h>


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
  phabkg->randgen     =NULL;

  // Initialize values.
  phabkg->nbins=0;

  // Initialize GSL random number generator.
  // TODO Replace this by another random number generator.
  struct timeb time_struct;
  ftime(&time_struct);
  srand((unsigned int)time_struct.millitm-(unsigned int)time_struct.time);
  phabkg->randgen=gsl_rng_alloc(gsl_rng_ranlux);
  gsl_rng_set(phabkg->randgen, rand()*((int)time_struct.time+(unsigned int)time_struct.millitm));


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
		  phabkg->channel, phabkg->channel, &anynul, status);
    fits_read_col(fptr, TFLOAT, crate, 1, 1, phabkg->nbins, 
		  phabkg->distribution, phabkg->distribution, &anynul, status);
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
    gsl_rng_free((*phabkg)->randgen);

    free(*phabkg);
    *phabkg=NULL;
  }
}


long* PHABkgGetEvents(const PHABkg* const phabkg, 
		      /** Regarded time interval in [s]. */
		      const double interval, 
		      /** Detector area in [m^2]. */
		      const float area,
		      unsigned int* const nevts, 
		      int* const status)
{
  long* phas=NULL;

  // Check if everything is set up properly.
  if (NULL==phabkg) {
    SIXT_ERROR("no PHA background model specified");
    *status=EXIT_FAILURE;
    return(phas);
  }

  if (NULL==phabkg->distribution) {
    SIXT_ERROR("no PHA background spectrum loaded");
    *status=EXIT_FAILURE;
    return(phas);
  }
    
  // Determine the on average expected event number.
  // Note that the rates are given in [counts/s/bin/cm^2], whereas
  // the illuminted detector area is given in [m^2]. Therefore we 
  // need a conversion factor of 1.e4.
  double mu=phabkg->distribution[phabkg->nbins-1]*area*1.e4*interval;

  // Determine the number of background events according to 
  // Poisson statistics.
  *nevts=gsl_ran_poisson(phabkg->randgen, mu);

  // Allocate memory.
  phas=(long*)malloc((*nevts)*sizeof(long));
  CHECK_NULL_RET(phas, *status, 
		 "memory allocation for PHA values of background events failed", 
		 phas);

  // Determine random PHA values according to the spectral distribution.
  unsigned int ii;
  for (ii=0; ii<*nevts; ii++) {
    double r=sixt_get_random_number(status)*phabkg->distribution[phabkg->nbins-1];
    CHECK_STATUS_RET(*status, phas);
    
    // Perform a binary search to obtain the corresponding detector channel.
    long min=0;
    long max=phabkg->nbins-1;
    long mid;
    while (max > min) {
      mid=(min+max)/2;
      if (phabkg->distribution[mid] < r) {
	min=mid+1;
      } else {
	max=mid;
      }
    }

    phas[ii]=phabkg->channel[min];
  }

  return(phas);
}
