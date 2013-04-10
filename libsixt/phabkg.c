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
  phabkg->fov_diameter=0.0;

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


unsigned int PHABkgGetEvents(const PHABkg* const phabkg, 
			     /** Regarded time interval in [s]. */
			     const double interval, 
			     const GenPixGrid* const pixgrid,
			     long** phas,
			     int** x,
			     int** y,
			     int* const status)
{
  // Check if everything is set up properly.
  if (NULL==pixgrid) {
    SIXT_ERROR("no pixel grid specified");
    *status=EXIT_FAILURE;
    return(0);
  }

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
    
  // Determine the on average expected event number.
  // Note that the rates are given in [counts/s/bin/cm^2], whereas
  // the illuminted detector area is given in [m^2]. Therefore we 
  // need a conversion factor of 1.e4.
  double mu=
    phabkg->distribution[phabkg->nbins-1]*
    pixgrid->xwidth*pixgrid->xdelt*
    pixgrid->ywidth*pixgrid->ydelt*
    1.e4*interval;

  // Determine the number of background events according to 
  // Poisson statistics.
  unsigned int nevts=gsl_ran_poisson(phabkg->randgen, mu);

  // Allocate memory. Note that too much memory might be allocated, because
  // some background events are rejected because of their location outside
  // the FoV.
  *phas=(long*)malloc(nevts*sizeof(long));
  CHECK_NULL_RET(*phas, *status, 
		 "memory allocation for PHA values of background events failed", 0);
  *x=(int*)malloc(nevts*sizeof(int));
  CHECK_NULL_RET(*x, *status, 
		 "memory allocation for pixel indices of background events failed", 0);
  *y=(int*)malloc(nevts*sizeof(int));
  CHECK_NULL_RET(*y, *status, 
		 "memory allocation for pixel indices of background events failed", 0);

  // Determine random PHA and pixel values according to the spectral distribution.
  unsigned int nacc=0;
  unsigned int ii;
  double cosrota=cos(pixgrid->rota);
  double sinrota=sin(pixgrid->rota);
  for (ii=0; ii<nevts; ii++) {

    // Determine the pixel indices.
    int xi=(int)(sixt_get_random_number(status)*pixgrid->xwidth);
    CHECK_STATUS_RET(*status, 0);
    int yi=(int)(sixt_get_random_number(status)*pixgrid->ywidth);
    CHECK_STATUS_RET(*status, 0);

    // Check if this pixel lies within the regarded region.
    if (phabkg->fov_diameter>0.0) {
      // If a FoV radius is defined, we have to check, whether the 
      // selected pixel is within the FoV.
      double xr=(xi-pixgrid->xrpix+1.0)*pixgrid->xdelt;
      double yr=(yi-pixgrid->yrpix+1.0)*pixgrid->ydelt;
      if (pow(xr*cosrota+yr*sinrota+pixgrid->xrval,2.0)+
	  pow(-xr*sinrota+yr*cosrota+pixgrid->yrval,2.0)
	  >
	  pow(phabkg->fov_diameter*0.5, 2.0)) {
	continue;
      }
    }

    (*x)[nacc]=xi;
    (*y)[nacc]=yi;

    // Determine the PHA value.
    double r=sixt_get_random_number(status)*phabkg->distribution[phabkg->nbins-1];
    CHECK_STATUS_RET(*status, 0);

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
    (*phas)[nacc]=phabkg->channel[min];

    nacc++;
  }

  return(nacc);
}
