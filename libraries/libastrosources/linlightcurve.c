#include "linlightcurve.h"


LinLightCurve* getLinLightCurve(long nvalues, int* status)
{
  // Create an empty LinLightCurve object by allocating the required memory.
  LinLightCurve* lc=(LinLightCurve*)malloc(sizeof(LinLightCurve));
  if (NULL==lc) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for LinLightCurve!\n", 
		   EXIT_FAILURE);
    return(NULL);
  }
  lc->a = (double*)malloc(nvalues*sizeof(double));
  if (NULL==lc->a) {
    free(lc);
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for LinLightCurve!\n", 
		   EXIT_FAILURE);
    return(NULL);
  }
  lc->b = (double*)malloc(nvalues*sizeof(double));
  if (NULL==lc->b) {
    free(lc->a);
    free(lc);
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for LinLightCurve!\n", 
		   EXIT_FAILURE);
    return(NULL);
  }
  lc->u = (double*)malloc(nvalues*sizeof(double));
  if (NULL==lc->u) {
    free(lc->a);
    free(lc->b);
    free(lc);
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for LinLightCurve!\n", 
		   EXIT_FAILURE);
    return(NULL);
  }
  // Set the number of data points in the light curve.
  lc->nvalues = nvalues;  

  return(lc);
}



int initConstantLinLightCurve(LinLightCurve* lc, double mean_rate, double t0,
			      double step_width)
{
  if (NULL==lc) { 
    // The LinLightCurve object already must exist.
    return(EXIT_FAILURE);
  }
  if (1>lc->nvalues) { 
    // The LinLightCurve object must provide memory for at least 1 data point.
    return(EXIT_FAILURE);
  }
  // Set basic properties of the light curve.
  lc->t0 = t0;
  lc->step_width = step_width;

  // Set the data points of the light curve.
  long index;
  for (index=0; index<lc->nvalues; index++) {
    lc->a[index] = 0.;
    lc->b[index] = mean_rate;
    lc->u[index] = 1.-exp(-mean_rate*step_width);
  }

  return(EXIT_SUCCESS);
}



/*
LinLightCurve* getTimmerLinLightCurve(double mean_rate, double sigma, double t0, 
				long nbins, double bin_width, int* status)
{
  // Get an empty, uninitialized LinLightCurve object.
  LinLightCurve* lc = getLinLightCurve(nbins, status);
  if (EXIT_SUCCESS!=*status) return(NULL);

  // Set basic properties of the light curve.
  lc->t0 = t0;
  lc->bin_width = bin_width;

  // TODO Create light curve data according to the algorithm proposed by Timmer & Koenig.


  return(lc);
}
*/



void freeLinLightCurve(LinLightCurve* lightcurve)
{
  if (NULL!=lightcurve) {
    if (NULL!=lightcurve->a) {
      free(lightcurve->a);
      lightcurve->a=NULL;
    }
    if (NULL!=lightcurve->b) {
      free(lightcurve->b);
      lightcurve->b=NULL;
    }
    if (NULL!=lightcurve->u) {
      free(lightcurve->u);
      lightcurve->u=NULL;
    }
    lightcurve->nvalues=0;
  }
}



double getPhotonTime(LinLightCurve* lc, double time)
{
  if (1>lc->nvalues) {
    // Invalid LinLightCurve!
    return(-1.);
  }

  // The LinLightCurve contains data points, so the
  // general algorithm proposed by Klein & Roberts has to 
  // be applied.

  // Step 1 in the algorithm.
  double u = get_random_number();
  
  // Determine the respective index k of the light curve.
  long k = (long)((time-lc->t0)/lc->step_width);
  while (k < lc->nvalues) {
    // Step 2 in the algorithm is not required, as the u_k are already stored as
    // auxiliary data in the LinLightCurve object.
    // Step 3 in the algorithm.
    if (u <= lc->u[k]) {
      if (0. != lc->a[k]) {
	return((-lc->b[k] + sqrt(pow(lc->b[k],2.) + pow(lc->a[k]*time,2.) + 
				 2.*lc->a[k]*lc->b[k]*time - 2.*lc->a[k]*log(1.-u)))
	       /lc->a[k]);
      } else { // a_k == 0
	return(time-log(1.-u)/lc->b[k]);
      }

    } else {
      // Step 4 (u > u_k).
      u = (u-lc->u[k])/(1-lc->u[k]);
      k++;
      time = lc->t0 + k*lc->step_width;
    }
  }

  // The range of the light curve was exceeded.
  return(-1.);
}


