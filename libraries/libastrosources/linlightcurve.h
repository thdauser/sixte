#ifndef LINLIGHTCURVE_H
#define LINLIGHTCURVE_H 1

#include "sixt.h"
#include "sixt_random.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Piece-wise linear light curve giving the average source photon rate 
 * for a particular X-ray source. */
typedef struct {
  /** Number of data points in the light curve. */
  long nvalues;

  /** Width of the interval between the individual data points ([s]). */
  double step_width;

  /** Time of the first light curve data point ([s]). */
  double t0;

  /** Piece-wise linear light curve data points ([photons/s]). 
   * The value a_k represents the gradient of the light curve between 
   * the time t_k (= t0 + k*step_width) and t_{k+1}. 
   * The value b_k represents the contants contribution at t_k. */
  double* a, *b;

  /** Auxiliary data.
   * The values u_k are required in the algorithm proposed by Klein & Roberts
   * and are calculated in advance in order to accelerate the calculation
   * of the photon creation times. */
  double* u;

} LinLightCurve;


//////////////////////////////////////////////////////////////////////////
//   Function declarations
//////////////////////////////////////////////////////////////////////////


/** Constructor allocating memory for a LinLightCurve object.
 * The routine allocates memory required for a LinLightCurve object with
 * given number of data points. The values of the individual data points are not 
 * initialized, i.e., may have arbitrary values. */
LinLightCurve* getLinLightCurve(long nvalues, int* status);

/** Initialization routine creating a light curve with a constant photon rate. 
 * The light curve object already has to be allocated in advance by the
 * constructor getLinLightCurve().
 * The light curve starts at t0 ([s]) and has the specified length ([s]). 
 * The return value of the function is either EXIT_SUCCESS or EXIT_FAILURE. */
int initConstantLinLightCurve(LinLightCurve* lc, double mean_rate, double t0, 
			      double step_width);

/** Constructor generating a light curve according to the algorithm proposed
 * by Timmer & Koenig (1995). 
 * The function first allocates the memory required for the LinLightCurve
 * object. Then it generates the light curve data using the Timmer & Koenig
 * randomization process with the given parameters. */
//LinLightCurve* getTimmerLinLightCurve(double mean_rate, double sigma, double t0, 
//				long nvalues, double step_width, int* status);

/** Destructor. */
void freeLinLightCurve(LinLightCurve*);

/** Determine next photon generation time from current point of time and given
 * light curve. 
 * The time is obtained from a time-varying Poisson arrival process generator
 * proposed by Klein & Roberts (1984). This paper describes a Poisson 
 * generator for piecewise-linear light curves. */
double getPhotonTime(LinLightCurve*, double time);


#endif /* LINLIGHTCURVE_H */
