#ifndef _ERODETBKGRNDGEN_H_
#define _ERODETBKGRNDGEN_H_ 1

#include <math.h>
#include <sys/timeb.h>

#include "sixt.h"
#include "simput.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "fitsio.h"

/** needed for seeding the random number generator with
 * a value that changes sufficiently fast */
struct timeb time_struct;

/** structure which contains all the input data from
 * the background simulation file */
typedef struct eroBackgroundInput {
  fitsfile *inputfptr;
  long numrows;
  double interval;
  double intervalsum;
  gsl_rng *randgen;
  
  char *timecolname;
  char *energycolname;
  char *xcolname;
  char *ycolname;

  int timecolnum;
  int energycolnum;
  int xcolnum;
  int ycolnum;
  int numevents;
  
  double eventsperinterval;
  double *hit_xpos;
  double *hit_ypos;
  double *hit_time;
  double *hit_energy;

  size_t *eventlist;
} eroBackgroundInput;

/** structure which contains the output data and
 * which will be passed to the calling routine */
typedef struct eroBackgroundOutput {
  int numevents;
  int numhits;
  double *hit_xpos;		// unit as specified in input file
  double *hit_ypos;
  double *hit_time;
  double *hit_energy;
} eroBackgroundOutput;

/** information extracted from a lightcurve which allows to modify
 * the background rate with an arbitrary function */
typedef struct eroBackgroundRateFct {
  long numelements;

  double* time;
  float* rate;

  double starttime;
  double* currenttime;
  double intervalsum;
  float* currentrate;
  float currentslope;
} eroBackgroundRateFct;

struct rateCurrentInterval {
  long numelements;
  double interval;

  float* rate;
  long ratesize;
};

void eroBkgSetRateFct(const char* const filename, int* const status);

/** return a randomly chosen eventlist for the given interval;
 * the number of events in the list is determined by poisson statistics
 * while the choice of the events themselves is based on a flat random
 * distribution. This function requires eroBkgInitialize to be called
 * once before being ready to deliver data. As soon as no more data will
 * be needed one should call the function eroBkgCleanUp in order to
 * release memory and close the input files. */
eroBackgroundOutput* eroBkgGetBackgroundList(const double interval);

/** open the simulation data file and initialize the random number
 * generator and the main structure */
void eroBkgInitialize(const char* const filename, int* const status);

/** free memory of passed eroBackgroundOutput structure */
void eroBkgFree(eroBackgroundOutput* struct_to_free);

/** try to cleanup everything and close the input file */
void eroBkgCleanUp(int* const status);

/** determine how many events occur in the whole background simulation */
int calcEvents(const double* const hit_time, const long numrows);

/** calculate event rate out of interval and number of events */
double calcEventRate(const double* hit_time,
                            const long numrows,
                            const int numevents,
                            const double interval);

void fillEventList(const double* const hit_time, const long numrows, size_t* const eventlist);

#endif   /* _ERODETBKGRNDGEN_H_ */
