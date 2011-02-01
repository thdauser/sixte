#ifndef _ERODETBKGRNDGEN_H_
#define _ERODETBKGRNDGEN_H_ 1

#include <math.h>
#include <sys/timeb.h>

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "fitsio.h"

enum Units {
  SECONDS = 1,
  MILLISECONDS = 1000,
  MICROSECONDS = 1000000
};

/** unit of interval */
static int INTERVAL_UNIT = MILLISECONDS;

/** needed for seeding the random number generator with
 * a value that changes sufficiently fast */
struct timeb time_struct;

/** structure which contains all the input data from
 * the background simulation file */
typedef struct eroBackgroundInput {
  fitsfile *inputfptr;
  long numrows;
  int interval;			// unit as specified above with INTERVAL_UNIT
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


/** return a randomly chosen eventlist for the given interval;
 * the number of events in the list is determined by poisson statistics
 * while the choice of the events themselves is based on a flat random
 * distribution. This function requires eroBkgInitialize to be called
 * once before being ready to deliver data. As soon as no more data will
 * be needed one should call the function eroBkgCleanup in order to
 * release memory and close the input files. */
eroBackgroundOutput* eroBkgGetBackgroundList(int interval);

/** open the simulation data file and initialize the random number
 * generator and the main structure */
void eroBkgInitialize(const char *filename, int *status);

/** free memory of passed eroBackgroundOutput structure */
void eroBkgFree(eroBackgroundOutput *struct_to_free);

/** try to cleanup everything and close the input file */
void eroBkgCleanUp(int *status);

/** determine how many events occur in the whole background simulation */
int calcEvents(double *hit_time, long numrows);

/** calculate event rate out of interval and number of events */
double calcEventRate(double *hit_time,
                     long numrows,
                     int numevents,
                     int interval);

#endif   /* _ERODETBKGRNDGEN_H_ */
