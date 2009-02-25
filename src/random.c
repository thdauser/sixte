#include <math.h>
#include <limits.h>

#include "random.h"
#include "headas_rand.h"


// This routines Creates random numbers using the HEAdas random number generator.
double get_random_number()
{
  // return value out of [0,1)
  return((double)HDmtRand()/ULONG_MAX);
}


/*
// creates a random number in the interval [0,1] (!!!!!!!!) and should therefore not be used !!!
double get_random_number()
{
  const int a = 1664525;
  const int b = 1013904223;
  static int x = 678549;   // start-value

  // generate next integer number
  x = a * x + b;

  // return value out of [0,1)
  return(fabs((double)x/INT_MIN));
}
*/



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns a random value on the basis of an exponential distribution with a given average distance.
// Here this function is used to calculate the temporal differences between individual photons from a source.
double rndexp(double avgdist)
{
  double rand = get_random_number();
  if (rand < 1.E-15) {
    rand = 1.E-15;
  }
  
  return(-log(rand)*avgdist);
}



