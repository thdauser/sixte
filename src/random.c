#include <math.h>
#include <limits.h>

#include "random.h"
#include "headas_rand.h"



//////////////////////////////////////////////
inline double get_random_number()
{
  // Return a value out of the interval [0,1):
  return((double)HDmtRand()/ULONG_MAX);
}



/////////////////////////////////////////////////////////////////////////
inline double rndexp(double avgdist)
{
  double rand = get_random_number();
  if (rand < 1.E-15) {
    rand = 1.E-15;
  }
  
  return(-log(rand)*avgdist);
}



