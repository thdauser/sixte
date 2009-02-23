#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_rand.h"
#include "create_rnd_source.h"

#define TOOLSUB test_distrndsources_main
#include "headas_main.c"


#define NSOURCES 1000000

// validation checks, wheter the sources are distributed isotropically on the unit sphere
// and checks the powerlaw distribution (log N - log S)
#define VALIDATION 1
#define VTHRESH 0.9
#define MAX_BINS 500000


////////////////////////////////////////////
int test_distrndsources_main()
{
  double dec, rasc;
  double xc, yc, zc;
  double intensity;
  int counter;
  const double powerlaw_index = -1.6;   // powerlaw index for logN-logS distribution


#ifdef VALIDATION
  int xcount = 0, ycount = 0, zcount = 0;
  int j, bin_index;
  int bins[MAX_BINS];

  for (j=0;j<MAX_BINS;j++) {
    bins[j] = 0;
  }
#endif

  // initialize random number generator
  HDmtInit(1);

  printf("# astronomical coordinates\t\t\tcartesian coordinates\t\t\tintensity\n# declination\tright ascension\tr=1.0\t\tx\t\ty\t\tz\t\tS\n");

  for (counter = 0; counter<NSOURCES; counter++) {
    // create source at random position on the unit sphere with intensity according to powerlaw distribution
    create_rnd_source_position(&rasc, &dec);
    intensity = get_rnd_intensity(powerlaw_index);

    // re-calculate cartesian coordiantes
    xc = cos(dec) * cos(rasc);
    yc = cos(dec) * sin(rasc);
    zc = sin(dec);

#ifdef VALIDATION
    if (fabs(xc)>VTHRESH) xcount++;
    if (fabs(yc)>VTHRESH) ycount++;
    if (fabs(zc)>VTHRESH) zcount++;

    bin_index=(int)(intensity/0.1);
    for(j=0; (j<bin_index) && (j<MAX_BINS); j++) {
      bins[j]++;
    }
#endif

#ifndef VALIDATION
    // print the astronomical and cartesian coordinates to standard output
    printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dec,rasc,1.0,xc,yc,zc,intensity);
#endif
  }

#ifdef VALIDATION
//  printf("%d\t%d\t%d\n", xcount,ycount,zcount);

  for (j=0;j<MAX_BINS;j++) {
    printf("%lf\t%d\n", j*0.1, bins[j]);
  }
#endif

  // release HEADAS random number generator
  HDmtFree();

  return(EXIT_SUCCESS);
}
