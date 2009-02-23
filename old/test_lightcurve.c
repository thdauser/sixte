#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_randist.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define MAXMSG 256

#include "vector.h"
#include "photon.h"
#include "sources.h"
#include "fits_ctlg.h"
#include "random.h"
#include "global_constants.h"

#define TOOLSUB test_lightcurve_main
#include "headas_main.c"


int test_lightcurve_main()
{
  struct source_cat_entry src;
  src.lightcurve = NULL;
  const gsl_rng_type *T;   // type of GSL random number generator
  gsl_rng *r;              // pointer to GSL random number generator


  // initialize random number generator
  HDmtInit(1);
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);


  // call the lightcurve creation routine (same as used by measurement)
  create_lightcurve(&src, 0., r);


  // plot the lightcurve to stdout
  int count;
  for (count=0; count<N_LIGHTCURVE_BINS; count++) {
    printf("%lf %lf\n", src.lightcurve[count].t, src.lightcurve[count].rate);
  }       


  // free random number generators
  HDmtFree();
  gsl_rng_free(r);

  return (EXIT_SUCCESS);
}
