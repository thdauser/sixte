#include "tesgenimpacts.h"

#define TOOLSUB tesgenimpacts_main
#include "headas_main.c"

// TODO: Proper input error handling (are times, energies etc consistent)
/////////////////////////////////////

/** main procedure */
int tesgenimpacts_main() {

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("tesgenimpacts");
  set_toolversion("0.01");

  do { // Beginning of the ERROR handling loop (will at
       // most be run once).

    // Parameter input
    tesgenimppars par;
    tesgenimpacts_getpar(&par, &status);
    CHECK_STATUS_BREAK(status);

    // open the piximpfile
    PixImpFile* pixilf=openNewPixImpFile(par.pixFilename, "VIRTUAL", "VIRTUAL",
                                         "NONE", "NONE", "NONE", "NONE", "NONE",
                                         55000., 0., par.tstart, par.tstop, 'y', &status);
    // arrays for impact energies and times
    double Earray[par.nPhotons];
    double tarray[par.nPhotons];

    if (par.mode == MODE_CONST) {
    // mode 1: Generate photons every dtau. Stop when you're over tstop
      // populate the Earray and tarray
      int n;
      for (n=0; n<par.nPhotons; n++) {
        Earray[n] = par.EConst;
        tarray[n] = par.tstart + (n+1)*par.dtau; 
      }

    
    } else if (par.mode == MODE_LIN) {
    // mode 2: Generate nsamples photons
    // with linearly increasing energies
      // populate the Earray and tarray
      int n;
      for (n=0; n<par.nPhotons; n++){
        Earray[n] = par.Emin + (par.Emax - par.Emin) * (n)/(par.nPhotons-1);
        tarray[n] = par.tstart + (n+1)*par.dtau; 
      }


    } else if (par.mode == MODE_RAND || par.mode == MODE_SINRAND) {
      // rng initialization
      gsl_rng *rng;
      rng = gsl_rng_alloc(gsl_rng_taus);

      // if seed is 0, do NOT use the rng's default seed (per gsl), but
      // initialize from the system clock
      if (par.seed==0) {
#if defined( __APPLE__) && defined(__MACH__)
        struct timeval tv;
        gettimeofday(&tv,NULL);
        seed=1000000*tv.ts_sec+tv.tv_usec;
#else
        struct timespec tp;
        clock_gettime(CLOCK_MONOTONIC,&tp);
        par.seed=tp.tv_nsec;
#endif
      }
      gsl_rng_set(rng,par.seed);

      int n;
      if (par.mode == MODE_RAND) { 
        // probability distribution in t is linear
        for (n = 0; n <par.nPhotons; n++) {
        // populate the Earray and tarray
        // since we're using gsl_rng_uniform, we're technically never getting
        // Emax and tstop exactly.
          Earray[n] = gsl_rng_uniform(rng) *(par.Emax - par.Emin) + par.Emin;
          tarray[n] = gsl_rng_uniform(rng) * (par.tstop - par.tstart) + par.tstart;
        }

      } else if (par.mode == MODE_SINRAND) {
        // probability distribution in t is sinusoidal
        // normalization, such that the distribution times mu, integrated from tstart to tstop = 1
        double mu = 1 / (sin_int(par.tstop, par.offset, par.amplitude, par.f, par.shift) - sin_int(par.tstart, par.offset, par.amplitude, par.f, par.shift));

        // acceptance / rejection method for tarray
        for (n = 0; n <par.nPhotons; n++) {
          // populate the Earray and tarray
          // since we're using gsl_rng_uniform, we're technically never getting
          // Emax and tstop exactly.
          Earray[n] = gsl_rng_uniform(rng) *(par.Emax - par.Emin) + par.Emin;
          tarray[n] = sin_get_impact_time(rng, mu*par.offset, mu*par.amplitude, par.f, par.tstart, par.tstop, par.shift);
        }
      }
      // sort the arrival times ascending
      qsort(tarray, par.nPhotons, sizeof(double), cmp);
      // since energies were selected at random, we don't need to sort them along
      // free the rng
      gsl_rng_free(rng);
    }

    // add impacts to impactlist
    int n;
    for (n=0; n<par.nPhotons; n++) {
      // construct the impact
      PixImpact impact;
      impact.pixID = par.pixid;
      impact.time = tarray[n];
      impact.energy = Earray[n];
      struct Point2d fakepoint = {0., 0.};
      impact.detposition = fakepoint;
      impact.pixposition = fakepoint;
      impact.ph_id = n;
      impact.src_id = 1;
      // write it into the file
      addImpact2PixImpFile(pixilf, &impact, &status);
    }
    // File has been filled, close it
    freePixImpFile(&pixilf, &status);
  } while(0); // END of the error handling loop

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }

}

// comparision function, for qsort
int cmp(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return 1;
  return 0;
}

// sine probability distribution. Not normalized!
double sin_dist(double t, double offset, double amplitude, double f, double shift)
{
  return offset + amplitude * (sin(2*M_PI*f*t + shift));
}

// indefinite integral of sin_dist
double sin_int(double t, double offset, double amplitude, double f, double shift)
{
  return offset * t + amplitude * (cos(2*M_PI*f*t + shift) / (2*M_PI*f));
}

// get a photon impact time from the sine distribution using
// an acceptance-rejection method
double sin_get_impact_time(gsl_rng *rng, double offset, double amplitude, double f, double tstart, double tstop, double shift)
{
  // Note: the input distribution should already be normalized
  // i.e. multiply offset and amplitude by mu
  double T = tstop - tstart;
  // g is a uniform density from tstart to tstop
  // so its value is 1/T
  double g = 1/T;
  // M must be chosen such that M*g(t) >= sin_dist(t) for all t
  double M = T*(offset + amplitude);
  // y and z are uniform random numbers from [tstart, tstop) and [0,1)
  // respectively
  double z = 10.; // value is much too high, just for beginning loop
  double y = tstart;
  while ( z >= sin_dist(y, offset, amplitude, f, shift) / (M*g) ) {
    z = gsl_rng_uniform(rng);
    y = tstart + gsl_rng_uniform(rng)*(tstop - tstart);
  }
  return y;
}

// parameter input function
void tesgenimpacts_getpar(tesgenimppars *par, int *status){

  query_simput_parameter_string("PixImpList",&(par->pixFilename), status);
  query_simput_parameter_int("PixID", &(par->pixid), status);

  char *modestr;
  query_simput_parameter_string("mode", &(modestr), status);
  if (strcmp(modestr, "const") == 0){
    par->mode = MODE_CONST;
  } else if (strcmp(modestr, "lin") == 0){
    par->mode = MODE_LIN;
  } else if (strcmp(modestr, "rand") == 0){
    par->mode = MODE_RAND;
  } else if (strcmp(modestr, "sin") == 0){
    par->mode = MODE_SINRAND;
  } else {
    printf("ERROR: Mode must be either const, lin, rand or sin!\n");
    exit(1);
  }

  query_simput_parameter_int("PixID", &(par->pixid), status);
  par->pixid -= 1; // internal functions have a different offset;
  
  query_simput_parameter_double("tstart", &(par->tstart), status);
  query_simput_parameter_double("tstop", &(par->tstop), status);
  if (par->mode == MODE_CONST || par->mode == MODE_LIN){
    query_simput_parameter_double("dtau", &(par->dtau), status);
    par->dtau *= 1e-3; // ms -> s
    if (par->mode == MODE_CONST){
      // calculate nr of photons internally
      par->nPhotons = 0;
      double time = par->tstart + par->dtau;
      while (time < par->tstop){
        par->nPhotons += 1;
        time += par->dtau;
      }
    }
    
  } 
  if (par->mode == MODE_CONST) {
    query_simput_parameter_double("EConst", &(par->EConst), status);
  }
  if (par->mode == MODE_LIN || par->mode == MODE_RAND || par->mode == MODE_SINRAND){
    query_simput_parameter_double("EMin", &(par->Emin), status);
    query_simput_parameter_double("EMax", &(par->Emax), status);
    query_simput_parameter_int("PhotSamples", &(par->nPhotons), status);
    if (par->mode == MODE_RAND || par->mode == MODE_SINRAND){
      query_simput_parameter_long("Seed", &(par->seed), status);
      if (par->mode == MODE_SINRAND) {
        query_simput_parameter_double("sinOffset", &(par->offset), status);
        query_simput_parameter_double("sinAmplitude", &(par->amplitude), status);
        query_simput_parameter_double("sinFreq", &(par->f), status);
        query_simput_parameter_double("phaseShift", &(par->shift), status);
        par->shift = par->shift * M_PI / 180; // degrees->radiants
      }
    } 
  }
}

