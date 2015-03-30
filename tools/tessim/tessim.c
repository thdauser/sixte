#include <assert.h>
#include <parinput.h>
#include "tessim.h"
#include <stdlib.h>


#define TOOLSUB tessim_main
#include "headas_main.c"

void tessim_getpar(tespxlparams *par, int *properties, int *status);

int tessim_main() {
  // error status
  int status=EXIT_SUCCESS;
  
  // register HEATOOL
  set_toolname("tessim");
  set_toolversion("0.01");


  do { // start of error handling loop
    // read parameters using PIL
    tespxlparams par;
    int properties;


    tessim_getpar(&par,&properties,&status);
    CHECK_STATUS_BREAK(status);

    // initialize the TES pixel
    tesparams *tes=tes_init(&par,&status);
    CHECK_STATUS_BREAK(status);


    // output information about simulation
    tes_print_params(tes);

    // exit program if we only want to display the properties
    // of the TES
    if (properties) {
      status=EXIT_SUCCESS; // NEED TO CHECK HERE
      break;
    }

    // Run the simulation
    tes_propagate(tes,par.tstop,&status);
    CHECK_STATUS_BREAK(status);

    // if we're done: free all allocated memory and deallocate the TES
    // (not really necessary here, but let's be anally retentive)
    tes_free(tes);
    free(tes);

    free(par.ID);
    free(par.impactlist);
    free(par.streamfile);

  } while(0); // end of error handling loop

  if (EXIT_SUCCESS==status) {
    headas_chat(3,"finished successfully!\n\n");
    return(EXIT_SUCCESS);
  }
  return(EXIT_FAILURE);
}

void tessim_getpar(tespxlparams *par, int *properties, int *status) {
  query_simput_parameter_string("PixelID", &(par->ID), status);
  query_simput_parameter_file_name("PixImpList", &(par->impactlist), status);
  query_simput_parameter_file_name("Streamfile", &(par->streamfile), status);
  query_simput_parameter_double("tstart", &(par->tstart), status);
  query_simput_parameter_double("tstop", &(par->tstop), status);
  query_simput_parameter_double("sample_rate",&(par->sample_rate),status);
  assert(par->sample_rate>0);
  query_simput_parameter_double("T_start", &(par->T_start), status);
  assert(par->T_start>0);
  query_simput_parameter_double("Tb", &(par->Tb), status);
  assert(par->Tb>0);
  query_simput_parameter_double("R0", &(par->R0), status);
  assert(par->R0>0);
  query_simput_parameter_double("RL", &(par->RL), status);
  assert(par->RL>=0);
  query_simput_parameter_double("alpha", &(par->alpha), status);
  assert(par->alpha>0);
  query_simput_parameter_double("beta", &(par->beta), status);
  query_simput_parameter_double("Lin", &(par->Lin), status);
  assert(par->Lin>0);
  query_simput_parameter_double("n", &(par->n), status);
  query_simput_parameter_double("imin", &(par->imin), status);
  query_simput_parameter_double("imax", &(par->imax), status);

  if (par->imax<=par->imin) {
    SIXT_ERROR("imax MUST be larger imin\n");
    *status=EXIT_FAILURE;
    return;
  }

  query_simput_parameter_bool("simnoise", &par->simnoise, status);

  long seed; // needed because par->seed is an unsigned long
  query_simput_parameter_long("seed", &seed, status);
  par->seed=labs(seed); 

  query_simput_parameter_bool("propertiesonly",&properties,status);

  query_simput_parameter_bool("clobber", &par->clobber, status);
  
}
