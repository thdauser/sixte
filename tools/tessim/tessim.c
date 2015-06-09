#include <assert.h>
#include <parinput.h>
#include "tessim.h"
#include <stdlib.h>


// define TES_DATASTREAM to use the datastream interface
// define TES_TESRECORD to use the (better) TESRECORD interface
//#define TES_DATASTREAM
#define TES_TESRECORD



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
      if (strcmp(par.streamfile,"none")!=0) {
	fitsfile *fptr;
	fits_create_file_clobber(&fptr,par.streamfile,par.clobber,&status); 
	fits_create_tbl(fptr,BINARY_TBL,0,0,NULL,NULL,NULL,"TESDATASTREAM",&status);
	printf("XXX%iXXX\n",status);
	tes_fits_write_params(fptr,tes,&status);
	fits_close_file(fptr,&status);
	CHECK_STATUS_BREAK(status);
      }
      status=EXIT_SUCCESS; // NEED TO CHECK HERE
      break;
    }

    // initialize photon provider
    tes->photoninfo=tes_init_impactlist(par.impactlist,&status);
    CHECK_STATUS_BREAK(status);
    tes->get_photon=&tes_photon_from_impactlist; // note: function pointer!

    // initialize stream writer
#ifdef TES_DATASTREAM
    tes->streaminfo=tes_init_datastream(par.tstart,par.tstop,tes,&status);
    CHECK_STATUS_BREAK(status);
    tes->write_to_stream=&tes_append_datastream; // note: function pointer!
    tes->write_photon=NULL; 
#endif
#ifdef TES_TESRECORD
    tes->streaminfo=tes_init_tesrecord(par.tstart,par.tstop,tes,
				       par.streamfile,par.impactlist,
				       par.clobber,
				       ((tes_impactfile_info *) (tes->photoninfo))->keywords,
				       &status);
    CHECK_STATUS_BREAK(status);
    tes->write_to_stream=&tes_append_tesrecord;
    tes->write_photon=&tes_append_photon_tesrecord;
#endif

    // Run the simulation
    printf("RUNNING!\n");
    tes_propagate(tes,par.tstop,&status);
    printf("DONE!\n");
    CHECK_STATUS_BREAK(status);

    // write the file and free all allocated memory
#ifdef TES_DATASTREAM
    tes_close_datastream(tes,par.streamfile,par.impactlist,(tes_datastream_info *) tes->streaminfo,
			((tes_impactfile_info *) (tes->photoninfo))->keywords,&status);
    tes_free_datastream((tes_datastream_info **) (tes->streaminfo), &status);
#endif
#ifdef TES_TESRECORD
    tes_close_tesrecord(tes,&status);
    tes_free_tesrecord((tes_record_info **) &tes->streaminfo, &status);
#endif

    tes_free_impactlist((tes_impactfile_info **) &tes->photoninfo, &status);
    tes_free(tes);
    free(tes);

    free(par.type);
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
  query_simput_parameter_string("PixType", &(par->type), status);

  int fromfile=0;

  
  query_simput_parameter_int("PixID",&(par->id),status);
  query_simput_parameter_file_name("PixImpList", &(par->impactlist), status);
  query_simput_parameter_file_name("Streamfile", &(par->streamfile), status);
  query_simput_parameter_double("tstart", &(par->tstart), status);
  query_simput_parameter_double("tstop", &(par->tstop), status);

  if (strncmp(par->type,"file:",5)==0) {
    char *file=strdup(par->type+5);
    tes_fits_read_params(file,par,status);
    CHECK_STATUS_VOID(*status);
    fromfile=1;
  }

  if (!fromfile){
    query_simput_parameter_double("sample_rate",&(par->sample_rate),status);
    assert(par->sample_rate>0);

    query_simput_parameter_double("Ce",&(par->Ce1),status);
    assert(par->Ce1>0);
    par->Ce1*=1e-12; // pJ/K -> J/K

    query_simput_parameter_double("Gb",&(par->Gb1),status);
    assert(par->Gb1>0);
    par->Gb1*=1e-12; // pW/K -> W/K

    query_simput_parameter_double("T_start", &(par->T_start), status); //mK
    assert(par->T_start>0);
    par->T_start*=1e-3; // mK -> K
    query_simput_parameter_double("Tb", &(par->Tb), status); // mK
    assert(par->Tb>0);
    par->Tb*=1e-3; // mK->K
    query_simput_parameter_double("R0", &(par->R0), status); // mOhm
    assert(par->R0>0);
    par->R0*=1e-3; // mOhm->Ohm

    query_simput_parameter_double("I0", &(par->I0), status); // muA
    assert(par->I0>0);
    par->I0*=1e-6; // muA->A

    query_simput_parameter_double("RL", &(par->RL), status); // mOhm
    assert(par->RL>=0);
    par->RL*=1e-3; // mOhm->Ohm
    query_simput_parameter_double("alpha", &(par->alpha), status);
    assert(par->alpha>0);
    query_simput_parameter_double("beta", &(par->beta), status);
    query_simput_parameter_double("Lin", &(par->Lin), status); //nH
    assert(par->Lin>0);
    par->Lin*=1e-9;
    query_simput_parameter_double("n", &(par->n), status);
    query_simput_parameter_double("imin", &(par->imin), status);
    query_simput_parameter_double("imax", &(par->imax), status);

    if (par->imax<=par->imin) {
      SIXT_ERROR("imax MUST be larger imin\n");
      *status=EXIT_FAILURE;
      return;
    }

    if (&par->simnoise) {
      query_simput_parameter_double("m_excess",&par->m_excess,status);
    } else {
      par->m_excess=0.; // switch explicitly off if not simulating noise
    }
    assert(par->m_excess>=0);
  }

  // read for all input parameter possibilities
  query_simput_parameter_bool("simnoise", &par->simnoise, status);

  long seed; // needed because par->seed is an unsigned long
  query_simput_parameter_long("seed", &seed, status);
  par->seed=labs(seed); 

  query_simput_parameter_bool("propertiesonly",properties,status);

  query_simput_parameter_bool("clobber", &par->clobber, status);
  
}
