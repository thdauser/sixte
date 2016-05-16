#include <assert.h>
#include <parinput.h>
#include "tessim.h"
#include <stdlib.h>


#define TOOLSUB tessim_main
#include "headas_main.c"

//
// THINGS TO FIX:
// * parfile seems not to update properly (email megan, 22 Aug)
//

void tessim_getpar(tespxlparams *par, int *properties, int *status);

int tessim_main() {
  // error status
  int status=EXIT_SUCCESS;
  
  // register HEATOOL
  set_toolname("tessim");
  set_toolversion("0.02");

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

    SixtStdKeywords *keywords=duplicateSixtStdKeywords(((tes_impactfile_info *) tes->photoninfo)->keywords,&status);

    // initialize stream writer
    if (strcmp(par.trigger,"stream")==0) {
      tes->streaminfo=tes_init_tesrecord(par.tstart,par.tstop,tes,
					 par.streamfile,par.impactlist,
					 par.clobber,
					 keywords,
					 &status);
      CHECK_STATUS_BREAK(status);
      tes->write_to_stream=&tes_append_tesrecord;
      tes->write_photon=&tes_append_photon_tesrecord;
    } else {
	double threshold;
	unsigned int npts;
	unsigned int suppress;
	int strategy=-1;
	if (strncmp(par.trigger,"movavg:",7)==0) {
	  strategy=TRIGGER_MOVAVG;
	  
	  // threshold is float after movavg
	  sscanf(par.trigger,"movavg:%u:%lf:%u",&npts,&threshold,&suppress);

	  printf("\nChoosing movavg trigger strategy: %u points with threshold %f\n",npts,threshold);
	  printf("Trigger suppression interval: %u\n\n",suppress);
	  
	  assert(npts<20);
	  assert(threshold>0.0);
	  
	} else {
	  if (strncmp(par.trigger,"diff:",5)==0) {
	    strategy=TRIGGER_DIFF;
	    // threshold is float after diff
	    sscanf(par.trigger,"diff:%u:%lf:%u",&npts,&threshold,&suppress);
	    
	    printf("\nChoosing differential trigger of grade %u and threshold %f\n",npts,threshold);
	    printf("Trigger suppression interval: %u\n\n",suppress);
	    assert(threshold>0.0);
	    
	    unsigned int points[]={2,3,4,5,6,7,8,9,10};
	    int found=0;
	    for (int i=0; i<9; i++) {
	      if (points[i]==npts) {
		found=1;
	      }
	    }
	    if (found==0) {
	      fprintf(stderr,"Error: Grade of differential trigger must be one of ");
	      char comma=' ';
	      for (unsigned int i=0; i<sizeof(points)/sizeof(unsigned int); i++) {
		printf("%c %ui",comma,points[i]);
		comma=',';
	      }
	      exit(1);
	    }
	  } else {
	    fprintf(stderr,"Trigger strategy must be one of stream, movavg or diff\n");
	    exit(1);
	  }
	}

	tes->streaminfo=tes_init_trigger(par.tstart,par.tstop,tes,strategy,
					 par.preBufferSize,
					 par.triggerSize,threshold,npts,suppress,
					 par.streamfile,
					 par.impactlist,
					 par.clobber,
					 keywords,
					 &status);
	CHECK_STATUS_BREAK(status);
	tes->write_to_stream=&tes_append_trigger;
	tes->write_photon=NULL;
    }

    // Run the simulation

    tes_propagate(tes,par.tstop,&status);
    if (tes->progressbar!=NULL) {
      progressbar_finish(tes->progressbar);
    }
    CHECK_STATUS_BREAK(status);

    // write the file and free all allocated memory
    if (strcmp(par.trigger,"stream")==0) {
      tes_close_tesrecord(tes,&status);
      tes_free_tesrecord((tes_record_info **) &tes->streaminfo, &status);
    } else {
      tes_close_trigger(tes,&status);
      tes_free_trigger((tes_trigger_info **) &tes->streaminfo, &status);
    }

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
    query_simput_parameter_bool("acbias", &par->acdc, status);

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

    if (par->acdc) {
      query_simput_parameter_double("Rparasitic", &(par->Rpara), status); // mOhm
      assert(par->Rpara>=0);
      par->Rpara*=1e-3; // mOhm->Ohm
      par->RL=0;
    } else {
      query_simput_parameter_double("RL", &(par->RL), status); // mOhm
      assert(par->RL>=0);
      par->RL*=1e-3; // mOhm->Ohm
      par->Rpara=0.;
    } 

    if (par->acdc) {
      query_simput_parameter_double("TTR", &(par->TTR), status); 
      assert(par->TTR>=0);
    } else {
      par->TTR=0.;
    }

    query_simput_parameter_double("alpha", &(par->alpha), status);
    assert(par->alpha>0);
    query_simput_parameter_double("beta", &(par->beta), status);

    if (par->acdc) {
      query_simput_parameter_double("Lfilter", &(par->Lfilter), status); //muH
      assert(par->Lfilter>0);
      par->Lfilter*=1e-6;
      par->Lin=0.;
    } else {
      query_simput_parameter_double("Lin", &(par->Lin), status); //nH
      assert(par->Lin>0);
      par->Lin*=1e-9;
      par->Lfilter=0.;
    }
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

  // get optional effective voltage bias
  query_simput_parameter_double("V0",&(par->V0),status);
  par->V0*=1e-6;

  query_simput_parameter_bool("simnoise", &par->simnoise, status);

  // trigger strategy
  par->trigger=malloc(MAXMSG*sizeof(char));
  query_simput_parameter_string_buffer("triggertype",par->trigger,MAXMSG,status);
  if (strcmp(par->trigger,"stream")==0) {
    par->triggerSize=0;
    par->preBufferSize=0;
  } else {
    long dummy;
    query_simput_parameter_long("triggersize",&dummy,status);
    par->triggerSize=(unsigned long) dummy;
    query_simput_parameter_long("prebuffer",&dummy,status);
    // To Do: proper error handling for this condition
    par->preBufferSize=(unsigned long) dummy;
    assert(par->triggerSize > par->preBufferSize);
  } 

  // read for all input parameter possibilities
  query_simput_parameter_bool("simnoise", &par->simnoise, status);

  long seed; // needed because par->seed is an unsigned long
  query_simput_parameter_long("seed", &seed, status);
  par->seed=labs(seed); 

      

  query_simput_parameter_bool("propertiesonly",properties,status);

  int progressbar;
  query_simput_parameter_bool("progressbar",&progressbar,status);
  if (progressbar) {
    par->progressbar=progressbar_new("simulating",
				     (unsigned long) ((par->tstop-par->tstart)*PROGRESSBAR_FACTOR));
  } else {
	par->progressbar=NULL;
  }


  query_simput_parameter_bool("clobber", &par->clobber, status);

}
