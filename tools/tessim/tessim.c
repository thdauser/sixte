#include <assert.h>
#include <parinput.h>
#include "tessim.h"
#include <stdlib.h>


#define TOOLSUB tessim_main
//   #include "headas_main.c"
#include "sixt_main.c"

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
	    if (strncmp(par.trigger,"noise",5)==0) {
	      strategy=TRIGGER_NOISE;
	      npts=0;
	      suppress=par.triggerSize;
	      threshold=0.;
	      printf("\nChoosing noise output with record length %u\n",npts);
	    } else {
	      fprintf(stderr,"Trigger strategy must be one of stream, movavg, diff, or noise\n");
	      exit(1);
	    }
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
    free(par.trigger);

    freeSixtStdKeywords(keywords);
    
  } while(0); // end of error handling loop

  if (EXIT_SUCCESS==status) {
    headas_chat(3,"finished successfully!\n\n");
    return(EXIT_SUCCESS);
  }
  return(EXIT_FAILURE);
}


int query_argc=-1;
char **query_cmdpars=NULL;

//
// remember all parameter arguments on the command line
// (i.e. split individual argv at the equal sign if present and
// store them)
//
void sixt_init_query_commandline(int *status) {
  query_cmdpars=malloc(sixt_argc*sizeof(char *));
  query_argc=sixt_argc;
  CHECK_NULL_VOID(query_cmdpars,*status,"No memory to cache commandline");
  query_cmdpars[0]=strdup(sixt_argv[0]); // copy program name
  CHECK_NULL_VOID(query_cmdpars,*status,"No memory to cache argv[0]");
  for (int i=1;i<sixt_argc;i++) {
    char *equal=strstr(sixt_argv[i],"=");
    if (equal==NULL) {
      // "="-sign not found, duplicate the whole string and hope for the best
      query_cmdpars[i]=strdup(sixt_argv[i]);
    } else {
      query_cmdpars[i]=strndup(sixt_argv[i],equal-sixt_argv[i]);
    }
    CHECK_NULL_VOID(query_cmdpars[i],*status,"No memory to cache argv");
  }
}

//
// return true if parameter parname was given on the commandline
// (ignoring the case of the string)
//
int sixt_par_on_commandline(char *parname) {
  for (int i=1;i<query_argc;i++) {
    if (strcasecmp(parname,query_cmdpars[i])!=0) {
      return 1;
    }
  }
  return 0;
}

//
// release memory allocated by sixt_init_query_commandline
//
void sixt_free_query_commandline() {
  if (query_argc<0) {
    return;
  }
  for (int i=0;i<query_argc;i++) {
    free(query_cmdpars[i]);
  }
  free(query_cmdpars);
  query_cmdpars=NULL;
}


void tessim_getpar(tespxlparams *par, int *properties, int *status) {
  query_simput_parameter_string("PixType", &(par->type), status);

  int fromcmd=1;
  
  query_simput_parameter_int("PixID",&(par->id),status);
  query_simput_parameter_file_name("PixImpList", &(par->impactlist), status);
  query_simput_parameter_file_name("Streamfile", &(par->streamfile), status);
  query_simput_parameter_double("tstart", &(par->tstart), status);
  query_simput_parameter_double("tstop", &(par->tstop), status);

  if (strncmp(par->type,"file:",5)==0) {
    char *file=strdup(par->type+5);
    tes_fits_read_params(file,par,status);
    CHECK_STATUS_VOID(*status);
    free(file);
    fromcmd=0;
    sixt_init_query_commandline(status);
    CHECK_STATUS_VOID(*status);
  }

  //
  // The following calls read the parameter value in two cases:
  // a) parameters have NOT been read from a FITS file
  // or
  // b) parameters have been read from a FITS file, BUT they are given
  //    on the command line
  // The code is fast since it makes use of short-circuiting of the ||-operator
  //
  // We leave the asserts outside of the ifs in order to get some
  // sanity checking for the parameters from the FITS file, i.e.,
  // sixt_par_on_commandline is called ONLY if we got parameters from
  // the command line and have to check whether they are overridden
  //
  
  if (fromcmd || sixt_par_on_commandline("acbias")) {
    query_simput_parameter_bool("acbias", &par->acdc, status);
  }

  if (fromcmd || sixt_par_on_commandline("sample_rate")) {
    query_simput_parameter_double("sample_rate",&(par->sample_rate),status);
  }
  assert(par->sample_rate>0);

  if (fromcmd || sixt_par_on_commandline("Ce")) {
    query_simput_parameter_double("Ce",&(par->Ce1),status);
    par->Ce1*=1e-12; // pJ/K -> J/K
  }
  assert(par->Ce1>0);

  if (fromcmd || sixt_par_on_commandline("Gb")) {
    query_simput_parameter_double("Gb",&(par->Gb1),status);
    par->Gb1*=1e-12; // pW/K -> W/K
  }
  assert(par->Gb1>0);
  
  if (fromcmd || sixt_par_on_commandline("T_start")) {
    query_simput_parameter_double("T_start", &(par->T_start), status); //mK
    par->T_start*=1e-3; // mK -> K
  }
  assert(par->T_start>0);
  
  if (fromcmd || sixt_par_on_commandline("Tb")) {
    query_simput_parameter_double("Tb", &(par->Tb), status); // mK
    par->Tb*=1e-3; // mK->K
  }
  assert(par->Tb>0);
  
  if (fromcmd || sixt_par_on_commandline("R0")) {
    query_simput_parameter_double("R0", &(par->R0), status); // mOhm
    par->R0*=1e-3; // mOhm->Ohm
  }
  assert(par->R0>0);
  
  if (fromcmd || sixt_par_on_commandline("I0")) {
    query_simput_parameter_double("I0", &(par->I0), status); // muA
    par->I0*=1e-6; // muA->A
  }
  assert(par->I0>0);
  
  if (par->acdc) {
    if (fromcmd || sixt_par_on_commandline("Rparasitic")) {
      query_simput_parameter_double("Rparasitic", &(par->Rpara), status); // mOhm
      par->Rpara*=1e-3; // mOhm->Ohm
      par->RL=0;
    }
    assert(par->Rpara>=0);
  } else {
    if (fromcmd || sixt_par_on_commandline("RL")) {
      query_simput_parameter_double("RL", &(par->RL), status); // mOhm
      par->RL*=1e-3; // mOhm->Ohm
      par->Rpara=0.;
    }
    assert(par->RL>=0);
  } 

  if (par->acdc) {
    if (fromcmd || sixt_par_on_commandline("TRR")) {
      query_simput_parameter_double("TTR", &(par->TTR), status); 
      assert(par->TTR>=0);
    }
  } else {
    par->TTR=0.;
  }

  if (fromcmd || sixt_par_on_commandline("alpha")) {
    query_simput_parameter_double("alpha", &(par->alpha), status);
  }
  assert(par->alpha>0);
  if (fromcmd || sixt_par_on_commandline("beta")) {
    query_simput_parameter_double("beta", &(par->beta), status);
  }
  
  if (par->acdc) {
    if (fromcmd || sixt_par_on_commandline("Lfilter")) {
      query_simput_parameter_double("Lfilter", &(par->Lfilter), status); //muH
      par->Lfilter*=1e-6;
      par->Lin=0.;
    }
    assert(par->Lfilter>0);
  } else {
    if (fromcmd || sixt_par_on_commandline("Lin")) {
      query_simput_parameter_double("Lin", &(par->Lin), status); //nH
      par->Lin*=1e-9;
      par->Lfilter=0.;
    }
    assert(par->Lin>0);
  }

  if (fromcmd || sixt_par_on_commandline("n")) {
    query_simput_parameter_double("n", &(par->n), status);
  }
  if (fromcmd || sixt_par_on_commandline("imin")) {
    query_simput_parameter_double("imin", &(par->imin), status);
  }
  if (fromcmd || sixt_par_on_commandline("imax")) {
    query_simput_parameter_double("imax", &(par->imax), status);
  }

  if (par->imax<=par->imin) {
    SIXT_ERROR("imax MUST be larger imin\n");
    *status=EXIT_FAILURE;
    return;
  }

  if (fromcmd || sixt_par_on_commandline("simnoise")) {
    query_simput_parameter_bool("simnoise", &par->simnoise, status);
  }

  if (&par->simnoise) {
    if (fromcmd || sixt_par_on_commandline("m_excess")) {
      query_simput_parameter_double("m_excess",&par->m_excess,status);
    }
  } else {
    par->m_excess=0.; // switch explicitly off if not simulating noise
  }
  assert(par->m_excess>=0);

  // get optional effective voltage bias
  // NB: NOT stored in the FITS file (why?)
  query_simput_parameter_double("V0",&(par->V0),status);
  par->V0*=1e-6;


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

  sixt_free_query_commandline();
  
}
