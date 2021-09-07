/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

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

void tessim_getpar(tespxlparams *par, AdvDet **det, int *properties, int *status);


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

    // advanced detector structure
    AdvDet* det=NULL;

    tessim_getpar(&par,&det,&properties,&status);
    CHECK_STATUS_BREAK(status);

    if (det->npix>1 && par.dobbfb){
      SIXT_ERROR("Multipixel simulation and BBFB cannot be simulated together yet");
      return EXIT_FAILURE;
    }

    /** LOOP FOR MULTI TESSIM START */
    tespxlparams loop_par; // parameter to account for different pixel parameters during loop
    for (int ii=0;ii<det->npix;ii++) {
      // initialize the TES pixel

      // Initialization for the 1 pixel case
      if (det->npix==1){
        loop_par = par;
        det->pix[ii].tes=tes_init(&par,&status);
        tesparams *tes = det->pix[ii].tes;
        CHECK_STATUS_BREAK(status);

        // output information about simulation
        tes_print_params(tes);

        // only for npix=1
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
          return(EXIT_SUCCESS);
        }
      } else {
        // get file+extension name for this pixel
        char* filename = (char*)malloc(strlen(det->tes_type_file)+strlen(det->pix[ii].tes_type)+3);
        strcpy(filename, det->tes_type_file);
        // +3: zero-terminator, "[" and "]"
        strncat(filename,"[",1);
        strcat(filename,det->pix[ii].tes_type);
        strncat(filename,"]",1);
        // read loop_par
        tes_fits_read_params(filename,&loop_par,&status);
        free(filename);
        // Copy relevant parameters from par to loop_par
        loop_par.tstart = par.tstart;
        loop_par.tstop = par.tstop;
        loop_par.acdc=par.acdc;
        loop_par.seed = par.seed;
        loop_par.simnoise = par.simnoise;
        loop_par.m_excess = par.m_excess;
        loop_par.thermal_bias = par.thermal_bias;
        loop_par.readoutMode = par.readoutMode;
        loop_par.twofluid = par.twofluid;
        loop_par.stochastic_integrator = par.stochastic_integrator;
        loop_par.frame_hit=0; //Setting to default false //TODO add to multi-tessim
        loop_par.dobbfb=0; // Disabled for now in multitessim

        // assign the actual pixid, input and output filename
        loop_par.id = det->pix[ii].pindex + 1; // id's start with 0 in advdet, but not in tessim...
        // turn the pixid into a string
        int id_string_length=snprintf(NULL,0,"%d",loop_par.id);
        char* id_string=malloc(id_string_length+1);
        snprintf(id_string, id_string_length+1, "%d", loop_par.id);
        // impact file
        char* impactbuffer=(char*)malloc(strlen(par.impactlist)+id_string_length+10);
        strcpy(impactbuffer, par.impactlist);
        strncat(impactbuffer, "[PIXID==",8);
        strncat(impactbuffer, id_string, id_string_length);
        strncat(impactbuffer, "]", 1);
        loop_par.impactlist = strdup(impactbuffer);
        free(impactbuffer);

        // output file - in the xml case, the input is just a prefix!
        char* streamfilebuffer= (char*)malloc(strlen(par.streamfile)+id_string_length+10);
        strcpy(streamfilebuffer,par.streamfile);
        strncat(streamfilebuffer, "_pix", 4);
        strncat(streamfilebuffer, id_string, id_string_length);
        strncat(streamfilebuffer, ".fits",5);
        loop_par.streamfile = strdup(streamfilebuffer);
        free(streamfilebuffer);

        free(id_string);
        // initialize the tes
        det->pix[ii].tes=tes_init(&loop_par,&status);
        CHECK_STATUS_BREAK(status);
      }

      tesparams *tes = det->pix[ii].tes;
      // initialize photon provider
      tes->photoninfo=tes_init_impactlist(loop_par.impactlist,&status);
      CHECK_STATUS_BREAK(status);
      tes->get_photon=&tes_photon_from_impactlist; // note: function pointer!

      SixtStdKeywords *keywords=duplicateSixtStdKeywords(((tes_impactfile_info *) tes->photoninfo)->keywords,&status);

      // initialize stream writer
      if (strcmp(par.trigger,"stream")==0) {
        tes->streaminfo=tes_init_tesrecord(par.tstart,par.tstop,tes,TESRECORD_BUFFERSIZE/det->npix,
                                           loop_par.streamfile,loop_par.impactlist,
                                           par.clobber,par.write_error,
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

            // print only once
            if (ii==0){
              printf("\nChoosing movavg trigger strategy: %u points with threshold %f\n",npts,threshold);
              printf("Trigger suppression interval: %u\n\n",suppress);
            }

            assert(npts<20);
            assert(threshold>0.0);

          } else {
            if (strncmp(par.trigger,"diff:",5)==0) {
              strategy=TRIGGER_DIFF;
              // threshold is float after diff
              sscanf(par.trigger,"diff:%u:%lf:%u",&npts,&threshold,&suppress);
              // print only once
              if (ii==0){

                printf("\nChoosing differential trigger of grade %u and threshold %f\n",npts,threshold);
                printf("Trigger suppression interval: %u\n\n",suppress);
              }
              assert(threshold>0.0);

              unsigned int points[]={2,3,4,5,6,7,8,9,10};
              int found=0;
              for (int i=0; i<9; i++) {
                if (points[i]==npts) {
                  found=1;
                }
              }
              if (found==0) {
                // exit condition, so it'll only print once in any case
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
                // print only once
                if (ii==0){
                  printf("\nChoosing noise output with record length %u\n",npts);
                }
              } else {
                if (strncmp(par.trigger,"impact",6)==0) {
                  assert(par.impactlist != NULL); // Can't trigger on impact without an impactlist
                  strategy=TRIGGER_IMPACT;
                  threshold=0.;
                  npts=0;
                  suppress=par.triggerSize;
                } else {
                  fprintf(stderr,"Trigger strategy must be one of stream, movavg, diff, noise or impact\n");
                  exit(1);
                }
              }
            }
          }

          tes->streaminfo=tes_init_trigger(par.tstart,par.tstop,tes,strategy,
                                           par.preBufferSize,
                                           par.triggerSize,threshold,npts,suppress,
                                           loop_par.streamfile,
                                           loop_par.impactlist,
                                           par.clobber,
                                           keywords,
                                           par.write_doubles,
                                           &status);
          CHECK_STATUS_BREAK(status);
          tes->write_to_stream=&tes_append_trigger;
          tes->write_photon=NULL;
      }
      freeSixtStdKeywords(keywords);

      // Initialize BBFB loop
      // For architecture reasons, this is for the moment at tes level whereas it should actually be at channel level
      // but this structure is not completely available yet
      if (tes->dobbfb) {
        tes->bbfb_info = init_timedomain_bbfb(par.bbfb_tclock,par.bbfb_delay,par.carrier_frequency,par.bbfb_gbw,0,
        		par.bias_leakage*tes->I0_start,par.bias_leakage_phase,par.fb_leakage,par.fb_leakage_lag,tes->I0_start,&status);
        tes->apply_bbfb = &run_timedomain_bbfb_loop;
      }

    }
    /** LOOP FOR MULTI TESSIM END */





    /** SET UP CROSSTALK */
    if (det->npix > 1 && par.doCrosstalk && par.acdc) {
      // Initialize the electrical FDM crosstalk
      if (det->readout_channels==NULL){
        det->readout_channels = get_readout_channels(det,&status);
        if (status!= EXIT_SUCCESS){
          SIXT_ERROR("failed when loading the readout channels");
          return EXIT_FAILURE;
        }
      }

      // get the pixel resonance frequencies
      get_resonance_frequencies(det, &status);
      if (status!= EXIT_SUCCESS){
          SIXT_ERROR("failed when loading the resonance frequencies");
          return EXIT_FAILURE;
      }

      // set up readout FDMSystem for each readout channel
      for (int ii=0; ii<det->readout_channels->num_channels; ii++){
        // assuming same TTR for all pixels
        double TTR = det->readout_channels->channels[ii].pixels[0]->tes->TTR;
        init_FDMSystem(&(det->readout_channels->channels[ii]), det->L_Common, det->C_Common, TTR, &status);
      }

      if (status!=EXIT_SUCCESS){
        SIXT_ERROR("failed when initializing electrical crosstalk");
        return EXIT_FAILURE;
      }
      headas_chat(0,"Initialized FDM Crosstalk\n");
    }



    // set up the progress bar, but only for one pixel
    if (par.showprogress) {
      det->pix[0].tes->progressbar=progressbar_new("simulating",
                                       (unsigned long) ((par.tstop-par.tstart)*PROGRESSBAR_FACTOR));
    }
    // Run the simulation
    tes_propagate(det,par.tstop,&status);

    /** LOOP FOR MULTI TESSIM START */
    for (int ii=0;ii<det->npix;ii++) {
      // get this pixel's tesparams pointer
      tesparams *tes=det->pix[ii].tes;

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
    }

    /** LOOP FOR MULTI TESSIM END */

    destroyAdvDet(&det);
    free(par.type);
    free(par.impactlist);
    free(par.streamfile);
    free(par.trigger);
    if (det->npix >1){
      free(loop_par.type);
      free(loop_par.impactlist);
      free(loop_par.streamfile);
    }
//    freeSixtStdKeywords(keywords);

  } while(0); // end of error handling loop

  //printf("%s %d %s","Paso 5: status=",status,"\n");
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
    if (strcasecmp(parname,query_cmdpars[i])==0) {
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

int cmd_query_simput_parameter_file_name(int fromcmd,char *name, char **field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_file_name(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_string(int fromcmd,char *name, char **field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_string(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_string_buffer(int fromcmd,char *name, char * const field, int buflen, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_string_buffer(name,field,buflen,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_int(int fromcmd,char *name, int *field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_int(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_long(int fromcmd,char *name, long *field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_long(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_float(int fromcmd,char *name, float *field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_float(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_double(int fromcmd,char *name, double *field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_double(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}

int cmd_query_simput_parameter_bool(int fromcmd,char *name, int *field, int *status) {
  if (fromcmd || sixt_par_on_commandline(name)) {
    query_simput_parameter_bool(name,field,status);
    if (!fromcmd) {
      printf("Overriding parameter %s from command line\n",name);
    }
    return 1;
  }
  return 0;
}


void tessim_getpar(tespxlparams *par, AdvDet **det, int *properties, int *status) {
  query_simput_parameter_string("PixType", &(par->type), status);

  int fromcmd=1;

  query_simput_parameter_file_name("PixImpList", &(par->impactlist), status);
  query_simput_parameter_file_name("Streamfile", &(par->streamfile), status);
  query_simput_parameter_double("tstart", &(par->tstart), status);
  query_simput_parameter_double("tstop", &(par->tstop), status);

  // if par->type starts with "xml:", load the advdet
  if (strncmp(par->type,"xml:",4)==0) {
    char *file=strdup(par->type+4);
    (*det)  = loadAdvDet(file,status);
    CHECK_STATUS_VOID(*status);
    free(file);
    fromcmd=1;
    sixt_init_query_commandline(status);
    CHECK_STATUS_VOID(*status);
    cmd_query_simput_parameter_bool(fromcmd,"doCrosstalk", &(par->doCrosstalk), status);
  } else {
    // create a new advdet, setting nr of pixels to 1
    *det = newAdvDet(status);
    CHECK_STATUS_VOID(*status);
    (*det)->npix=1;
    (*det)->pix=newAdvPix(status);
    CHECK_STATUS_VOID(*status);
    // there is no crosstalk
    par->doCrosstalk=0;
  }



  if ((*det)->npix==1){
  query_simput_parameter_int("PixID",&(par->id),status);
    if (strncmp(par->type,"file:",5)==0) {
      char *file=strdup(par->type+5);
      //strcpy(det->tes_type_file,file);
      tes_fits_read_params(file,par,status);
      CHECK_STATUS_VOID(*status);
      free(file);
      fromcmd=0;
      sixt_init_query_commandline(status);
      CHECK_STATUS_VOID(*status);
    }
  } else {
    // set par->id to a dummy value
    par->id = 0;
  }

  //
  // The following calls read the parameter value in two cases:
  // a) parameters have NOT been read from a FITS file
  // or
  // b) parameters have been read from a FITS file, BUT they are given
  //    on the command line
  // The code is fast since the cmd_... routines make use of short-circuiting of the ||-operator
  cmd_query_simput_parameter_bool(fromcmd,"acbias", &par->acdc, status);
  cmd_query_simput_parameter_double(fromcmd,"sample_rate",&(par->sample_rate),status);

  assert(par->sample_rate>0);

  // trigger strategy
  par->trigger=malloc(MAXMSG*sizeof(char));
  query_simput_parameter_string_buffer("triggertype",par->trigger,MAXMSG,status);
  if (strcmp(par->trigger,"stream")==0) {
    par->triggerSize=0;
    par->preBufferSize=0;
    par->write_doubles=0;
  } else {
    long dummy;
    query_simput_parameter_long("triggersize",&dummy,status);
    par->triggerSize=(unsigned long) dummy;
    query_simput_parameter_long("prebuffer",&dummy,status);
    // To Do: proper error handling for this condition
    par->preBufferSize=(unsigned long) dummy;
    assert(par->triggerSize > par->preBufferSize);
    query_simput_parameter_bool("write_doubles", &par->write_doubles, status);
  }

  long seed; // needed because par->seed is an unsigned long
  query_simput_parameter_long("seed", &seed, status);
  par->seed=labs(seed);

  query_simput_parameter_bool("propertiesonly",properties,status);

  query_simput_parameter_bool("progressbar",&par->showprogress,status);

  query_simput_parameter_bool("clobber", &par->clobber, status);

  // query parameter for calculating I0 via the thermal balance
  query_simput_parameter_bool("thermalBias", &par->thermal_bias, status);

  // modifying pixel parameters from command line is only allowed for one-pixel runs
  if ((*det)->npix==1){
    if (cmd_query_simput_parameter_double(fromcmd,"Ce",&(par->Ce1),status) ) {
      par->Ce1*=1e-12; // pJ/K -> J/K
    }
    assert(par->Ce1>0);

    if (cmd_query_simput_parameter_double(fromcmd,"Gb",&(par->Gb1),status) ) {
      par->Gb1*=1e-12; // pW/K -> W/K
    }
    assert(par->Gb1>0);

    if (cmd_query_simput_parameter_double(fromcmd,"T_start", &(par->T_start), status) ) {
      par->T_start*=1e-3; // mK -> K
    }
    assert(par->T_start>0);

    if (cmd_query_simput_parameter_double(fromcmd,"Tb", &(par->Tb), status) ) {
      par->Tb*=1e-3; // mK->K
    }
    assert(par->Tb>0);

    if (cmd_query_simput_parameter_double(fromcmd,"R0", &(par->R0), status) ) {
      par->R0*=1e-3; // mOhm->Ohm
    }
    assert(par->R0>0);

    if (cmd_query_simput_parameter_double(fromcmd, "bias", &(par->bias), status)){
        par->bias*=1e-2; // %->decimal
    }
    assert(par->bias>=0.); //Not 0
    assert(par->bias<=1.); //Below Rn

    par->Pload=0;
    if (cmd_query_simput_parameter_double(fromcmd, "Pload", &(par->Pload), status)){
        par->Pload*=1e-12; // pW->W
    }
    printf("Pload %.15f\n", par->Pload);

    if (cmd_query_simput_parameter_double(fromcmd,"I0", &(par->I0), status) ) {
      // if this is set simultaneously to the thermal balance, this throws an error
      if (par->thermal_bias){
        SIXT_ERROR("Conflict: I0 is both given as parameter and should be automatically calculated");
        *status=EXIT_FAILURE;
        return;
      }
      par->I0*=1e-6; // muA->A
    }
    assert(par->I0>0);

    // get optional effective voltage bias
    if (cmd_query_simput_parameter_double(fromcmd,"V0", &(par->V0), status) ) {
      // if this is set simultaneously to the thermal balance, this throws an error
      if (par->thermal_bias){
        SIXT_ERROR("Conflict: V0 is both given as parameter and should be automatically calculated");
        *status=EXIT_FAILURE;
        return;
      }
      par->V0*=1e-6; // muV->V
    }

    if (par->acdc) {
      if (cmd_query_simput_parameter_double(fromcmd,"Rparasitic", &(par->Rpara), status)) {
        par->Rpara*=1e-3; // mOhm->Ohm
      }
      par->RL=0;
      assert(par->Rpara>=0);
    } else {
      if (cmd_query_simput_parameter_double(fromcmd,"RL", &(par->RL), status)) {
        par->RL*=1e-3; // mOhm->Ohm
      }
      par->Rpara=0.;
      assert(par->RL>=0);
    }

    if (par->acdc) {
      cmd_query_simput_parameter_double(fromcmd,"TTR", &(par->TTR), status);
      assert(par->TTR>=0);
    } else {
      par->TTR=0.;
    }

    cmd_query_simput_parameter_double(fromcmd,"alpha", &(par->alpha), status);
    assert(par->alpha>0);
    cmd_query_simput_parameter_double(fromcmd,"beta", &(par->beta), status);

    if (par->acdc) {
      if (cmd_query_simput_parameter_double(fromcmd,"Lfilter", &(par->Lfilter), status)){
        par->Lfilter*=1e-6;
      }
      par->Lin=0.;
      assert(par->Lfilter>0);
    } else {
      if (cmd_query_simput_parameter_double(fromcmd,"Lin", &(par->Lin), status)) {
        par->Lin*=1e-9;
      }
      par->Lfilter=0.;
      assert(par->Lin>0);
    }

    cmd_query_simput_parameter_double(fromcmd,"n", &(par->n), status);
    cmd_query_simput_parameter_double(fromcmd,"imin", &(par->imin), status);
    cmd_query_simput_parameter_double(fromcmd,"imax", &(par->imax), status);

    if (par->imax<=par->imin) {
      SIXT_ERROR("imax MUST be larger imin\n");
      *status=EXIT_FAILURE;
      return;
    }

  }
  query_simput_parameter_bool("simnoise", &par->simnoise, status);
  query_simput_parameter_double("squidnoise",&par->squid_noise,status);
  query_simput_parameter_double("bias_noise", &par->bias_noise, status);

  if (&par->simnoise) {
    cmd_query_simput_parameter_double(fromcmd,"m_excess",&par->m_excess,status);
  } else {
    par->m_excess=0.; // switch explicitly off if not simulating noise
  }
  assert(par->m_excess>=0);

  query_simput_parameter_bool("twofluid", &par->twofluid, status);
  //TODO Change keyword name when new models become available
  // readout mode is only not 'total' for crosstalk
  query_simput_parameter_bool("stochastic_integrator", &par->stochastic_integrator, status);

  //Handling the frame impacts
  query_simput_parameter_bool("frame_hit", &par->frame_hit, status);
  query_simput_parameter_file_name("frame_hit_file", &(par->frame_hit_file), status);
  query_simput_parameter_double("frame_hit_time", &(par->frame_hit_time), status);

  if (par->frame_hit!=0){
      if (NULL==par->frame_hit_file){
          SIXT_ERROR("An input file is needed for the frame hit!");
          *status=EXIT_FAILURE;
          return;
      }
  }

  if (par->frame_hit!=0){
	  printf("Reading frame hit template from %s, impact time %f s \n", par->frame_hit_file, par->frame_hit_time);
  }

  if (!par->doCrosstalk){
    par->readoutMode=READOUT_TOTAL;
  } else {
    char mode[MAXMSG];
    query_simput_parameter_string_buffer("readoutMode", mode, MAXMSG, status);
    if (strncmp(mode,"total",5) == 0) {
      par->readoutMode = READOUT_TOTAL;
    } else if (strncmp(mode,"I",1) == 0) {
      par->readoutMode = READOUT_ICHANNEL;
    } else if (strncmp(mode,"Q",1) == 0) {
      par->readoutMode = READOUT_QCHANNEL;
    } else {
      SIXT_ERROR("Output mode must be one of 'total', 'I' or 'Q'");
      *status=EXIT_FAILURE;
      return;
    }
  }

  // read bbfb parameters
  query_simput_parameter_bool("dobbfb", &par->dobbfb, status);
  query_simput_parameter_bool("decimation_filter", &par->decimation_filter, status);
  query_simput_parameter_double("bbfb_tclock", &par->bbfb_tclock, status);
  query_simput_parameter_double("carrier_frequency", &par->carrier_frequency, status);
  query_simput_parameter_int("bbfb_delay", &par->bbfb_delay, status);
  query_simput_parameter_double("bbfb_gbw", &par->bbfb_gbw, status);
  query_simput_parameter_double("bias_leakage", &par->bias_leakage, status);
  query_simput_parameter_double("bias_leakage_phase", &par->bias_leakage_phase, status);
  if (!par->simnoise && par->bias_leakage_phase==-999){
	par->bias_leakage_phase=0;
	printf("\nWARNING: bias_leakage_phase is supposed to be random but noiseless simulation requested -> phase arbitrarily fixed to 0\n");
  }
  query_simput_parameter_double("fb_leakage", &par->fb_leakage, status);
  query_simput_parameter_int("fb_leakage_lag", &par->fb_leakage_lag, status);
  if (!par->simnoise && par->fb_leakage_lag==-999){
	par->fb_leakage_lag=0;
	printf("\nWARNING: fb_leakage_lag is supposed to be random but noiseless simulation requested -> lag arbitrarily fixed to 0\n");
  }
  query_simput_parameter_double("M_in", &par->M_in, status);
  query_simput_parameter_bool("write_error", &par->write_error, status);


  sixt_free_query_commandline();

}
