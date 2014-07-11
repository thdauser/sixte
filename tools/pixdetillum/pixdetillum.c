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


   Copyright 2014 Thorsten Brand, FAU
*/



#include "pixdetillum.h"

/////////////////////////////////
/** random number generator */



/** Use the C rand() function to determine a random number between 0
    and 1. */
static inline double getRndNum(int* const status) 
{
  double r=(double)rand()/((double)RAND_MAX+1.0);
  assert(r<1.0);
  return(r);

  // Status variable is not needed.
  (void)(*status);
}

//////////////////////////////////
/** main procedure */
int pixdetillum_main(){
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Detector data structure.
  AdvDet *det=NULL;
  
  // Pixel impact list
  PixImpFile* plf=NULL;
  
  // Error status.
  int status=EXIT_SUCCESS;
  
  long ii;
  
  // Register HEATOOL:
  set_toolname("pixdetillum");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at 
       // most be run once).
       
    // Read parameters using PIL library.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    headas_chat(3, "initialize ...\n");
    
    // Set the rng seed
    sixt_init_rng(par.seed, &status);
    
    // Load detector information
    det=loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    
    if(par.Nactive==-1){
      par.Nactive=det->npix;
      par.nlo=0;
      par.nhi=det->npix-1;
    }
    
    // Determine the event list output file.
    char piximplist_filename[MAXFILENAME];
    strcpy(piximplist_filename, par.PixImpList);
    
    // Create origin impactlist_filename
    char impactlist_filename[]="pixdetillum";
    
    // Open the output file
    plf=openNewPixImpFile(piximplist_filename,
			   par.telescop,
			   par.instrume,
			   par.filter,
			   par.ancrfile,
			   par.respfile,
			   par.XMLFile,
			   impactlist_filename,
			   par.mjdref,
			   par.timezero,
			   par.tstart,
			   par.tstop,
			   par.clobber,
			   &status);
    CHECK_STATUS_BREAK(status);
    
    // Generate random impacts
    
    double time=0.;
    ii=0;
      
    while(time<=par.tstop){
      // Determine random impact time
      double u=sixt_get_random_number(&status);
      time=time-log(1.-u)/par.rate/par.Nactive;
      // determine random pixel
      long pixid=(long)(sixt_get_random_number(&status)*par.Nactive)+(long)par.nlo;
      // determine photon attributes
      PixImpact piximp;
      piximp.pixID=pixid;
      piximp.time=time;
      piximp.energy=(float)par.energy;
      piximp.ph_id=ii;
      piximp.src_id=0;
      double x=sixt_get_random_number(&status);
      double y=sixt_get_random_number(&status);
      piximp.pixposition.x=(-0.5+x)*det->pix[pixid].width;
      piximp.pixposition.y=(-0.5+y)*det->pix[pixid].height;
      piximp.detposition.x=piximp.pixposition.x+det->pix[pixid].sx;
      piximp.detposition.y=piximp.pixposition.y+det->pix[pixid].sy;
      
      // save impact      
      addImpact2PixImpFile(plf, &piximp, &status);
      ii++;
    }
    
  } while(0); // END of the error handling loop.
  
  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");
  freePixImpFile(&plf, &status);
 
  destroyAdvDet(&det);
  
  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }  
  
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.

  status=ape_trad_query_string("PixImpList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the pixel impact list");
    return(status);
  } 
  strcpy(par->PixImpList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the XML file");
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Telescop", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the telescope");
    return(status);
  } 
  strcpy(par->telescop, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Instrume", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the instrument");
    return(status);
  } 
  strcpy(par->instrume, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("Filter", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the filter");
    return(status);
  } 
  strcpy(par->filter, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("ancrfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the ancrfile");
    return(status);
  } 
  strcpy(par->ancrfile, sbuffer);
  free(sbuffer); 

  status=ape_trad_query_string("respfile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the respfile");
    return(status);
  } 
  strcpy(par->respfile, sbuffer);
  free(sbuffer); 
  
  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstart", &par->tstart);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstart parameter");
    return(status);
  }
  
  status=ape_trad_query_double("tstop", &par->tstop);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tstop parameter");
    return(status);
  }
  
  status=ape_trad_query_double("timezero", &par->timezero);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the timezero parameter");
    return(status);
  }
  
  status=ape_trad_query_double("mjdref", &par->mjdref);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the mjdref parameter");
    return(status);
  }
  
  status=ape_trad_query_double("rate", &par->rate);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the rate parameter");
    return(status);
  }
  
  status=ape_trad_query_double("energy", &par->energy);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy parameter");
    return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  int seed=0;
  status=ape_trad_query_int("Seed", &seed);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the seed.");
    return(status);
  }
  
  if(seed==-1){
    // Initialize with system time
    struct timeval tv;
    gettimeofday(&tv,NULL);
    par->seed=1000000*tv.tv_sec+tv.tv_usec;
  }else{
    par->seed=(unsigned long int)seed;
  }
  printf("Seed=%ld\n", par->seed);
    
  char *pix=NULL;
  status=ape_trad_query_string("pixels", &pix);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the pixel parameter");
    return(status);
  }
  if(!strcmp(pix, "all")){
    par->Nactive=-1;
  } else {
    char *ptr=strchr(pix, '-');
    if(ptr!=NULL){
      sscanf(pix, "%d - %d", &(par->nlo), &(par->nhi));
      par->nlo--;
      par->nhi--;
      par->Nactive=par->nhi-par->nlo+1;
      printf("Activate pixels %d-%d (%d pixels).\n", par->nlo+1, par->nhi+1, par->Nactive);
    } else {
      sscanf(pix, "%d", &(par->nlo));
      par->nlo--;
      par->nhi=par->nlo;
      par->Nactive=1;
    }
    if(par->Nactive<0 || par->nlo<0){
      status=EXIT_FAILURE;
      SIXT_ERROR("Input values for active pixels corrupted.");
      return(status);
    }
  }
  free(pix);

  return(status);
}