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


   Copyright 2014 Jelle de Plaa, SRON, Thorsten Brand, FAU
*/

#include "pulsetemplimport.h"

////////////////////////////////////
/** Main procedure. */
int pulsetemplimport_main() {
    
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Template structure
  TESProfiles *prof=NULL;
  
  // Error status.
  int status=EXIT_SUCCESS;
  
  int ii;
  
  // Register HEATOOL:
  set_toolname("pulsetemplimport");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at 
       // most be run once).
       
    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);
    prof=newTESProfiles(&status);
    CHECK_STATUS_BREAK(status);
    
    // allocate memory for data structures
    prof->Nv=1;
    prof->version=(char**)malloc(prof->Nv*sizeof(char*));
    if(prof->version==NULL){
      status=EXIT_FAILURE;
      CHECK_STATUS_BREAK(status);
    }
    
    for(ii=0; ii<prof->Nv; ii++){
      prof->version[ii]=(char*)malloc(9*sizeof(char));
      if(prof->version[ii]==NULL){
	status=EXIT_FAILURE;
	CHECK_STATUS_BREAK(status);
      }
    }
    
    strcpy(prof->version[0], par.version);
    
    prof->profiles=(TESProfilesEntries*)malloc(prof->Nv*sizeof(TESProfilesEntries));
    if(prof->profiles==NULL){
      status=EXIT_FAILURE;
      CHECK_STATUS_BREAK(status);
    }
    newTESProfilesEntries(prof->profiles);
    
    // read ascii table with pulse profile
    read_ascii_pulse(par.ascii, 
		      prof->profiles,
		      par.energy,
		      &status);
    CHECK_STATUS_BREAK(status);
    
    // write pulse profile into fits table
    
    writeTESProfiles(par.filename, 
			prof->version[0], 
			&prof->profiles[0], 
			0,
			"Modified with pulsetemplimport (SIXTE tool).",
			&status);
    CHECK_STATUS_BREAK(status);
    
  } while(0); // END of the error handling loop.
 
  
  return status;
  
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS; 

  // Read all parameters via the ape_trad_ routines.
  status=ape_trad_query_string("Filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the Filename");
    return(status);
  } 
  strcpy(par->filename, sbuffer);
  free(sbuffer);
   
  status=ape_trad_query_string("ASCIIFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the ASCII Filename");
    return(status);
  } 
  strcpy(par->ascii, sbuffer);
  free(sbuffer);
  
  status=ape_trad_query_string("Version", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the version name");
    return(status);
  }
  strcpy(par->version, sbuffer);
  free(sbuffer);
  
  
  status=ape_trad_query_double("energy", &par->energy);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy parameter");
    return(status);
  }
  
  return status;
}

void read_ascii_pulse(char *ascii, 
		      TESProfilesEntries *prof,
		      double energy,
		      int* const status){
 
  FILE *in=fopen(ascii, "r");
  if(in==NULL){
    SIXT_ERROR("Failed opening input ASCII file.");
    *status=EXIT_FAILURE;
    return ;
  }
  
  double t=0., p=0.;
  prof->Nt=0;
  prof->NE=1;
  prof->energy=(double*)malloc(sizeof(double));
  if(prof->energy==NULL){
    SIXT_ERROR("Malloc error.");
    *status=EXIT_FAILURE;
    return ;
  }
  prof->energy[0]=energy;
  
  prof->adc_value=(double**)malloc(sizeof(double*));
  if(prof->adc_value==NULL){
    SIXT_ERROR("Malloc error.");
    *status=EXIT_FAILURE;
    return ;
  }
  prof->adc_value[0]=NULL;
  prof->time=NULL;
  
  
  while(fscanf(in, "%le %le", &t, &p)!=EOF){
    prof->Nt++;
    prof->adc_value[0]=(double*)realloc(prof->adc_value[0], prof->Nt*sizeof(double));
    if(prof->adc_value[0]==NULL){
      SIXT_ERROR("Malloc error.");
      *status=EXIT_FAILURE;
      return ;
    }
    
    prof->time=(double*)realloc(prof->time, prof->Nt*sizeof(double));
    if(prof->time==NULL){
      SIXT_ERROR("Malloc error.");
      *status=EXIT_FAILURE;
      return ;
    }
    
    prof->time[prof->Nt-1]=t;
    prof->adc_value[0][prof->Nt-1]=p;
       
  }
  
  fclose(in);
}