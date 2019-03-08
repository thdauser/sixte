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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "pulsetemplgen.h"

////////////////////////////////////
/** Main procedure. */
int pulsetemplgen_main() {

  // Containing all programm parameters read by PIL.
  struct Parameters par;

  // Template structure
  TESProfiles* ptemp=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("pulsetemplgen");
  set_toolversion("0.05");

  do { // Beginning of the ERROR handling loop (will at
       // most be run once).

    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);

    genTESProfile(&par.pinp, &ptemp, &status);
    CHECK_STATUS_BREAK(status);

    int ii;
    for(ii=0; ii<ptemp->Nv; ii++){
      writeTESProfiles(par.filename,
			ptemp->version[ii],
			&ptemp->profiles[ii],
			par.clobber,
			par.history,
			&status);
      CHECK_STATUS_BREAK(status);
    }
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  destroyTESProfiles(ptemp);
  freePar(&par);

  return status;

}

void freePar(struct Parameters *par){

  if(par->pinp.version!=NULL){
    if(par->nver!=0){
      int ii;
      for(ii=0; ii<par->nver; ii++){
	if(par->pinp.version[ii]!=NULL){
	  free(par->pinp.version[ii]);
	  par->pinp.version[ii]=NULL;
	}
      }
    }
    free(par->pinp.version);
    par->pinp.version=NULL;
  }
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // set pointer to NULL
  par->pinp.version=NULL;
  par->pinp.energies=NULL;

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

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  par->nver=1;
  par->pinp.version=(char**)malloc(sizeof(char*));
  if(par->pinp.version==NULL){
    status=EXIT_FAILURE;
    SIXT_ERROR("Malloc error");
    return(status);
  }
  int ii;
  for(ii=0; ii<par->nver; ii++){
    par->pinp.version[ii]=(char*)malloc(9*sizeof(char));
    if(par->pinp.version[ii]==NULL){
      status=EXIT_FAILURE;
      SIXT_ERROR("Malloc error");
      return(status);
    }
  }
  status=ape_trad_query_string("Version", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the version name");
    return(status);
  }
  strcpy(par->pinp.version[0], sbuffer);
  free(sbuffer);


  status=ape_trad_query_double("freq", &par->pinp.freq);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the freq parameter");
    return(status);
  }

  status=ape_trad_query_double("energy_low", &par->energy_low);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy_low parameter");
    return(status);
  }

  status=ape_trad_query_double("energy_high", &par->energy_high);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy_high parameter");
    return(status);
  }

  status=ape_trad_query_int("energy_steps", &par->energy_steps);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the energy_steps parameter");
    return(status);
  }

  par->pinp.energies=(double*)malloc(par->energy_steps*sizeof(double));
  if(par->pinp.energies==NULL){
    status=EXIT_FAILURE;
    SIXT_ERROR("Malloc error");
    return(status);
  }
  for(ii=0; ii<par->energy_steps; ii++){
    par->pinp.energies[ii]=par->energy_low+ii*(par->energy_high-par->energy_low)/(double)(par->energy_steps-1);
  }
  par->pinp.ne=par->energy_steps;

  status=ape_trad_query_double("trise", &par->pinp.trise);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the trise parameter");
    return(status);
  }

  status=ape_trad_query_double("tfall", &par->pinp.tfall);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the tfall parameter");
    return(status);
  }

  status=ape_trad_query_long("nsamp", &par->pinp.nsamp);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the nsamp parameter");
    return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  return status;
}
