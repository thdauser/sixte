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

#include "tesproftemplates.h"

void destroyTESProfilesEntries(TESProfilesEntries* prof){
  
  if(prof->adc_value!=NULL){
    
    int ii;
    for(ii=0; ii<prof->NE; ii++){
      if(prof->adc_value[ii]!=NULL){
	free(prof->adc_value[ii]);
	prof->adc_value[ii]=NULL;
      }
    }
    free(prof->adc_value);
    prof->adc_value=NULL;
  }
  
  if(prof->energy!=NULL){
    free(prof->energy);
    prof->energy=NULL;
  }
  
  if(prof->time!=NULL){
    free(prof->time);
    prof->time=NULL;
  }
  
  prof->NE=0;
  prof->Nt=0;
}

void destroyTESProfiles(TESProfiles* prof){
  
  int ii;
  if (prof!=NULL){
    if(prof->version!=NULL){
      for(ii=0; ii<prof->Nv; ii++){
	if(prof->version[ii]!=NULL){
	  free(prof->version[ii]);
	  prof->version[ii]=NULL;
	}
      }
      free(prof->version);
      prof->version=NULL;
    }
    if(prof->profiles!=NULL){
      for(ii=0; ii<prof->Nv; ii++){
	destroyTESProfilesEntries(&(prof->profiles[ii]));
      }
      free(prof->profiles);
      prof->profiles=NULL;
    }
    prof->Nv=0;
  }
}

void newTESProfilesEntries(TESProfilesEntries *prof){
  
  // Initialize values and pointers
  prof->Nt=0;
  prof->NE=0;
  prof->time=NULL;
  prof->energy=NULL;
}

TESProfiles* newTESProfiles(int* const status){
  
  TESProfiles* prof=(TESProfiles*)malloc(sizeof(TESProfiles));
  if (NULL==prof) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TESProfiles failed");
    return(prof);
  }
  
  // Initialize values and pointers
  prof->Nv=0;
  prof->version=NULL;
  prof->profiles=NULL;
  
  return prof;
}

int testTESProfilesExist(TESProfiles *prof, char *version)
{
  if(prof->Nv==0){
    return 0;
  }
  int ii;
  for(ii=0; ii<prof->Nv; ii++){
    if(!strcmp(version, prof->version[ii])){
      return 1;
    }
  }
  return 0;
}

void readTESProfiles(char *filename, 
			     char *version, 
			     TESProfiles *prof, 
			     int* const status)
{
  
  // Open fits file
  fitsfile *fptr=NULL;  
  if (fits_open_table(&fptr, filename, READONLY, status)){
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "File '%s' can not be opened.", filename);
    SIXT_ERROR(msg);
    CHECK_STATUS_VOID(*status);
  }
  
  // Move to version extension
  if (fits_movnam_hdu(fptr, BINARY_TBL, version, 0, status)){
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Extension %s in file '%s' can not be opened.", 
	   version, filename);
    SIXT_ERROR(msg);
    CHECK_STATUS_VOID(*status);
  }
  
  // Determine number of energy steps
  int ncols=0, nenergy;
  fits_get_num_cols(fptr, &ncols, status);
  CHECK_STATUS_VOID(*status);
  nenergy=ncols-1;
  
  // Determine number of time steps
  long ntime=0;
  fits_get_num_rows(fptr, &ntime, status);
  CHECK_STATUS_VOID(*status);
  
  // Allocate memory for the new version
  prof->profiles=(TESProfilesEntries*)realloc(prof->profiles, 
		    (1+prof->Nv)*(sizeof(TESProfilesEntries)));
  if (NULL==prof->profiles){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for CalorimeterProfiles failed");
    CHECK_STATUS_VOID(*status);
  }
  newTESProfilesEntries(&(prof->profiles[prof->Nv]));
  CHECK_STATUS_VOID(*status);
  
  // Alocate memory for the arrays
  prof->profiles[prof->Nv].Nt=ntime;
  prof->profiles[prof->Nv].NE=nenergy;
  
  prof->profiles[prof->Nv].time=(double*)malloc(ntime*sizeof(double));
  if(prof->profiles[prof->Nv].time==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for time array failed");
    CHECK_STATUS_VOID(*status);
  }
  
  prof->profiles[prof->Nv].energy=(double*)malloc(nenergy*sizeof(double));
  if(prof->profiles[prof->Nv].energy==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for energy array failed");
    CHECK_STATUS_VOID(*status);
  }
  
  prof->profiles[prof->Nv].adc_value=(double**)malloc(nenergy*sizeof(double*));
  if(prof->profiles[prof->Nv].adc_value==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for adc_value array failed");
    CHECK_STATUS_VOID(*status);
  }
  int ii;
  for(ii=0; ii<nenergy; ii++){
    prof->profiles[prof->Nv].adc_value[ii]=NULL;
  }
  for(ii=0; ii<nenergy; ii++){
    prof->profiles[prof->Nv].adc_value[ii]=(double*)malloc(ntime*sizeof(double));
    if(prof->profiles[prof->Nv].adc_value[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for adc_value array failed");
      CHECK_STATUS_VOID(*status);
    }
  }
  
  // Read FITS columns
  
  // Check if first column is time column
  char keyname[9], colname[9], comment[MAXMSG];
  sprintf(keyname, "TTYPE1");
  fits_read_key(fptr, TSTRING, keyname, colname, comment, status);
  CHECK_STATUS_VOID(*status);
  if(strcmp(colname, "TIME")){
    *status=EXIT_FAILURE;
    SIXT_ERROR("First column in input file is not a time column.");
    CHECK_STATUS_VOID(*status);
  }
  
  // Read time column
  int anynul=0;
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, ntime, prof->profiles[prof->Nv].time, 
		prof->profiles[prof->Nv].time, &anynul, status);
  CHECK_STATUS_VOID(*status);
  
  // Read the pulse columns
  for(ii=0; ii<nenergy; ii++){
    // Read the column name and determine the energy value
    sprintf(keyname, "TTYPE%d", ii+2);
    fits_read_key(fptr, TSTRING, keyname, colname, comment, status);
    CHECK_STATUS_VOID(*status);
    char *ptr=NULL;
    if(colname[0]!='E'){
      *status=EXIT_FAILURE;
      SIXT_ERROR("Column in input file is not a energy column.");
      CHECK_STATUS_VOID(*status);
    }
    ptr=&(colname[1]);
    prof->profiles[prof->Nv].energy[ii]=((double)atoi(ptr))/1000.; // Convert from eV to keV
    printf("Found a pulse template for E=%.3lfkeV.\n", prof->profiles[prof->Nv].energy[ii]);
    
    // Read pulse profile
    fits_read_col(fptr, TDOUBLE, 2+ii, 1, 1, ntime, prof->profiles[prof->Nv].adc_value[ii], 
		prof->profiles[prof->Nv].adc_value[ii], &anynul, status);
    CHECK_STATUS_VOID(*status);
  }
  
  //allocate memory for version array and print version string into last entry
  prof->version=(char**)realloc(prof->version, (1+prof->Nv)*sizeof(char*));
  if(prof->version==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for version array failed");
    CHECK_STATUS_VOID(*status);
  }
  prof->version[prof->Nv]=(char*)malloc((1+strlen(version))*sizeof(char));
  if(prof->version[prof->Nv]==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for version array failed");
    CHECK_STATUS_VOID(*status);
  }
  sprintf(prof->version[prof->Nv], "%s", version);  
  
  prof->Nv++;
  
  fits_close_file(fptr, status);
  CHECK_STATUS_VOID(*status);
  
}

int findTESProfileVersionIndex(TESProfiles* prof,
			       char *version)
{
  
  // Look through the version array and return the index
  // If the string matches
  int ii;
  for(ii=0; ii<prof->Nv; ii++){
    if(!strcmp(version, prof->version[ii])){
      return ii;
    }
  }
  // If the loop was left without finding the version,
  // a -1 is returned.
  return -1;
}

int findTESProfileEnergyIndex(TESProfiles* prof, 
			      int version, 
			      double energy)
{
  
  // Look through the energy array and find the index
  // for which the energy of the profile belongs to the
  // largest energy smaller than the photon energy
  int ii, best_match=prof->profiles[version].NE+1;
  
  for(ii=0; ii<prof->profiles[version].NE; ii++){
    if(prof->profiles[version].energy[ii]<energy){
      best_match=ii-1;
      break;
    }
  }
  // If the best match lies at the borders of the array,
  // make sure that it takes the border and not outside
  if(best_match<0){
    best_match=0;
  }else if(best_match>prof->profiles[version].NE-1){
    best_match=prof->profiles[version].NE-1;
  }
  return (best_match);
}

