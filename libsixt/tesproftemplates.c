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


   Copyright 2014 Jelle de Plaa & Thorsten Brand, FAU
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

int genTESProfile(TESTemplateInput* pinp, TESProfiles** ptemp, int* const status) {
    
  int i,j,k;
  double norm;
  
  
  if ((*ptemp)==NULL) {
    (*ptemp)=newTESProfiles(status);
  }
  
  /** For now, we create one version for all pixels */
  (*ptemp)->Nv=1;
  (*ptemp)->version=(char**)malloc((*ptemp)->Nv*sizeof(char*));
  if((*ptemp)->version==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("Malloc error.");
    return *status;
  }
  
  /* Allocate memory for different versions */
  (*ptemp)->profiles = (TESProfilesEntries*) malloc((*ptemp)->Nv * sizeof(TESProfilesEntries));
  if((*ptemp)->profiles==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("Malloc error.");
    return *status;
  }
  
  for (i=0;i<(*ptemp)->Nv;i++) {
    // copy version name
    (*ptemp)->version[i]=(char*)malloc((strlen(pinp->version[i])+1)*sizeof(char));
    if((*ptemp)->version[i]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("Malloc error.");
      return *status;
    }
    strcpy((*ptemp)->version[i], pinp->version[i]);
  
    /* Initialize profile pointer in entry */
    newTESProfilesEntries(&(*ptemp)->profiles[i]);
    
    /* Set dimensions of profiles arrays */
    (*ptemp)->profiles[i].Nt=pinp->nsamp;
    (*ptemp)->profiles[i].NE=pinp->ne;
    
    /* Allocate memory for energy and time arrays */
    (*ptemp)->profiles[i].energy = (double*)malloc((*ptemp)->profiles[i].NE * sizeof(double));
    if((*ptemp)->profiles[i].energy==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("Malloc error.");
      return *status;
    }
    (*ptemp)->profiles[i].time = (double*)malloc((*ptemp)->profiles[i].Nt * sizeof(double));
    if((*ptemp)->profiles[i].time==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("Malloc error.");
      return *status;
    }
    
    /* Fill energy arrays with values */
    for (j=0;j<(*ptemp)->profiles[i].NE;j++) {
      (*ptemp)->profiles[i].energy[j]=pinp->energies[j];
    } 
    
    /* Fill time arrays with values */
    for (j=0;j<(*ptemp)->profiles[i].Nt;j++) {
      (*ptemp)->profiles[i].time[j]=(double) j / pinp->freq;
    }
    
    /* Allocate adc_value arrays */
    (*ptemp)->profiles[i].adc_value= (double**) malloc((*ptemp)->profiles[i].NE*sizeof(double*));
    if((*ptemp)->profiles[i].adc_value==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("Malloc error.");
      return *status;
    }
    for (j=0;j<(*ptemp)->profiles[i].NE;j++) {
      (*ptemp)->profiles[i].adc_value[j]=NULL;
      (*ptemp)->profiles[i].adc_value[j]=(double*)malloc((*ptemp)->profiles[i].Nt*sizeof(double));
      if((*ptemp)->profiles[i].adc_value[j]==NULL){
	*status=EXIT_FAILURE;
	SIXT_ERROR("Malloc error.");
      return *status;
      }
    }
    
    /* Calculate profiles */
    for (j=0;j<(*ptemp)->profiles[i].NE;j++) {
      /* Normalize based on area (in s) calculated by analytic integration of pulse profile */
      norm = pinp->tfall - (pinp->tfall * pinp->trise)/(pinp->tfall + pinp->trise); 

      /* Calculate exponential pulse */
      for (k=0;k<(*ptemp)->profiles[i].Nt;k++) {
	(*ptemp)->profiles[i].adc_value[j][k]=ExponentialPulse(&(*ptemp)->profiles[i].time[k],&pinp->trise,&pinp->tfall);
	(*ptemp)->profiles[i].adc_value[j][k]=(*ptemp)->profiles[i].adc_value[j][k]/norm;
      }
	
    }
       
  }
  return *status;
} 


double ExponentialPulse(double *t, double *trise, double *tfall) {

  double result;
     
  result = (1.-exp(-*t/ *trise))*exp(-*t / *tfall);

  return result;
}

int createTESProfilesFile(char *filename, 
			   const char clobber,
			   char history,
			   int* const status)
{
  
  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, *status);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(filename);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", filename);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      CHECK_STATUS_RET(*status, *status);
    }
  }
  
  // Create the new empty file
  fitsfile *fptr=NULL;
  
  fits_create_file(&fptr, filename, status);
  CHECK_STATUS_RET(*status, *status);
  
  // Write the neccessary keywords to Primary extension
  int logic=(int)'T';
  int bitpix=8;
  int naxis=0;
  fits_update_key(fptr, TLOGICAL, "SIMPLE", &(logic), NULL, status);
  fits_update_key(fptr, TINT, "BITPIX", &(bitpix), NULL, status);
  fits_update_key(fptr, TINT, "NAXIS", &(naxis), NULL, status);
  CHECK_STATUS_RET(*status, *status);
  
  // Write history to primary header
  if (history) {
    HDpar_stamp(fptr, 1, status);
    CHECK_STATUS_RET(*status,*status);
  }
  
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, *status);

  return *status;
  
}

int writeTESProfiles(char *filename, 
			     char *version, 
			     TESProfilesEntries *prof, 
			     const char clobber,
			     char history,
			     int* const status)
{
  
  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  CHECK_STATUS_RET(*status, *status);
  if (0==exists) {
    createTESProfilesFile(filename, clobber, history, status);
    CHECK_STATUS_RET(*status, *status);
  }
  
  fitsfile *fptr=NULL;
  
  fits_open_file(&fptr, filename, READWRITE, status);
  CHECK_STATUS_RET(*status, *status);
  
  // Try to move to the version-extension
  if(fits_movnam_hdu(fptr, BINARY_TBL, version, 0, status)){
    // That version does not exist in the file yet, create it
    *status=EXIT_SUCCESS;
    char *name[]={"TIME"};
    char *type[]={"D"};
    char *unit[]={"s"};
    fits_create_tbl(fptr, BINARY_TBL, prof->Nt, 1, name,type,unit, version, status);
    CHECK_STATUS_RET(*status, *status);
    
    // Write the TIME column
    fits_write_col(fptr, TDOUBLE, 1, 1, 1, prof->Nt, 
		   prof->time, status);
    CHECK_STATUS_RET(*status, *status);
  }else{
    // Extension is already there, check if TIME column is compatible
    long nrows=0;
    int tcol=0;
    
    // First read number of rows and compare to length of time array
    fits_get_num_rows(fptr, &nrows, status);
    CHECK_STATUS_RET(*status, *status);
    
    if(nrows!=prof->Nt){
      // Length does not match to already existing table
      *status=EXIT_FAILURE;
      SIXT_ERROR("The specified file already contains the same pulse template version but with different length.");
      CHECK_STATUS_RET(*status, *status);
    }
    
    // Read the TIME column
    fits_get_colnum(fptr, CASEINSEN, "TIME", &tcol, status);
    if(status!=EXIT_SUCCESS){
      // TIME column does not exist
      SIXT_ERROR("The specified file already contains the same pulse template version but without a TIME column.");
      CHECK_STATUS_RET(*status, *status);
    }
    
    double nulval=0.;
    double tarray[nrows];
    int anynul=0;
    
    fits_read_col(fptr, TDOUBLE, tcol, 1, 1, nrows, 
		  &nulval, tarray, &anynul, status);
    CHECK_STATUS_RET(*status, *status);
    
    long ii;
    // loop over all time steps and compare time values
    for(ii=0; ii<nrows; ii++){
      if(tarray[ii]!=prof->time[ii]){
	*status=EXIT_FAILURE;
	SIXT_ERROR("The specified file already contains the same pulse template version but with different TIME vector.");
	CHECK_STATUS_RET(*status, *status);
      }
    }
  }
    
  // Now fptr should point to a valid template file
  // Insert the adc columns one by one in correct position
  
  int jj;
  
  for(jj=0; jj<prof->NE; jj++){
    InsertTESProfADCCol(fptr, prof->Nt, prof->energy[jj], 
			prof->adc_value[jj], status);
    CHECK_STATUS_RET(*status, *status);
  }
  
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, *status);
  
  return *status;
  
}

int InsertTESProfADCCol(fitsfile *fptr, 
			long nt, 
			double energy, 
			double *adc, 
			int* const status)
{
 
  // Construct new column name
  char ename[9];
  sprintf(ename, "E%07.0lf", energy*1000.);
  
  // Read number of existing columns
  int ncols;
  fits_get_num_cols(fptr, &ncols, status);
  CHECK_STATUS_RET(*status, *status);
  
  int col=2;
  
  // Determine place for new column
  if(ncols==1){
    // Only TIME column is present until now, place ADC column behind
    col=2;
  }else{
    // There are already ADC columns, look for correct place
    int ii;
    for(ii=0; ii<ncols-1; ii++){
      // Read column number
      char colname[9], keyname[9], comment[MAXMSG];
      sprintf(keyname, "TTYPE%d", ii+2);
      fits_read_key(fptr, TSTRING, keyname, colname, comment, status);
      CHECK_STATUS_RET(*status, *status);
      
      // Determine energy of ADC column
      char *ptr=NULL;
      double colen;
      if(colname[0]!='E'){
	*status=EXIT_FAILURE;
	SIXT_ERROR("Column in input file is not a energy column.");
	CHECK_STATUS_RET(*status, *status);
      }
      ptr=&(colname[1]);
      colen=((double)atoi(ptr))/1000.;
      if(colen<energy){
	col++;
      }else if(colen==energy){
	*status=EXIT_FAILURE;
	SIXT_ERROR("Input file already contains a template for identical energy.");
	CHECK_STATUS_RET(*status, *status);
      }else{
	break;
      }
    }
  }
  
  // Now insert a new column
  fits_insert_col(fptr, col, ename, "D", status);
  CHECK_STATUS_RET(*status, *status);
  
  // Write the adc-array to the new column
  fits_write_col(fptr, TDOUBLE, col, 1, 1, nt, adc, status);
  CHECK_STATUS_RET(*status, *status);
  
  return *status;
}
