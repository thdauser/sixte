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


   Copyright 2014 Philippe Peille, IRAP
*/

#include "testriggerfile.h"

/** Constructor. Returns a pointer to an empty TesTriggerFile data
    structure. */
TesTriggerFile* newTesTriggerFile(int* const status) {
  TesTriggerFile* file=(TesTriggerFile*)malloc(sizeof(TesTriggerFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesTriggerFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr    =NULL;

  // Initialize values.
  file->nrows	   =0;
  file->row  	   =0;
  file->timeCol    =1;
  file->trigCol    =2;
  file->pixIDCol   =3;
  file->ph_idCol   =4;

  return(file);
}

/** Destructor. */
void freeTesTriggerFile(TesTriggerFile** const file, int* const status){
  if (NULL!=*file) {
    (*file)->nrows	      =0;
    (*file)->row  	      =0;
    (*file)->timeCol          =0;
    (*file)->trigCol          =0;
    (*file)->ph_idCol         =0;
    (*file)->pixIDCol         =0;
    
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed TesStream list file (containing %ld rows).\n", 
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}

/** Create and open a new TesTriggerFile. */
TesTriggerFile* opennewTesTriggerFile(const char* const filename,
				  char* const telescop,
				  char* const instrume,
				  char* const filter,
				  char* const ancrfile,
				  char* const respfile,
				  char* const xmlfile,
				  char* const impactlist,
				  const double mjdref,
				  const double timezero,
				  const double tstart,
				  const double tstop,
				  int triggerSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status){

  TesTriggerFile* file = newTesTriggerFile(status);
  CHECK_STATUS_RET(*status, file);

  int exists;
  char buffer[MAXFILENAME];
  sprintf(buffer,"%s",filename);
  fits_file_exists(buffer, &exists, status);
  CHECK_STATUS_RET(*status,file);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(buffer);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", buffer);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      CHECK_STATUS_RET(*status,file);
    }
  }
  fits_create_file(&file->fptr,buffer, status);
  CHECK_STATUS_RET(*status,file);
  int logic=(int)'T';
  int bitpix=8;
  int naxis=0;
  fits_update_key(file->fptr, TLOGICAL, "SIMPLE", &(logic), NULL, status);
  fits_update_key(file->fptr, TINT, "BITPIX", &(bitpix), NULL, status);
  fits_update_key(file->fptr, TINT, "NAXIS", &(naxis), NULL, status);
  fits_update_key(file->fptr, TINT, "TRIGGSZ", &(triggerSize), "Number of samples in triggers", status);
  fits_update_key(file->fptr, TINT, "PREBUFF", &(preBufferSize), "Number of samples before start of pulse", status);
  double deltat = 1./sampleFreq;
  fits_update_key(file->fptr, TDOUBLE, "DELTAT", &(deltat), "Time resolution of data stream", status);
  sixt_add_fits_stdkeywords(file->fptr, 1, telescop, instrume, filter,
			    ancrfile, respfile, mjdref, timezero, 
			    tstart, tstop, status);
  CHECK_STATUS_RET(*status,file);
  
  //Write XML into header
  char comment[MAXMSG];
  sprintf(comment, "XMLFILE: %s", xmlfile);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status,file);
  
  //Write pixel impact file into header
  sprintf(comment, "PIXFILE: %s", impactlist);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status,file);

   // Create table
  int tlen=9;
  
  char *ttype[4];
  char *tform[4];
  char *tunit[4];
  
  int ii;
  
  for(ii=0; ii<4; ii++){
    ttype[ii]=(char*)malloc(tlen*sizeof(char));
    if(ttype[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for ttype failed");
      CHECK_STATUS_RET(*status,file);
    }
    tform[ii]=(char*)malloc(tlen*sizeof(char));
    if(tform[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for tform failed");
      CHECK_STATUS_RET(*status,file);
    }
    tunit[ii]=(char*)malloc(tlen*sizeof(char));
    if(tunit[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for tunit failed");
      CHECK_STATUS_RET(*status,file);
    }
  }
  CHECK_STATUS_RET(*status,file);
  
  //Create first column TIME
  sprintf(ttype[0], "TIME");
  sprintf(tform[0], "1D");
  sprintf(tunit[0], "s");
  
  //Create second column ADC
  sprintf(ttype[1], "ADC");
  sprintf(tform[1], "%dU",triggerSize);
  sprintf(tunit[1], "ADC");

  //Create third column PIXID
  sprintf(ttype[2], "PIXID");
  sprintf(tform[2], "1J");
  sprintf(tunit[2], "");
  CHECK_STATUS_RET(*status,file);

  //Create fourth column PH_ID
  sprintf(ttype[3], "PH_ID");
  sprintf(tform[3], "1PJ(%d)",MAXIMPACTNUMBER);
  sprintf(tunit[3], "");
  CHECK_STATUS_RET(*status,file);

  char extName[9];
  sprintf(extName,"RECORDS");
  fits_create_tbl(file->fptr, BINARY_TBL, 0, 4,
		  ttype, tform, tunit,extName, status);
  //Add keywords to other extension
  fits_update_key(file->fptr, TINT, "TRIGGSZ", &(triggerSize), "Number of samples in triggers", status);
  fits_update_key(file->fptr, TINT, "PREBUFF", &(preBufferSize), "Number of samples before start of pulse", status);
  fits_update_key(file->fptr, TDOUBLE, "DELTAT", &(deltat), "Time resolution of data stream", status);
  CHECK_STATUS_RET(*status,file);

  return(file);

}


void writeTriggerFileWithImpact(TESDataStream* const stream,
				char* const tesTriggerFilename,char* const telescop,
				char* const instrume,char* const filter,
				char* const ancrfile,char* const respfile,
				char* const xmlfile,char* const impactlist,
				const double mjdref,const double timezero,
				double tstart,double tstop,const int triggerSize,
				const int preBufferSize,const double sampleFreq,
				const char clobber,const int pixlow,const int Npix,
				float monoen,int* const status){

  ////////////////////////////////
  //Open output file
  ////////////////////////////////
  TesTriggerFile* outputFile = opennewTesTriggerFile(tesTriggerFilename,
						     telescop,
						     instrume,
						     filter,
						     ancrfile,
						     respfile,
						     xmlfile,
						     impactlist,
						     mjdref,
						     timezero,
						     tstart,
						     tstop,
						     triggerSize,
						     preBufferSize,
						     sampleFreq,
						     clobber,
						     status);
  CHECK_STATUS_VOID(*status);
    
  ////////////////////////////////
  //Open Impact file
  ////////////////////////////////
  PixImpFile* impfile=openPixImpFile(impactlist, READONLY, status);
  CHECK_STATUS_VOID(*status);
  double tstartImp=0;
  double tstopImp=0;
  char comment[MAXMSG];
  fits_read_key(impfile->fptr, TDOUBLE, "TSTART", &tstartImp, comment, status);
  fits_read_key(impfile->fptr, TDOUBLE, "TSTOP", &tstopImp, comment, status);
  CHECK_STATUS_VOID(*status);

  //Check if tstart/tstop are compatible with pix impact file and correct if necessary
  double tstartTES = tstart;
  printf("Pix impact file reaches from %lfs-%lfs .\n", tstartImp, tstopImp);
  if(tstartImp>tstart){
    if(tstartImp>tstop){
      SIXT_ERROR("Impact file tstart is larger than end of TES ADC data -> abort");
      *status=EXIT_FAILURE;
      CHECK_STATUS_VOID(*status);
    }
    puts("Impact file tstart is larger than in TES ADC data.");    
    tstart=tstartImp;
  }
  if(tstopImp<tstop){
    if(tstopImp<tstart){
      SIXT_ERROR("Impact file tstop is smaller than start of TES ADC data -> abort");
      *status=EXIT_FAILURE;
      CHECK_STATUS_VOID(*status);
    }
    puts("Impact file tstop is smaller than in TES ADC data.");
    tstop=tstopImp;
  }
  printf("Simulate from %lfs-%lfs .\n", tstart, tstop);



  //Allocation of array containing the triggers being constructed
  uint16_t** adc_value = (uint16_t**)malloc(Npix*sizeof(uint16_t*));
  if(adc_value==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for adc_value failed.");
    CHECK_STATUS_VOID(*status);
  }
  int ii;
  int jj;
  //Current time
  double t = tstart;
  //Current impact
  PixImpact impact;
  // Status of getNextImpactFromPixImpFile
  int piximpstatus=0;
  // Number of triggers per pixel
  int* numberTrigger = (int*)malloc(Npix*sizeof(int));
  // Number of simulated pulses per pixel
  int* numberSimulated = (int*)malloc(Npix*sizeof(int));
  // Pixel iterator
  int pixNumber=0;
  // Current time stamp
  double time;
  // Current pixID
  long pixID;
  //Array to save at which point we are while constructiong a trigger
  int* positionInTrigger = (int*)malloc(Npix*sizeof(int)); 

  //Initializations
  for(ii=0; ii<Npix; ii++){
    adc_value[ii]=(uint16_t*)malloc(triggerSize*sizeof(uint16_t));
    positionInTrigger[ii]=-1;
    numberTrigger[ii]=0;
    numberSimulated[ii]=0;
    if(adc_value[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for adc_value failed.");
      CHECK_STATUS_VOID(*status);
    }
    for (jj=0;jj<triggerSize;jj++){
      adc_value[ii][jj]=0;
    }
  }

  
 

  /*Old pulse counter -> not relliable
  //Array to save start times of pulses that may be recorded and not counted otherwise
  double** save_times = (double**)malloc(Npix*sizeof(double*));//if more than 10 triggers in a pixel are to be saved that way, this means the CR is probably too high
  int* nb_saved_times = (int*)malloc(Npix*sizeof(int));//number of start times in the array that are actually still worth looking at
  for (ii =0;ii<Npix;ii++){
    save_times[ii]=(double*)malloc(10*sizeof(double));
    nb_saved_times[ii]=0;
    for (int ibuffer=0;ibuffer<10;ibuffer++){
      save_times[ii][ibuffer]=0;
    }
    }*/

  //Loop over times
  int Nt = (int)((tstop-tstart)*sampleFreq-preBufferSize);
  for (int tstep=0;tstep < Nt;tstep++){
    /* Get first pulse in correct time frame from the impact file */
    if (tstep==0) {
      do {
	piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
	CHECK_STATUS_VOID(*status);
      } while(impact.time<(tstart+(double)preBufferSize/sampleFreq));
    }
   
    /* If the pulse occurs in a time bin away from preBufferSize samples of the current time, is in an active pixel, */ 
    /* and we are not recording, put flag to record*/
    while ((piximpstatus>0) &&(impact.time>=t+(double)preBufferSize/sampleFreq)&&(impact.time<t+(((double)preBufferSize+1.0)/sampleFreq))) {
      if((impact.pixID>=pixlow) && (impact.pixID<(pixlow+Npix))){
	if (positionInTrigger[impact.pixID-pixlow]==-1) {
	  positionInTrigger[impact.pixID-pixlow]=0;
	  //Check if pulses are gonna be in the pre-buffer
	  //for (ii=0;ii<nb_saved_times[impact.pixID-pixlow];ii++){
	  //  if (t<save_times[impact.pixID-pixlow][ii]){
	  //    numberTrigger[impact.pixID-pixlow]++;
	  // }
	  //}
	}
      //if (impact.time<(t+(double)(triggerSize-positionInTrigger[impact.pixID-pixlow])/sampleFreq)){ //if the pulse if gonna be in the trigger currently being recorded
      //	  numberTrigger[impact.pixID-pixlow]++;
      //	}
      //	else { //save the start time in case another pulse triggers and this pulse is in the preBuffer part
      //	  if (nb_saved_times[impact.pixID-pixlow]>10){
      //    *status=EXIT_FAILURE;
      //    SIXT_ERROR("Too many pulses happened in a pre-buffer. The count rate is probably too high.");
      //    CHECK_STATUS_VOID(*status);
      //  }
      //  save_times[impact.pixID-pixlow][nb_saved_times[impact.pixID-pixlow]] = impact.time;
      //  nb_saved_times[impact.pixID-pixlow]++;
      //}
      //numberSimulated[impact.pixID-pixlow]++;
      }
      CHECK_STATUS_VOID(*status);
      piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
      CHECK_STATUS_VOID(*status);
    }
    
    // Record the ADC values in the triggers if it is necessary
    
    for (pixNumber=0;pixNumber<Npix;pixNumber++) {
      if ((positionInTrigger[pixNumber]>=0) && (positionInTrigger[pixNumber]<triggerSize)) {
	adc_value[pixNumber][positionInTrigger[pixNumber]] = stream->adc_value[(int)round((t-tstartTES)*sampleFreq)][pixNumber];
	positionInTrigger[pixNumber]++;
      }
      //If Trigger is full, record in output FITS file
      if (positionInTrigger[pixNumber]==triggerSize){
	outputFile->nrows++;
	outputFile->row++;
	time = t-(triggerSize-1)/sampleFreq;
	pixID = pixNumber+pixlow+1;
	fits_write_col(outputFile->fptr, TDOUBLE, outputFile->timeCol, 
		       outputFile->row, 1, 1, &time, status);
	fits_write_col(outputFile->fptr, TUSHORT, outputFile->trigCol, 
		       outputFile->row, 1, triggerSize, (adc_value[pixNumber]), status);
	fits_write_col(outputFile->fptr, TLONG, outputFile->pixIDCol, 
		       outputFile->row, 1, 1, &pixID, status);
	//Reset adc_value array and positionInTrigger
	for (jj=0;jj<triggerSize;jj++){
	  adc_value[pixNumber][jj]=0;
	}
	positionInTrigger[pixNumber]=-1;
	CHECK_STATUS_VOID(*status);
      }
    }
    fflush(stdout);
    t=t+1./sampleFreq;
  }

  //Count number of pulses that should be visible
  double* outputTimeCol = (double*)malloc(outputFile->nrows*sizeof(double*)); //Time column of the trigger file
  long* outputPixIDCol = (long*)malloc(outputFile->nrows*sizeof(long*)); //PixID column of the trigger file
  int* currentTimeIndex = (int*)malloc(Npix*sizeof(int)); //Current index of time column number for each pixel
  long** currentImpactArray = (long**)malloc(Npix*sizeof(long*)); //Array containing for each pixel the PH_ID in the current trigger
  int* currentImpactNumber = (int*)malloc(Npix*sizeof(int*)); //Number of impacts in the current trigger for each pixel
  unsigned char* eofArray = (unsigned char*)malloc(Npix*sizeof(unsigned char*));//Array containing the EOF signal for each pixel
  int anynul=0;

  //Read time column in output file
  fits_read_col(outputFile->fptr, TDOUBLE, outputFile->timeCol, 
		1,1,outputFile->nrows, NULL, outputTimeCol, &anynul,status);
  CHECK_STATUS_VOID(*status);

  //Read pixID column in output file
  fits_read_col(outputFile->fptr, TLONG, outputFile->pixIDCol, 
		1,1,outputFile->nrows, NULL, outputPixIDCol, &anynul,status);
  CHECK_STATUS_VOID(*status);

  //Initialize arrays
  int timeColIterator = 0;//iterator over the time column
  for (pixNumber=0;pixNumber<Npix;pixNumber++) {
    //Look for first record of corresponding pixID
    timeColIterator=0;
    while (outputPixIDCol[timeColIterator]!=(pixNumber+pixlow+1)) {
      timeColIterator+=1;
      if (timeColIterator==outputFile->nrows){
	timeColIterator=-1;
	break;
      }
    }
    if (timeColIterator!=-1) {
      currentTimeIndex[pixNumber] = timeColIterator;
    }
    
    //Allocate impact array
    currentImpactArray[pixNumber] = (long*)malloc(MAXIMPACTNUMBER*sizeof(long));
    for (int jj=0;jj<MAXIMPACTNUMBER;jj++){
      currentImpactArray[pixNumber][jj]=0;
    }
    
    //Set impact number for each pixel
    currentImpactNumber[pixNumber] = 0;
    //Set eof indicator for each pixel
    eofArray[pixNumber]=1;
  }
  
  //Reinitialize impfile
  impfile->row=0;
  //Get first impact after tstart
  do {
	piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
	CHECK_STATUS_VOID(*status);
  } while(impact.time<tstart);

  //Reinitialize trigger file
  outputFile->row=1;

  //Iterate over the impacts
  while ((piximpstatus>0) && (impact.time<tstop)){
    if ((impact.pixID>=pixlow) && (impact.pixID<(pixlow+Npix))){
      numberSimulated[impact.pixID-pixlow]++;
      if ((impact.time>(outputTimeCol[currentTimeIndex[impact.pixID-pixlow]]+(double)triggerSize/sampleFreq)) &&
	  (eofArray[impact.pixID-pixlow])){
	//Change trigger -> save corresponding impacts and go to next
	fits_write_col(outputFile->fptr, TLONG, outputFile->ph_idCol, 
		       currentTimeIndex[impact.pixID-pixlow]+1, 1,currentImpactNumber[impact.pixID-pixlow],currentImpactArray[impact.pixID-pixlow], status);
	CHECK_STATUS_VOID(*status);
	//Look for next trigger of corresponding pixel
	if (currentTimeIndex[impact.pixID-pixlow]==outputFile->nrows-1) {
	  eofArray[impact.pixID-pixlow]=0;
	}
	else{
	  timeColIterator=currentTimeIndex[impact.pixID-pixlow]+1;
	  while (outputPixIDCol[timeColIterator]!=impact.pixID+1) {
	    timeColIterator+=1;
	    if (timeColIterator==outputFile->nrows){
	      eofArray[impact.pixID-pixlow]=0;
	      break;
	    }
	  }
	  if (eofArray[impact.pixID-pixlow]) {
	    currentTimeIndex[impact.pixID-pixlow]=timeColIterator;
	  }
	}
	currentImpactNumber[impact.pixID-pixlow]=0;
      }
      //If we have not reached the end of the impacts and are still inside the trigger
      if ((eofArray[impact.pixID-pixlow]) && 
	  (impact.time>(outputTimeCol[currentTimeIndex[impact.pixID-pixlow]]))){
	numberTrigger[impact.pixID-pixlow]++;
	currentImpactArray[impact.pixID-pixlow][currentImpactNumber[impact.pixID-pixlow]] = impact.ph_id;
	currentImpactNumber[impact.pixID-pixlow]++;
      }
    }
    piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
    CHECK_STATUS_VOID(*status);
  } 
  //Save keywords with first, last and N pixels
  int firstpix = pixlow+1;
  int lastpix = pixlow+Npix;
  int numberpix = Npix;
  fits_update_key(outputFile->fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in trigger file", status);
  fits_update_key(outputFile->fptr, TINT, "LASTPIX", &lastpix, "Last pixel in trigger file", status);
  fits_update_key(outputFile->fptr, TINT, "NPIX", &numberpix, "Number of pixels in trigger file", status);
  //Save keywords with number of counts
  char keyword[9];
  for (pixNumber=0;pixNumber<Npix;pixNumber++) {
    if (currentImpactNumber[pixNumber]>0){
      fits_write_col(outputFile->fptr, TLONG, outputFile->ph_idCol, 
		     currentTimeIndex[impact.pixID-pixlow]+1, 1,currentImpactNumber[pixNumber],currentImpactArray[pixNumber], status);
      CHECK_STATUS_VOID(*status);
    }
    sprintf(keyword,"NES%05d",pixNumber+pixlow+1);
    fits_update_key(outputFile->fptr, TINT, keyword, &(numberSimulated[pixNumber]), "Number of simulated pulses", status);
    sprintf(keyword,"NET%05d",pixNumber+pixlow+1);
    fits_update_key(outputFile->fptr, TINT, keyword, &(numberTrigger[pixNumber]), "Number of triggered pulses", status);
  }
  //Save monochromatic energy keyword
  fits_update_key(outputFile->fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);

  //Free memory
  freeTesTriggerFile(&(outputFile), status);
  freePixImpFile(&impfile, status);
  if(adc_value!=NULL){
    for(ii=0; ii<Npix; ii++){
      if(adc_value[ii]!=NULL){
	free(adc_value[ii]);
	adc_value[ii]=NULL;
      }
    }
    free(adc_value);
    adc_value=NULL;
  }
  free(positionInTrigger);
  positionInTrigger=NULL;
  free(outputTimeCol);
  outputTimeCol=NULL;
  free(outputPixIDCol);
  outputPixIDCol=NULL;
  free(currentTimeIndex);
  currentTimeIndex=NULL;
  if (currentImpactArray!=NULL){
    for(ii=0; ii<Npix; ii++){
      if(currentImpactArray[ii]!=NULL){
	free(currentImpactArray[ii]);
	currentImpactArray[ii]=NULL;
      }
    }
    free(currentImpactArray);
    currentImpactArray=NULL;
  }
  free(currentImpactNumber);
  currentImpactNumber=NULL;
  free(eofArray);
  eofArray=NULL;
  free(numberSimulated);
  numberSimulated=NULL;
  free(numberTrigger);
  numberTrigger=NULL;

}


