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

#include "teseventfile.h"

/** Constructor. Returns a pointer to an empty TesEventFile data
    structure. */
TesEventFile* newTesEventFile(int* const status) {
  TesEventFile* file=(TesEventFile*)malloc(sizeof(TesEventFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesEventFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr    =NULL;

  // Initialize values.
  file->nrows	=0;
  file->row  	=0;
  file->timeCol =1;
  file->evtCol  =2;
  file->pixID   =0;

  return(file);
}

/** Destructor. */
void freeTesEventFile(TesEventFile** const file, int* const status){
  if (NULL!=*file) {
    (*file)->nrows	      =0;
    (*file)->row  	      =0;
    (*file)->pixID	      =0;
    (*file)->timeCol          =0;
    (*file)->evtCol           =0;
    
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed TesStream list file (containing %ld rows).\n", 
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}

/** Create and open a new TesEventFile. */
TesEventFile* opennewTesEventFile(const char* const filename,
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
				  const int pixID,
				  int eventSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status){

  TesEventFile* file = newTesEventFile(status);
  CHECK_STATUS_RET(*status, file);

  int exists;
  char buffer[MAXFILENAME];
  sprintf(buffer, "%s_%d.fits", filename,pixID);
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
  fits_update_key(file->fptr, TINT, "EVENTSZ", &(eventSize), "Number of samples in events", status);
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
  
  char *ttype[2];
  char *tform[2];
  char *tunit[2];
  
  int ii;
  
  for(ii=0; ii<2; ii++){
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
  
  //Create second column PXL00001
  sprintf(ttype[1], "PXL%05d", pixID);
  sprintf(tform[1], "%dU",eventSize);
  sprintf(tunit[1], "ADC");
  CHECK_STATUS_RET(*status,file);
  
  char extName[9];
  sprintf(extName,"ADC%03d",pixID);
  fits_create_tbl(file->fptr, BINARY_TBL, 0, 2,
		  ttype, tform, tunit,extName, status);
  //Add keywords to other extension
  fits_update_key(file->fptr, TINT, "EVENTSZ", &(eventSize), "Number of samples in events", status);
  fits_update_key(file->fptr, TINT, "PREBUFF", &(preBufferSize), "Number of samples before start of pulse", status);
  fits_update_key(file->fptr, TDOUBLE, "DELTAT", &(deltat), "Time resolution of data stream", status);
  CHECK_STATUS_RET(*status,file);

  return(file);

}


void writeEvents2FITS(TesEventFile** outputFiles,TESDataStream* stream,PixImpFile* impfile,
		      int pixlow,int Npix,double tstart,double tstartTES,
		      double tstop,double sampleFreq,int eventSize,int preBufferSize,
		      float monoen,int* const status){
  
  //Allocation of array containing the events being constructed
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
  // Number of triggered events per pixel
  int* numberTrigger = (int*)malloc(Npix*sizeof(int));
  // Number of simulated events per pixel
  int* numberSimulated = (int*)malloc(Npix*sizeof(int));
  // Pixel iterator
  int pixNumber=0;
  // Current time stamp
  double time;
  //Array to save at which point we are while constructiong an event
  int* positionInEvent = (int*)malloc(Npix*sizeof(int)); 

  //Initializations
  for(ii=0; ii<Npix; ii++){
    adc_value[ii]=(uint16_t*)malloc(eventSize*sizeof(uint16_t));
    positionInEvent[ii]=-1;
    numberTrigger[ii]=0;
    numberSimulated[ii]=0;
    if(adc_value[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for adc_value failed.");
      CHECK_STATUS_VOID(*status);
    }
    for (jj=0;jj<eventSize;jj++){
      adc_value[ii][jj]=0;
    }
  }

  
 

  /*Old pulse counter -> not relliable
  //Array to save start times of pulses that may be recorded and not counted otherwise
  double** save_times = (double**)malloc(Npix*sizeof(double*));//if more than 10 events in a pixel are to be saved that way, this means the CR is probably too high
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
    /* Get first event in correct time frame from the impact file */
    if (tstep==0) {
      do {
	piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
	CHECK_STATUS_VOID(*status);
      } while(impact.time<(tstart+(double)preBufferSize/sampleFreq));
    }
   
    /* If the event occurs in a time bin away from preBufferSize samples of the current time, is in an active pixel, */ 
    /* and we are not recording, put flag to record*/
    while ((piximpstatus>0) &&(impact.time>=t+(double)preBufferSize/sampleFreq)&&(impact.time<t+(((double)preBufferSize+1.0)/sampleFreq))) {
      if((impact.pixID>=pixlow) && (impact.pixID<=(pixlow+Npix))){
	if (positionInEvent[impact.pixID-pixlow]==-1) {
	  positionInEvent[impact.pixID-pixlow]=0;
	  //Check if pulses are gonna be in the pre-buffer
	  //for (ii=0;ii<nb_saved_times[impact.pixID-pixlow];ii++){
	  //  if (t<save_times[impact.pixID-pixlow][ii]){
	  //    numberTrigger[impact.pixID-pixlow]++;
	  // }
	  //}
	}
      //if (impact.time<(t+(double)(eventSize-positionInEvent[impact.pixID-pixlow])/sampleFreq)){ //if the pulse if gonna be in the event currently being recorded
      //	  numberTrigger[impact.pixID-pixlow]++;
      //	}
      //	else { //save the start time in case another event triggers and this pulse is in the preBuffer part
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
    
    // Record the ADC values in the events if it is necessary
    
    for (pixNumber=0;pixNumber<Npix;pixNumber++) {
      if ((positionInEvent[pixNumber]>=0) && (positionInEvent[pixNumber]<eventSize)) {
	adc_value[pixNumber][positionInEvent[pixNumber]] = stream->adc_value[(int)((t-tstartTES)*sampleFreq)][pixNumber];
	positionInEvent[pixNumber]++;
      }
      //If Event is full, record in output FITS file
      if (positionInEvent[pixNumber]==eventSize){
	outputFiles[pixNumber]->nrows++;
	outputFiles[pixNumber]->row++;
	time = t-(eventSize-1)/sampleFreq;
	fits_write_col(outputFiles[pixNumber]->fptr, TDOUBLE, outputFiles[pixNumber]->timeCol, 
		       outputFiles[pixNumber]->row, 1, 1, &time, status);
	fits_write_col(outputFiles[pixNumber]->fptr, TUSHORT, outputFiles[pixNumber]->evtCol, 
		       outputFiles[pixNumber]->row, 1, eventSize, (adc_value[pixNumber]), status);
	//Reset adc_value array and positionInEvent
	for (jj=0;jj<eventSize;jj++){
	  adc_value[pixNumber][jj]=0;
	}
	positionInEvent[pixNumber]=-1;
	CHECK_STATUS_VOID(*status);
      }
    }
    fflush(stdout);
    t=t+1./sampleFreq;
  }

  //Count number of pulses that should be visible
  double** outputTimeCol = (double**)malloc(Npix*sizeof(double*));
  int* currentTimeIndex = (int*)malloc(Npix*sizeof(int));
  int anynul=0;
  for (pixNumber=0;pixNumber<Npix;pixNumber++) {
      outputTimeCol[pixNumber] = (double*)malloc(outputFiles[pixNumber]->nrows * sizeof(double));
      currentTimeIndex[pixNumber] = 0;
      //Read time column in output file
      fits_read_col(outputFiles[pixNumber]->fptr, TDOUBLE, outputFiles[pixNumber]->timeCol, 
		    1,1,outputFiles[pixNumber]->nrows, NULL, (outputTimeCol[pixNumber]), &anynul,status);
      CHECK_STATUS_VOID(*status);
  }
  //Reinitialize impfile
  impfile->row=0;
  //Get first impact
  piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
  CHECK_STATUS_VOID(*status);

  //Iterate over the impacts
  while ((piximpstatus>0) && (impact.time<tstop)){
    if ((impact.pixID>=pixlow) && (impact.pixID<=(pixlow+Npix))){
      numberSimulated[impact.pixID-pixlow]++;
      while ((impact.time>(outputTimeCol[impact.pixID-pixlow][currentTimeIndex[impact.pixID-pixlow]]+(double)eventSize/sampleFreq)) &&
	     (currentTimeIndex[impact.pixID-pixlow]<outputFiles[impact.pixID-pixlow]->nrows-1 )){
	       currentTimeIndex[impact.pixID-pixlow]++;
      }
      if (currentTimeIndex[impact.pixID-pixlow]<outputFiles[impact.pixID-pixlow]->nrows){//If we have reach the end of the impacts 
	if (impact.time>(outputTimeCol[impact.pixID-pixlow][currentTimeIndex[impact.pixID-pixlow]])) {
	  numberTrigger[impact.pixID-pixlow]++;
	}
      }
    }
    piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
    CHECK_STATUS_VOID(*status);
  } 
  //Save keywords with number of counts and monochromatic energy
  char keyword[9];
  for (pixNumber=0;pixNumber<Npix;pixNumber++) {
    sprintf(keyword,"NES%05d",pixNumber+1);
    fits_update_key(outputFiles[pixNumber]->fptr, TINT, keyword, &(numberSimulated[pixNumber]), "Number of simulated events", status);
    sprintf(keyword,"NET%05d",pixNumber+1);
    fits_update_key(outputFiles[pixNumber]->fptr, TINT, keyword, &(numberTrigger[pixNumber]), "Number of triggered events", status);
    fits_update_key(outputFiles[pixNumber]->fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);
  }

  //Free memory
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
  free(positionInEvent);
  positionInEvent=NULL;

}


