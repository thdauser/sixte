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
TesTriggerFile* newTesTriggerFile(int triggerSize,int* const status) {
  TesTriggerFile* file=(TesTriggerFile*)malloc(sizeof(TesTriggerFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesTriggerFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr    =NULL;

  // Initialize values.
  file->nrows	     =0;
  file->row  	     =1;
  file->timeCol      =1;
  file->trigCol      =2;
  file->pixIDCol     =3;
  file->ph_idCol     =4;
  file->trigger_size =triggerSize;
  file->delta_t=0;

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

  TesTriggerFile* file = newTesTriggerFile(triggerSize,status);
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

  //Free memory
  for(ii=0; ii<4; ii++){
    free(ttype[ii]);
    ttype[ii]=NULL;
    free(tform[ii]);
    tform[ii]=NULL;
    free(tunit[ii]);
    tunit[ii]=NULL;
  }

  return(file);

}

/** Save pixels, NES/NET and monoen keywords to the given FITS file */
void saveTriggerKeywords(fitsfile* fptr,int firstpix,int lastpix,int numberpix,float monoen,
		int* const numberSimulated,int* const numberTrigger,int* const status){
	fits_update_key(fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in trigger file", status);
	fits_update_key(fptr, TINT, "LASTPIX", &lastpix, "Last pixel in trigger file", status);
	fits_update_key(fptr, TINT, "NPIX", &numberpix, "Number of pixels in trigger file", status);
	fits_update_key(fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);
	CHECK_STATUS_VOID(*status);

	//Save keywords with number of counts
	char keyword[9];
	for (int pixNumber=0;pixNumber<numberpix;pixNumber++) {
		sprintf(keyword,"NES%05d",pixNumber+firstpix);
		fits_update_key(fptr, TINT, keyword, &(numberSimulated[pixNumber]), "Number of simulated pulses", status);
		sprintf(keyword,"NET%05d",pixNumber+firstpix);
		fits_update_key(fptr, TINT, keyword, &(numberTrigger[pixNumber]), "Number of triggered pulses", status);
	}
	CHECK_STATUS_VOID(*status);
}

/** Open an existing TesTriggerFile. */
TesTriggerFile* openexistingTesTriggerFile(const char* const filename,int* const status){
	TesTriggerFile* file = newTesTriggerFile(0,status);
	CHECK_STATUS_RET(*status, file);

	//Open record file in READONLY mode
	fits_open_file(&(file->fptr), filename, READONLY, status);
	CHECK_STATUS_RET(*status, file);

	//Move to the first hdu
	int hdu_type;
	fits_movabs_hdu(file->fptr,2, &hdu_type, status);
	CHECK_STATUS_RET(*status, file);

	//Get number of rows
	char comment[MAXMSG];
	fits_read_key(file->fptr, TINT, "NAXIS2", &(file->nrows), comment, status);
	CHECK_STATUS_RET(*status, file);

	//Get trigger_size
	fits_read_key(file->fptr, TINT, "TRIGGSZ", &(file->trigger_size), comment, status);
	CHECK_STATUS_RET(*status, file);

	//Get delta_t
	fits_read_key(file->fptr, TDOUBLE, "DELTAT", &(file->delta_t), comment, status);
	CHECK_STATUS_RET(*status, file);


	//Associate column numbers
	fits_get_colnum(file->fptr, CASEINSEN,"TIME", &(file->timeCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"ADC", &(file->trigCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"PIXID", &(file->pixIDCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"PH_ID", &(file->ph_idCol), status);
	CHECK_STATUS_RET(*status, file);

	return(file);

}

void triggerWithImpact(TESDataStream* const stream,TESGeneralParameters * par,
		TESInitStruct* init,float monoen,const char write_file,const char reconstruct,
		ReconstructInit* reconstruct_init,char* const tes_event_file,int event_list_size,
		int* const status){

	//Get parameters from structures
	char* const xmlfile = par->XMLFile;
	char* const impactlist = par->PixImpList;
	double tstart = par->tstart;
	double tstop = par->tstop;
	const int triggerSize = par->triggerSize;
	const int preBufferSize = par->preBufferSize;
	const double sampleFreq = init->det->SampleFreq;
	char clobber = par->clobber;
	const int pixlow = par->nlo;
	const int Npix = par->nhi-par->nlo+1;
	char* const tesTriggerFilename = par->tesTriggerFile;


	////////////////////////////////
	//Open output files depending on what we wish to do
	////////////////////////////////
	TesTriggerFile* outputFile = NULL;
	TesEventList* event_list = NULL;
	TesEventFile* out_event_file = NULL;
	if (write_file){
		outputFile = opennewTesTriggerFile(tesTriggerFilename,
				init->telescop,
				init->instrume,
				init->filter,
				init->ancrfile,
				init->respfile,
				xmlfile,
				impactlist,
				init->mjdref,
				init->timezero,
				tstart,
				tstop,
				triggerSize,
				preBufferSize,
				sampleFreq,
				clobber,
				status);
		CHECK_STATUS_VOID(*status);
	}
	if (reconstruct){
		//Build up TesEventList to recover the results of the reconstruction
		event_list = newTesEventList(status);
		allocateTesEventListTrigger(event_list,event_list_size,status);
		CHECK_STATUS_VOID(*status);

		//Open TesEventFile
		out_event_file = opennewTesEventFile(tes_event_file,
				init->telescop,
				init->instrume,
				init->filter,
				init->ancrfile,
				init->respfile,
				init->mjdref,
				init->timezero,
				tstart,
				tstop,
				clobber,
				status);
		CHECK_STATUS_VOID(*status);
	}
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



	//Allocation of array containing the TesRecords being constructed
	TesRecord ** records = malloc(Npix*sizeof(*records));
	if(records==NULL){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for records failed.");
		CHECK_STATUS_VOID(*status);
	}
	int ii;
	//Current time
	double t = tstart;
	long tlong = 0;
	//Current impact
	PixImpact impact;
	// Status of getNextImpactFromPixImpFile
	int piximpstatus=0;
	// Number of records found
	int nRecords = 0;
	// Number of triggers per pixel
	int* numberTrigger = malloc(Npix*sizeof(int));
	// Number of simulated pulses per pixel
	int* numberSimulated = malloc(Npix*sizeof(int));
	// Pixel iterator
	int pixNumber=0;
	//Array to save at which point we are while constructing a trigger
	int* positionInTrigger = malloc(Npix*sizeof(int));
	//Array to save the signal to record right after the end of a record
	//unsigned char* forceRecord = malloc(Npix*sizeof(*forceRecord));

	//Allocate phid lists
	PhIDList ** preBuffPhIDLists = malloc(Npix*sizeof(*preBuffPhIDLists)); //List containing the ph_ids of the current record of each pixel
	if(preBuffPhIDLists==NULL || numberTrigger == NULL || numberSimulated==NULL || positionInTrigger == NULL){// || forceRecord==NULL){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation during triggering initialization failed failed.");
		CHECK_STATUS_VOID(*status);
	}

	//Initializations
	for(ii=0; ii<Npix; ii++){
		records[ii] = newTesRecord(status);
		allocateTesRecord(records[ii],triggerSize,status);
		CHECK_STATUS_VOID(*status);
		positionInTrigger[ii]=-1;
		numberTrigger[ii]=0;
		numberSimulated[ii]=0;
		//forceRecord[ii]=0;
		preBuffPhIDLists[pixNumber] = newAllocatedPhIDList((int)((double)preBufferSize/triggerSize*MAXIMPACTNUMBER),1,status);
		CHECK_STATUS_VOID(*status);
	}

	//Loop over times
	int Nt = (int)((tstop-tstart)*sampleFreq);
	double old_impact_time=0;
	long old_ph_id=0;
	for (int tstep=0;tstep < Nt;tstep++){
		/* Get first pulse in correct time frame from the impact file */
		if (tstep==0) {
			do {
				piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
				CHECK_STATUS_VOID(*status);
			} while((impact.time<tstart) && piximpstatus );
		}

		/* If the pulse occurs in a time bin away from preBufferSize samples of the current time, is in an active pixel, */
		/* and we are not recording, put flag to record*/
		while ((piximpstatus>0) && (impact.time<t+(((double)preBufferSize+1.0)/sampleFreq)) && (impact.time < tstop)) {
			if((impact.pixID>=pixlow) && (impact.pixID<(pixlow+Npix))){
				if (positionInTrigger[impact.pixID-pixlow]==-1) {
					positionInTrigger[impact.pixID-pixlow]=0;

					//Check the wait list to see if some previous impacts may in fact be in this new record
					while(popPhID(preBuffPhIDLists[impact.pixID-pixlow],&old_ph_id,&old_impact_time,status)){
						if(t<old_impact_time){
							numberTrigger[impact.pixID-pixlow]++;
							appendPhID(records[impact.pixID-pixlow]->phid_list,old_ph_id,0.,status);
						}
					}
				}
				//If the pulse is gonna be in the currently recorded trigger, increase numberTrigger and save ph_id (no need to save time)
				if (impact.time< (t+(double)(triggerSize-positionInTrigger[impact.pixID-pixlow])/sampleFreq)){
					numberTrigger[impact.pixID-pixlow]++;
					appendPhID(records[impact.pixID-pixlow]->phid_list,impact.ph_id,0.,status);
				// If the pulse has a chance to fall in the missed part of the data, save its time and ph_id in case another "saves" it
				} else if (impact.time<(t+(double)(triggerSize-positionInTrigger[impact.pixID-pixlow]+preBufferSize)/sampleFreq)){
					appendPhID(preBuffPhIDLists[impact.pixID-pixlow],impact.ph_id,impact.time,status);
					//forceRecord[impact.pixID-pixlow]=1;
				}
				CHECK_STATUS_VOID(*status);
				numberSimulated[impact.pixID-pixlow]++;
				numberTrigger[impact.pixID-pixlow]++;
			}
			CHECK_STATUS_VOID(*status);
			piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
			CHECK_STATUS_VOID(*status);
		}

		// Record the ADC values in the triggers if it is necessary
		for (pixNumber=0;pixNumber<Npix;pixNumber++) {
			if ((positionInTrigger[pixNumber]>=0) && (positionInTrigger[pixNumber]<triggerSize)) {
				records[pixNumber]->adc_array[positionInTrigger[pixNumber]] = stream->adc_value[(int)round((t-tstartTES)*sampleFreq)][pixNumber];
				records[pixNumber]->adc_double[positionInTrigger[pixNumber]] = (double) records[pixNumber]->adc_array[positionInTrigger[pixNumber]];
				positionInTrigger[pixNumber]++;
			}
			//If Trigger is full, record in output FITS file
			if (positionInTrigger[pixNumber]==triggerSize){
				records[pixNumber]->time = t-(triggerSize-1)/sampleFreq;
				records[pixNumber]->pixid = pixNumber+pixlow+1;
				if(write_file){
					writeRecord(outputFile,records[pixNumber],status);
					CHECK_STATUS_VOID(*status);
					nRecords++;//count records
				}
				if(reconstruct){
					reconstructRecord(records[pixNumber],event_list,reconstruct_init,status);
					saveEventListToFile(out_event_file,event_list,records[pixNumber]->time,1./sampleFreq,records[pixNumber]->pixid,status);
					CHECK_STATUS_VOID(*status);

					//Reinitialize event list
					event_list->index=0;
					nRecords++;//count records
				}

				//Reinitialize for next record
				positionInTrigger[pixNumber]=-1;
				records[pixNumber]->phid_list->index=0;
				/*if(forceRecord[pixNumber]){
					positionInTrigger[pixNumber]=0;
					forceRecord[pixNumber] = 0;
				} else {
					positionInTrigger[pixNumber]=-1;
				}*/
			}
		}
		fflush(stdout);
		tlong++;
		t=tstart+tlong/sampleFreq;
	}

	//Record incomplete record if there is one
	for (pixNumber=0;pixNumber<Npix;pixNumber++) {
		if (positionInTrigger[pixNumber]>0){
			for (int j=positionInTrigger[pixNumber];j<triggerSize;j++) {
				records[pixNumber]->adc_array[j] = 0;
				records[pixNumber]->adc_double[j] = 0;
			}
			records[pixNumber]->time = t-(triggerSize-1)/sampleFreq;
			records[pixNumber]->pixid = pixNumber+pixlow+1;
			if(write_file){
				writeRecord(outputFile,records[pixNumber],status);
				CHECK_STATUS_VOID(*status);
				nRecords++;//count records
			}
			if(reconstruct){
				reconstructRecord(records[pixNumber],event_list,reconstruct_init,status);
				saveEventListToFile(out_event_file,event_list,records[pixNumber]->time,1./sampleFreq,records[pixNumber]->pixid,status);
				CHECK_STATUS_VOID(*status);
				nRecords++;//count records
			}
		}
	}

	//If there is no trigger, print WARNING. Still compute numberSimulated.
	if (nRecords==0) {
		puts("WARNING: No trigger found. Check in impact file that there does exist an event inside the simulation time in the given pixels");
		//Reinitialize impfile
		impfile->row=0;
		//Get first impact after tstart
		do {
			piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
			CHECK_STATUS_VOID(*status);
		} while(impact.time<tstart);

		//Iterate over the impacts
		while ((piximpstatus>0) && (impact.time<tstop)){
			if ((impact.pixID>=pixlow) && (impact.pixID<(pixlow+Npix))){
				numberSimulated[impact.pixID-pixlow]++;
			}
			piximpstatus=getNextImpactFromPixImpFile(impfile,&impact,status);
			CHECK_STATUS_VOID(*status);
		}
	}

	//Save keywords
	int firstpix = pixlow+1;
	int lastpix = pixlow+Npix;
	int numberpix = Npix;
	if(write_file){
		saveTriggerKeywords(outputFile->fptr,firstpix,lastpix,numberpix,monoen,
				numberSimulated,numberTrigger,status);
	}
	if(reconstruct){
		saveTriggerKeywords(out_event_file->fptr,firstpix,lastpix,numberpix,monoen,
				numberSimulated,numberTrigger,status);
	}

	//Free memory
	freeTesTriggerFile(&(outputFile), status);
	freePixImpFile(&impfile, status);
	free(numberSimulated);
	free(numberTrigger);
	free(positionInTrigger);
	freeTesEventFile(out_event_file,status);
	freeTesEventList(event_list);
	//free(forceRecord);
	for (int ii = 0 ; ii < Npix;ii++){
		freeTesRecord(&(records[ii]));
	}
	if(NULL!=preBuffPhIDLists){
		for(ii=0;ii<Npix;ii++){
			freePhIDList(preBuffPhIDLists[ii]);
		}
		free(preBuffPhIDLists);
	}

}

/** Populates a TesRecord structure with the next record */
int getNextRecord(TesTriggerFile* const file,TesRecord* record,int* const status){
	int anynul=0;
	if (NULL==file || NULL==file->fptr) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("No opened trigger file to read from");
		CHECK_STATUS_RET(*status,0);
	}


	if (file->row<=file->nrows) {
		fits_read_col(file->fptr, TUSHORT, file->trigCol,
					  file->row,1,file->trigger_size,0,record->adc_array, &anynul,status);
		CHECK_STATUS_RET(*status,0);

//		fits_read_col(file->fptr, TLONG, file->ph_idCol,
//					  file->row,1,MAXIMPACTNUMBER,0,record->phid_array, &anynul,status);
//		CHECK_STATUS_RET(*status,0);

		fits_read_col(file->fptr, TLONG, file->pixIDCol,
					  file->row,1,1,0,&(record->pixid), &anynul,status);
		CHECK_STATUS_RET(*status,0);

		fits_read_col(file->fptr, TDOUBLE, file->timeCol,
					  file->row,1,1,0,&(record->time), &anynul,status);
		CHECK_STATUS_RET(*status,0);

		for (int i=0 ; i < file->trigger_size ; i++) {
			record->adc_double[i]= (double) (record->adc_array[i]);
		}

		file->row++;
		return(1);
	} else {
		return(0);
	}


}

/** Writes a record to a file */
void writeRecord(TesTriggerFile* outputFile,TesRecord* record,int* const status){
	fits_write_col(outputFile->fptr, TDOUBLE, outputFile->timeCol,
			outputFile->row, 1, 1, &(record->time), status);
	fits_write_col(outputFile->fptr, TUSHORT, outputFile->trigCol,
			outputFile->row, 1, record->trigger_size,record->adc_array, status);
	fits_write_col(outputFile->fptr, TLONG, outputFile->pixIDCol,
			outputFile->row, 1, 1, &(record->pixid), status);
	fits_write_col(outputFile->fptr, TLONG, outputFile->ph_idCol,
			outputFile->row, 1,record->phid_list->index,record->phid_list->phid_array, status);
	outputFile->nrows++;
	outputFile->row++;
}

