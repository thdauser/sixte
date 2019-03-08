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


   Copyright 2015 Philippe Peille, IRAP
   Copyright 2016-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "testrigger.h"

/** Save pixels, NES/NET and monoen keywords to the given FITS file */
void saveTriggerKeywords(fitsfile* fptr,int firstpix,int lastpix,int numberpix,float monoen,
		int* const numberSimulated,int* const numberTrigger,int* const status){
	int old_firstpix,old_lastpix,old_numberpix;
	long net_tot,nes_tot;
	float old_monoen;
	char comment[MAXMSG];

	//Read old keywords
	fits_read_key(fptr, TINT, "FIRSTPIX", &old_firstpix, comment, status);
	fits_read_key(fptr, TINT, "LASTPIX", &old_lastpix, comment, status);
	fits_read_key(fptr, TINT, "NPIX", &old_numberpix, comment, status);
	fits_read_key(fptr, TFLOAT, "MONOEN", &old_monoen, comment, status);
	fits_read_key(fptr, TLONG, "NESTOT", &nes_tot, comment, status);
	fits_read_key(fptr, TLONG, "NETTOT", &net_tot, comment, status);


	//If this is actually the first time the header is set, simply save the given keywords
	if(old_monoen==-1 || old_firstpix==0 || old_lastpix==0 || old_numberpix==0 || nes_tot==0 || net_tot==0){
		if(!(old_monoen==-1 && old_firstpix==0 && old_lastpix==0 && old_numberpix==0 && nes_tot==0 && net_tot==0)){
			*status=EXIT_FAILURE;
			SIXT_ERROR("One of the trigger keywords of the given file has been changed without the others. File probably corrupted, abort");
			CHECK_STATUS_VOID(*status);
		}
		fits_update_key(fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in record file", status);
		fits_update_key(fptr, TINT, "LASTPIX", &lastpix, "Last pixel in record file", status);
		fits_update_key(fptr, TINT, "NPIX", &numberpix, "Number of pixels in record file", status);
		fits_update_key(fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);
	} else{
		if(firstpix<old_firstpix) fits_update_key(fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in record file", status);
		if(lastpix>old_lastpix) fits_update_key(fptr, TINT, "LASTPIX", &lastpix, "Last pixel in record file", status);
		old_numberpix = old_numberpix+numberpix;
		fits_update_key(fptr, TINT, "NPIX", &old_numberpix, "Number of pixels in record file", status);
		if(!(old_monoen && monoen && old_monoen==monoen)) {
			old_monoen=0;
			fits_update_key(fptr, TFLOAT, "MONOEN", &old_monoen, "Monochromatic energy of photons [keV]", status);
		}
	}
	CHECK_STATUS_VOID(*status);

	//Save keywords with number of counts
	char keyword[9];
	for (int pixNumber=0;pixNumber<numberpix;pixNumber++) {
		sprintf(keyword,"NES%05d",pixNumber+firstpix);
		fits_update_key(fptr, TINT, keyword, &(numberSimulated[pixNumber]), "Number of simulated pulses", status);
		nes_tot+=numberSimulated[pixNumber];
		sprintf(keyword,"NET%05d",pixNumber+firstpix);
		fits_update_key(fptr, TINT, keyword, &(numberTrigger[pixNumber]), "Number of triggered pulses", status);
		net_tot+=numberTrigger[pixNumber];
	}
	CHECK_STATUS_VOID(*status);
	fits_update_key(fptr, TLONG, "NESTOT", &nes_tot, "Total number of events simulated", status);
	fits_update_key(fptr, TLONG, "NETTOT", &net_tot, "Total number of events actually triggered", status);
	CHECK_STATUS_VOID(*status);
}

void copyTriggerKeywords(fitsfile* fptr,fitsfile* fptr2,int* const status){
	int firstpix,lastpix,numberpix;
	long net_tot,nes_tot,nes,net;
	float monoen;
	char comment[MAXMSG];
	char keyword[9];

	//Read keywords in first file and dump it in new one
	fits_read_key(fptr, TINT, "FIRSTPIX", &firstpix, comment, status);
	fits_update_key(fptr2, TINT, "FIRSTPIX", &firstpix, comment, status);
	CHECK_STATUS_VOID(*status);

	fits_read_key(fptr, TINT, "LASTPIX", &lastpix, comment, status);
	fits_update_key(fptr2, TINT, "LASTPIX", &lastpix, comment, status);
	CHECK_STATUS_VOID(*status);

	fits_read_key(fptr, TINT, "NPIX", &numberpix, comment, status);
	fits_update_key(fptr2, TINT, "NPIX", &numberpix, comment, status);
	CHECK_STATUS_VOID(*status);

	fits_read_key(fptr, TFLOAT, "MONOEN", &monoen, comment, status);
	fits_update_key(fptr2, TFLOAT, "MONOEN", &monoen, comment, status);
	CHECK_STATUS_VOID(*status);

	fits_read_key(fptr, TLONG, "NESTOT", &nes_tot, comment, status);
	if (*status==KEY_NO_EXIST){ // to allow the reading of old record files that did not have this keyword
		*status=EXIT_SUCCESS;
	} else{
		fits_update_key(fptr2, TLONG, "NESTOT", &nes_tot, comment, status);
	}
	CHECK_STATUS_VOID(*status);

	fits_read_key(fptr, TLONG, "NETTOT", &net_tot, comment, status);
	if (*status==KEY_NO_EXIST){ // to allow the reading of old record files that did not have this keyword
		*status=EXIT_SUCCESS;
	} else{
		fits_update_key(fptr2, TLONG, "NETTOT", &net_tot, comment, status);
	}
	CHECK_STATUS_VOID(*status);

	int status_to_check_key=EXIT_SUCCESS;
	for (int i=firstpix; i <=lastpix;i++){
		sprintf(keyword,"NES%05d",i);
		fits_read_key(fptr, TLONG, keyword, &nes, comment, &status_to_check_key);
		if(status_to_check_key==EXIT_SUCCESS){
			fits_update_key(fptr2, TLONG, keyword, &nes, comment, status);
		} else if (status_to_check_key!=KEY_NO_EXIST){
			*status = status_to_check_key;
		} else {
			status_to_check_key=EXIT_SUCCESS;
		}
		CHECK_STATUS_VOID(*status);

		sprintf(keyword,"NET%05d",i);
		fits_read_key(fptr, TLONG, keyword, &net, comment, &status_to_check_key);
		if(status_to_check_key==EXIT_SUCCESS){
			fits_update_key(fptr2, TLONG, keyword, &net, comment, status);
		} else if (status_to_check_key!=KEY_NO_EXIST){
			*status = status_to_check_key;
		} else {
			status_to_check_key=EXIT_SUCCESS;
		}
		CHECK_STATUS_VOID(*status);
	}

}

void triggerWithImpact(TESDataStream* const stream,TESGeneralParameters * par,
		TESInitStruct* init,float monoen,ReconstructInit* reconstruct_init,int event_list_size,
		const char identify,int* const status){

	//Get parameters from structures
	char* const impactlist = par->PixImpList;
	double tstart = init->tstart;
	double tstop = init->tstop;
	const int triggerSize = par->triggerSize;
	const int preBufferSize = par->preBufferSize;
	const double sampleFreq = init->det->SampleFreq;
	const int pixlow = par->nlo;
	const int Npix = par->nhi-par->nlo+1;


	////////////////////////////////
	//Open output files depending on what we wish to do
	////////////////////////////////
	TesEventList* event_list = NULL;
	SixtStdKeywords* keywords = buildSixtStdKeywords(init->telescop,init->instrume,init->filter,
			init->ancrfile,init->respfile,"NONE",init->mjdref,init->timezero,tstart,tstop,status);

	if (par->Reconstruct){
		//Build up TesEventList to recover the results of the reconstruction
		event_list = newTesEventList(status);
		allocateTesEventListTrigger(event_list,event_list_size,status);
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
		allocateTesRecord(records[ii],triggerSize,1./sampleFreq,identify,status);
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
							appendPhID(records[impact.pixID-pixlow]->phid_list,old_ph_id,old_impact_time,status);
						}
					}
				}
				//If the pulse is gonna be in the currently recorded trigger, increase numberTrigger and save ph_id (no need to save time)
				if (impact.time< (t+(double)(triggerSize-positionInTrigger[impact.pixID-pixlow])/sampleFreq)){
					numberTrigger[impact.pixID-pixlow]++;
					appendPhID(records[impact.pixID-pixlow]->phid_list,impact.ph_id,impact.time,status);
				// If the pulse has a chance to fall in the missed part of the data, save its time and ph_id in case another "saves" it
				} else if (impact.time<(t+(double)(triggerSize-positionInTrigger[impact.pixID-pixlow]+preBufferSize)/sampleFreq)){
					appendPhID(preBuffPhIDLists[impact.pixID-pixlow],impact.ph_id,impact.time,status);
					//forceRecord[impact.pixID-pixlow]=1;
				}
				CHECK_STATUS_VOID(*status);
				numberSimulated[impact.pixID-pixlow]++;
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
				nRecords++;//count records
				if(par->WriteRecordFile){
					writeRecord(init->record_file,records[pixNumber],status);
					CHECK_STATUS_VOID(*status);
				}
				if(par->Reconstruct){
					reconstructRecord(records[pixNumber],event_list,reconstruct_init,identify,status);
					saveEventListToFile(init->event_file,event_list,records[pixNumber]->time,1./sampleFreq,records[pixNumber]->pixid,status);
					CHECK_STATUS_VOID(*status);

					//Reinitialize event list
					event_list->index=0;
				}

				//Reinitialize for next record
				positionInTrigger[pixNumber]=-1;
				records[pixNumber]->phid_list->index=0;
				records[pixNumber]->phid_list->n_elements=0;
				/*if(forceRecord[pixNumber]){
					positionInTrigger[pixNumber]=0;
					forceRecord[pixNumber] = 0;
				} else {
					positionInTrigger[pixNumber]=-1;
				}*/
			}
		}
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
			records[pixNumber]->time = t-(positionInTrigger[pixNumber])/sampleFreq;
			records[pixNumber]->pixid = pixNumber+pixlow+1;
			nRecords++;//count records
			if(par->WriteRecordFile){
				writeRecord(init->record_file,records[pixNumber],status);
				CHECK_STATUS_VOID(*status);
			}
			if(par->Reconstruct){
				reconstructRecord(records[pixNumber],event_list,reconstruct_init,identify,status);
				saveEventListToFile(init->event_file,event_list,records[pixNumber]->time,1./sampleFreq,records[pixNumber]->pixid,status);
				CHECK_STATUS_VOID(*status);
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
	if(par->WriteRecordFile){
		saveTriggerKeywords(init->record_file->fptr,firstpix,lastpix,numberpix,monoen,
				numberSimulated,numberTrigger,status);
	}
	if(par->Reconstruct){
		saveTriggerKeywords(init->event_file->fptr,firstpix,lastpix,numberpix,monoen,
				numberSimulated,numberTrigger,status);
	}

	//Free memory
	freePixImpFile(&impfile, status);
	free(numberSimulated);
	free(numberTrigger);
	free(positionInTrigger);
	freeTesEventList(event_list);
	//free(forceRecord);
	for (int ii = 0 ; ii < Npix;ii++){
		freeTesRecord(&(records[ii]));
	}
	free(records);
	if(NULL!=preBuffPhIDLists){
		for(ii=0;ii<Npix;ii++){
			freePhIDList(preBuffPhIDLists[ii]);
		}
		free(preBuffPhIDLists);
	}
	freeSixtStdKeywords(keywords);

}
