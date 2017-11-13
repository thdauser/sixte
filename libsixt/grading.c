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


   Copyright 2016 Philippe Peille, IRAP; Thomas Dauser, ECAP; Edoardo Cucchetti, IRAP;
*/

#include "grading.h"

/** given grade1 and grade 2, make a decision about the high/mid/los res events **/
int makeGrading(long grade1,long grade2,AdvPix* pixel){

	for (int i=0;i<pixel->ngrades;i++){
		if((grade1 >= pixel->grades[i].gradelim_post) && (grade2>=pixel->grades[i].gradelim_pre)){
			return i;
		}
	}
	return -1;

}

/** calculate the grading in samples from the a given impact, and its previous and next impact **/
void calcGradingTimes(double sample_length, gradingTimeStruct pnt,long *grade1, long *grade2, int* status){

	if ((pnt.next<0) || (pnt.next - pnt.current > sample_length*DEFAULTGOODSAMPLE)){
		*grade1 = DEFAULTGOODSAMPLE;
	} else {
		*grade1 = floor ((pnt.next - pnt.current)/sample_length);
	}
	if ((pnt.previous<0) || (pnt.current - pnt.previous > sample_length*DEFAULTGOODSAMPLE)){
		*grade2 = DEFAULTGOODSAMPLE;
	} else {
		*grade2 = floor ((pnt.current - pnt.previous)/sample_length);
	}

	if ( (*grade1<0) || (*grade2<0) ){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Grading of impacts failed in calcGradingTimes");
		printf("\n *** times are previous=%e current=%e %e=next\n",pnt.previous,pnt.current,pnt.next);
		return;
	}
}

/** writes the grading to an existing piximpact file **/
void writeGrading2PixImpactFile(AdvDet *det,PixImpFile *piximpacfile,int *status){

	long grade1, grade2;
	const double sample_length = 1./(det->SampleFreq);

	PixImpact impact;

	pixGrade pnt[det->npix];
	for (int ii=0;ii<det->npix;ii++){
		pnt[ii].times = NULL;
	}

	int id = -1;
	while (getNextImpactFromPixImpFile(piximpacfile,&impact,status)){
		id = impact.pixID;
		if (pnt[id].times == NULL){
			pnt[id].times = (gradingTimeStruct*) malloc (sizeof(gradingTimeStruct));
			CHECK_NULL_VOID(pnt[id].times,*status,"malloc failed");

			pnt[id].times->previous = -1.0;
			pnt[id].times->current = impact.time;
			pnt[id].times->next = -1.0;
			pnt[id].totalenergy = impact.energy;

			pnt[id].row = piximpacfile->row;
		} else {
			// save the time
			pnt[id].times->next = impact.time;

			// calculate grading
			calcGradingTimes(sample_length,*(pnt[id].times),&grade1,&grade2,status);
			CHECK_STATUS_VOID(*status);

			// treat pileup case
			// TODO : add pileup limit as a pixel parameter (for the moment, equal to one sample)
			if (grade1==0){
				// add energy to total instead of replacing
				pnt[id].totalenergy+=impact.energy;
				grade1=-1;
				// add grading to file
				updateGradingPixImp(piximpacfile, pnt[id].row, grade1, grade2,0,status);

				// do not change grading times in this case as this is a pileup (this event will no be saved in the end)

			} else {
				// add grading to file
				updateGradingPixImp(piximpacfile, pnt[id].row, grade1, grade2,pnt[id].totalenergy,status);

				// move one ahead
				pnt[id].totalenergy=impact.energy;
				pnt[id].times->previous = pnt[id].times->current;
				pnt[id].times->current = pnt[id].times->next;
			}
			CHECK_STATUS_VOID(*status);

			// finishing moving ahead
			pnt[id].times->next = -1.0;
			pnt[id].row = piximpacfile->row;
		}
	}

	// now we need to clean the impact-array and write the remaining impact grade to the file
	// TODO: only for over the pixels that were hit
	for (int ii=0;ii<det->npix;ii++){
		if (pnt[ii].times != NULL){
			calcGradingTimes(sample_length,*(pnt[ii].times),&grade1,&grade2,status);
			updateGradingPixImp(piximpacfile, pnt[ii].row, grade1, grade2,pnt[ii].totalenergy,status);
			free(pnt[ii].times);
		}
	}

}


/** Process the impacts contained in the piximpacts file with the RMF method */
void processImpactsWithRMF(AdvDet* det,PixImpFile* piximpacfile,TesEventFile* event_file,int* const status){
	long channel;
	struct RMF* rmf;
	int grading;
	int grading_index;

	// initialize
	PixImpact impact;

	while(getNextImpactFromPixImpFile(piximpacfile,&impact,status)){

		//////////////////////////////////
		//Copied and adapted from addGenDetPhotonImpact
		//////////////////////////////////


		if (impact.grade1==-1){ // this is a pile-up event that was added to the following event, i.e. not detected
			grading = PILEUP;
			impact.energy=0.;
		} else { // this is an event that was actually detected
			grading_index = makeGrading(impact.grade1,impact.grade2,&(det->pix[impact.pixID]));

			if (grading_index>=0){

				rmf = det->pix[impact.pixID].grades[grading_index].rmf;
				grading = det->pix[impact.pixID].grades[grading_index].value;

				// Determine the measured detector channel (PI channel) according
				// to the RMF.
				// The channel is obtained from the RMF using the corresponding
				// HEAdas routine which is based on drawing a random number.
				returnRMFChannel(rmf, impact.totalenergy, &channel); //use total energy here to take pileup into account

				// Check if the photon is really measured. If the PI channel
				// returned by the HEAdas RMF function is '-1', the photon is not
				// detected. This should not happen as the rmf is supposedly
				// normalized
				if (channel<rmf->FirstChannel) {
					// flag as invalid (seemed better than discarding)
					char msg[MAXMSG];
					sprintf(msg,"Impact found as not detected due to failure during RMF based energy allocation for energy %g keV and phi_id %ld",impact.totalenergy,impact.ph_id);
					SIXT_WARNING(msg);
					impact.energy=0.;
					grading=INVGRADE;
				} else {

					// Determine the corresponding detected energy.
					// NOTE: In this simulation the collected charge is represented
					// by the nominal photon energy [keV], which corresponds to the
					// PI channel according to the EBOUNDS table.
					impact.energy=getEBOUNDSEnergy(channel, rmf, status);
					//printf("%f ",impact.energy);
					CHECK_STATUS_VOID(*status);
				}
			} else {
				impact.energy=0.;
				grading=INVGRADE;
			}
		}
		assert(impact.energy>=0.);
		// Add final impact to event file
		addRMFImpact(event_file,&impact,impact.grade1,impact.grade2,grading,0,0.,status);
		CHECK_STATUS_VOID(*status);

	}

}

/** Processes the impacts, including crosstalk and RMF energy randomization **/
void impactsToEvents(AdvDet *det,PixImpFile *piximpactfile,TesEventFile* event_file,int save_crosstalk, FILE* progressfile, int* const status){

	const double sample_length = 1./(det->SampleFreq);

	PixImpact impact;
	GradeProxy grade_proxys[det->npix];
	for (int ii=0;ii<det->npix;ii++){
		grade_proxys[ii].impact = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_VOID_STATUS(grade_proxys[ii].impact,*status);
		grade_proxys[ii].impact->energy=0.;
		grade_proxys[ii].xtalk_proxy = NULL;
		grade_proxys[ii].times = NULL;
		grade_proxys[ii].row = event_file->row;
		grade_proxys[ii].is_first=0;
		grade_proxys[ii].crosstalk_energy=0.;
		grade_proxys[ii].nb_crosstalk_influence=0;
	}
	int id = -1;

	//Get the number of impacts
	unsigned int total_length=piximpactfile->nrows;
	unsigned int ndone=0;
	unsigned int progress=0;
	if (NULL==progressfile) {
		headas_chat(2, "\r%.0lf %%", 0.);
		fflush(NULL);
	} else {
		rewind(progressfile);
		fprintf(progressfile, "%.2lf", 0.);
		fflush(progressfile);
	}

	// Iterate over impacts
	while (getNextImpactFromPixImpFile(piximpactfile,&impact,status)){
		id = impact.pixID;
		ndone+=1;

		// First thing: we update the crosstalk proxy of the impacts only with the type. Basically, we say if a given photon
		// hits pixel A, we know that pixels B,C,D... will have crosstalk but we do not know the energy (depends on their grading)
		// intermod crosstalk
		if (det->crosstalk_imod_table!=NULL){
			applyIntermodCrossTalk(grade_proxys,&impact,det,status);
			CHECK_STATUS_VOID(*status);
		}

		// thermal crosstalk
		if (det->pix[id].thermal_cross_talk !=NULL){
			//1st matrix given as only the indices of the pixels are needed (same for all grades obviously!)
			applyMatrixCrossTalk(det->pix[id].thermal_cross_talk,grade_proxys,&impact,det,status);
			CHECK_STATUS_VOID(*status);
		}
		// electrical crosstalk
		if (det->pix[id].electrical_cross_talk !=NULL){
			//1st matrix given as only the indices of the pixels are needed (same for all grades obviously!)
			applyMatrixEnerdepCrossTalk(&(det->pix[id].electrical_cross_talk[0]),grade_proxys,&impact,det,status);
			CHECK_STATUS_VOID(*status);
		}

		// Once all the crosstalk potential impacts were stored, go on and process the impact itself
		// Process impact and get its grade
		processGradedEvent(&(grade_proxys[id]),sample_length,&impact,det,event_file,0,save_crosstalk,-99,status); //Grade here given only for comparison
		CHECK_STATUS_VOID(*status);

		// Program progress output.
		while((unsigned int)((ndone*100./total_length)>progress)) {
			progress++;
			if (NULL==progressfile) {
				headas_chat(2, "\r%.0lf %%", progress*1.);
				fflush(NULL);
			} else {
				rewind(progressfile);
				fprintf(progressfile, "%.2lf", progress*1./100.);
				fflush(progressfile);
			}
		}
	}
	printf("\n");
	// now we need to clean the remaining events (maximum grade as no more 'real' impact on pixel)
	for (int ii=0;ii<det->npix;ii++){
		if (grade_proxys[ii].times != NULL){
			int is_crosstalk=0;
			if (grade_proxys[ii].impact->ph_id<0){
				is_crosstalk=1;
			}
			processGradedEvent(&(grade_proxys[ii]),sample_length,NULL,det,event_file,is_crosstalk,save_crosstalk,-99,status);
			free(grade_proxys[ii].times);
		}
		free(grade_proxys[ii].impact);
		freeCrosstalkProxy(&(grade_proxys[ii].xtalk_proxy));
	}
}

// 4d & 3d interpolation routine: expecting array with dimensions arr[][][][]:
//  - fac and iarr has the length of the dimension (ndim=4 in this case)
//  - always giving val[ii] = arr[iarr[0]]... and val[ii+1] = arr[iarr[0]+1]...
//  TODO: Make it a general routine
//double interp_lin_ndim(double**** arr, int* iarr, double* fac){
double interp_lin_ndim(void* inp_arr, int* iarr, double* fac, int ndim){

	double tmp_3d[8];

	if (ndim==4){ // ***** 4D ****

		double**** arr = (double****) inp_arr;
		for (int ii=0;ii<2;ii++){
			for (int jj=0;jj<2;jj++){
				for (int kk=0;kk<2;kk++){
					tmp_3d[ (ii*2+jj)*2 + kk ] =
							arr[iarr[0]+ii][iarr[1]+jj][iarr[2]+kk][iarr[3]]  *(1-fac[3]) +
							arr[iarr[0]+ii][iarr[1]+jj][iarr[2]+kk][iarr[3]+1]*(  fac[3]) ;

				}
			}
		}
	} else if (ndim==3){ // ***** 3D ****

		double*** arr = (double***) inp_arr;
		for (int ii=0;ii<2;ii++){
			for (int jj=0;jj<2;jj++){
				for (int kk=0;kk<2;kk++){
					tmp_3d[ (ii*2+jj)*2 + kk ] =
							arr[iarr[0]+ii][iarr[1]+jj][iarr[2]+kk];
				}
			}
		}
	} else{
		printf(" *** error: interpolation not implmentend for %i dimensions",ndim);
		return 0.0;
	}

	double tmp_2d[4];
	for (int ii=0;ii<2;ii++){
		for (int jj=0;jj<2;jj++){
					tmp_2d[ ii*2+jj  ] =
							tmp_3d[ (ii*2+jj)*2     ]*(1-fac[2]) +
							tmp_3d[ (ii*2+jj)*2 + 1 ]*(  fac[2]) ;
		}
	}

	double tmp_1d[2];
	for (int ii=0;ii<2;ii++){
					tmp_1d[ ii  ] =
							tmp_2d[ ii*2     ]*(1-fac[1]) +
							tmp_2d[ ii*2 + 1 ]*(  fac[1]) ;
	}

	double tmp_val;
	tmp_val =	tmp_1d[ 0 ]*(1-fac[0]) +
			    tmp_1d[ 1 ]*(  fac[0]) ;

	return tmp_val;

}


/** Apply matrix cross talk: create new events on concerned pixels if corresponding
    energy is above the detection threshold, affect previous event otherwise */
void applyIntermodCrossTalk(GradeProxy* grade_proxys,PixImpact* impact, AdvDet* det, int* const status){

    // This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
    crosstalk_impact.detposition.x = 0.;
    crosstalk_impact.detposition.y = 0.;
    crosstalk_impact.pixposition.x = 0.;
    crosstalk_impact.pixposition.y = 0.;

	// find out the channel of pixel the photon hit
	Channel* active_chan = det->pix[impact->pixID].channel;

	// loop over all pixels in the channel to see if there is a hit
	for (int ii=0; ii < active_chan->num_pixels; ii++){
		int active_ind = active_chan->pixels[ii]->pindex;

		// do not count crosstalk with ourself!
		if (active_ind==impact->pixID) continue;

		// now we just store the event in the proxy and go (analyzed later on)
		double df=get_imod_df(det->pix[active_ind].freq,det->pix[impact->pixID].freq,status);
		crosstalk_impact.energy = impact->energy; 			//We store the energy pixel
		crosstalk_impact.pixID = impact->pixID; 			//We store the impacts pixel
		crosstalk_impact.time = impact->time;            //  the event has the time of the perturber
		crosstalk_impact.ph_id = -impact->ph_id;         //  negative ID to show it's a crosstalk event
		crosstalk_impact.src_id = impact->src_id;
		crosstalk_impact.weight_index = 0;				//Useless info for this type of crosstalk
		//Now we add the event to the proxy of the grade proxy
		addCrosstalkEvent(&(grade_proxys[active_ind]),&crosstalk_impact,det,IMODCTK,df,status);
		CHECK_STATUS_VOID(*status);
	}
	//}
}


/** Apply matrix cross talk: create new events on concerned pixels if corresponding energy is above the detection threshold, affect previous event otherwise */
void applyMatrixCrossTalk(MatrixCrossTalk* cross_talk,GradeProxy* grade_proxys,PixImpact* impact,AdvDet* det, int* const status){

	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;

	// Iterate over affected pixels we just store the cross-talk impact and go
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){
		crosstalk_impact.energy = impact->energy; // We store the energy of the ctk event as it does not depend on grading
		crosstalk_impact.pixID = impact->pixID; //We store the perturber pixel
		crosstalk_impact.time = impact->time;
		crosstalk_impact.ph_id = -impact->ph_id;
		crosstalk_impact.src_id = impact->src_id;
		crosstalk_impact.weight_index = ii; //Store which of the neighbours to look into
		double df = get_imod_df(det->pix[cross_talk->cross_talk_pixels[ii]->pindex].freq,det->pix[impact->pixID].freq,status);
		addCrosstalkEvent(&(grade_proxys[cross_talk->cross_talk_pixels[ii]->pindex]),&crosstalk_impact,det,THERCTK,df,status);
		CHECK_STATUS_VOID(*status);
	}
}

/** Same as applyMatricCrosstalk , but now the weights are energy dependent */
void applyMatrixEnerdepCrossTalk(MatrixEnerdepCrossTalk* cross_talk,GradeProxy* grade_proxys, PixImpact* impact,AdvDet* det, int* const status){

	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;
	// Iterate over affected pixels, copy and go
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){
		crosstalk_impact.energy = impact->energy; //We store the energy pixel
		crosstalk_impact.pixID = impact->pixID; //We store the impact pixel "perturber"
		crosstalk_impact.time = impact->time;
		crosstalk_impact.ph_id = -impact->ph_id;
		crosstalk_impact.src_id = impact->src_id;
		crosstalk_impact.weight_index = ii; //Store the index of this given pixel to apply ctk later
		double df = get_imod_df(det->pix[cross_talk->cross_talk_pixels[ii]->pindex].freq,det->pix[impact->pixID].freq,status);
		addCrosstalkEvent(&(grade_proxys[cross_talk->cross_talk_pixels[ii]->pindex]),&crosstalk_impact,det,ELECCTK,df,status);
		CHECK_STATUS_VOID(*status);
	}
}


// We just store the event in the proxy and we shall treat it afterwards once we have the grading
void addCrosstalkEvent(GradeProxy* grade_proxy,PixImpact* impact, AdvDet* det, int type, double df, int* const status){
	// If needed we save the ctk. /!\IT WILL HAVE THE PERTURBER ENERGY AND INDEX!!!
	if(grade_proxy->xtalk_proxy==NULL){
		grade_proxy->xtalk_proxy = newCrosstalkProxy(status);
		CHECK_STATUS_VOID(*status);
	}
	if (grade_proxy->times==NULL){
		addCrosstalk2Proxy(grade_proxy->xtalk_proxy,-1.0,impact,det,type,df,status);
	}
	else{
		addCrosstalk2Proxy(grade_proxy->xtalk_proxy,grade_proxy->times->current,impact,det,type,df,status);
	}
}

/** Processes a graded event : update grading proxy and save previous event */
void processGradedEvent(GradeProxy* grade_proxy, const double sample_length,PixImpact* next_impact,
		AdvDet* det,TesEventFile* event_file, int is_crosstalk, int save_crosstalk, int grdcmp, int* const status){
	long grade1, grade2;
	int grading;
	int grading_index;
	struct RMF* rmf;
	long channel;
	int is_trigger=0; //Checks if a pileup occurs when computing cross-talk
	PixImpact* impact_to_save=grade_proxy->impact;

	// If this is the first event on the pixel, we initialise first
	if(grade_proxy->times==NULL){
		assert(next_impact!=NULL); // this function can only be called with NULL next_impact if there are saved impacts
		grade_proxy->times = (gradingTimeStruct*) malloc (sizeof(gradingTimeStruct));
		CHECK_NULL_VOID(grade_proxy->times,*status,"Malloc of first event failed");
		grade_proxy->times->previous = -1.0;
		grade_proxy->times->current = -1.0;
		grade_proxy->times->next = -1.0;
		grade_proxy->impact->pixID=next_impact->pixID;
		impact_to_save->time=-1.0;
	}

	//Very fist event, check if ctk has not triggered before
	if(grade_proxy->is_first==0 && is_crosstalk==0){
		//First we have to check if previous xtalk events do not trigger. If so, we store the triggered event and reenter the loop.
		computeAllCrosstalkInfluence(det,impact_to_save,grade_proxy->xtalk_proxy,grade_proxy,event_file,&(grade_proxy->crosstalk_energy),
								&(grade_proxy->nb_crosstalk_influence),next_impact->time,sample_length,&is_trigger,save_crosstalk,0,0,1,status);

		//If trigger there was, break and reprocess event with the pileup in grade proxy as current now.
		//Otherwise just copy impact normally
		if(is_trigger==0){
			grade_proxy->times->current = next_impact->time;
			grade_proxy->is_first=1; //We truly have our first event
			copyPixImpact(grade_proxy->impact,next_impact);
		}

	//If there is a trigger before first event (is_ctk==1), it is the first impact
	} else if(grade_proxy->is_first==0 && is_crosstalk==1){
		grade_proxy->times->current = next_impact->time;
		grade_proxy->is_first=1; //Ctk impact is the first
		copyPixImpact(grade_proxy->impact,next_impact);

	// We already have an other impact, we grade the previous and apply crosstalk depending on its grade
	} else {

		// Update grade proxy
		if (next_impact!=NULL){
			assert(grade_proxy->times->current>=0);
			// very rare case that a real crosstalk event is produced before (should never happen anymore)
			if ( (next_impact->time < grade_proxy->times->current) ){
				// switch the time of the events
				grade_proxy->times->next = grade_proxy->times->current;
				grade_proxy->times->current = next_impact->time;

				// also switch the impacts
				PixImpact* tmp_impact = next_impact;
				next_impact=impact_to_save;
				impact_to_save=tmp_impact;
			} else {
				grade_proxy->times->next = next_impact->time;
			}
		} else {
			grade_proxy->times->next = -1.0;
		}

		// Calculate grades
		calcGradingTimes(sample_length,*(grade_proxy->times),&grade1,&grade2,status);
		CHECK_STATUS_VOID(*status);
		// treat pileup case
		// TODO : add pileup limit as a pixel parameter (for the moment, equal to one sample)
		if(grade1==0){
			assert(next_impact!=NULL); // if next_impact is NULL, grade1 should be DEFAULTGOODSAMPLE
			next_impact->energy+=grade_proxy->impact->energy;
			impact_to_save->energy=0.;
			grade1=-1;
			grading=PILEUP;
			// do not change grading times in this case as this is a pile-up (this event will not be saved in the end)
			// We will apply the crosstalk to the total event at next time step
		} else {

			// We compute the grading of the event
			grading_index = makeGrading(grade1,grade2,&(det->pix[impact_to_save->pixID]));

			if (grading_index>=0){ //0 highest grade

				rmf = det->pix[impact_to_save->pixID].grades[grading_index].rmf;
				grading = det->pix[impact_to_save->pixID].grades[grading_index].value;

				// Get crosstalk influence from crosstalks closer to current time than max backward influence (others should have been dealt with) as function of the grading of the current pulse
				// -1 because grade.value returns starts at 1.

				computeAllCrosstalkInfluence(det,impact_to_save,grade_proxy->xtalk_proxy, grade_proxy,event_file,&(grade_proxy->crosstalk_energy),
						&(grade_proxy->nb_crosstalk_influence),grade_proxy->times->next,sample_length,&is_trigger,save_crosstalk,grading-1,0,1,status);

				//If there was a triggered event, we break from this loop, and reprocess the impact once more
				//with the new grade-proxy updated with the xt pileup (part 1)
				if (is_trigger==0){
					// Determine the measured detector channel (PI channel) according to the RMF.
					// The channel is obtained from the RMF using the corresponding HEAdas routine which is based on drawing a random number.
					returnRMFChannel(rmf, grade_proxy->impact->energy, &channel); //use total energy here to take pileup into account

					// Check if the photon is really measured. If the PI channel returned by the HEAdas RMF function is '-1', the photon is not
					// detected. This should not happen as the rmf is supposedly normalized
					if (channel<rmf->FirstChannel) {
						// flag as invalid (seemed better than discarding)
						if (!is_crosstalk && grade_proxy->impact->ph_id>=0){
							char msg[MAXMSG];
							sprintf(msg,"Impact found as not detected due to failure during RMF based energy allocation for energy %g keV (channel %ld) and phi_id %ld",
									grade_proxy->impact->energy,channel,grade_proxy->impact->ph_id);
							SIXT_WARNING(msg);
						}
						impact_to_save->energy=0.;
						grading=WRONGE;
					} else {

						// Determine the corresponding detected energy.
						// NOTE: In this simulation the collected charge is represented
						// by the nominal photon energy [keV], which corresponds to the
						// PI channel according to the EBOUNDS table.
						impact_to_save->energy=getEBOUNDSEnergy(channel, rmf, status);
						//printf("%f ",impact.energy);
						CHECK_STATUS_VOID(*status);
					}
				}
			} else { //Should the event be invalid for whichever reason, we do not treat it but flag it
				impact_to_save->energy=0.;
				grading=INVGRADE;
			}

			// Move ahead if it is worth it
			if (next_impact!=NULL && is_trigger==0){
				grade_proxy->times->previous = grade_proxy->times->current;
				grade_proxy->times->current = grade_proxy->times->next;
			}
		}

		//Should we have a crosstalk event changing the grading, flag it and go on (very rare!)
		if(is_trigger==0){
			if ((grdcmp != -99) && grdcmp!=grading-1){
				headas_chat(7," *** Grade num %i changed to %i because of crosstalk *** \n", grdcmp+1, grading);
				grading=GRADECHG; //The grade changed so we flag it
			}
			// Add processed event to event file if not pile-up!

			updateSignal(event_file,grade_proxy->row,impact_to_save->energy,grade1,grade2,grading,
					grade_proxy->nb_crosstalk_influence,grade_proxy->crosstalk_energy,status);
			if (*status!=EXIT_SUCCESS){
				SIXT_ERROR("updating signal energy in the event file failed");
				return;
			}

			// Finish moving ahead if worth it if no pile-up!
			if (next_impact!=NULL){
				grade_proxy->nb_crosstalk_influence = 0;
				grade_proxy->crosstalk_energy=0.;
				copyPixImpact(grade_proxy->impact,next_impact);
			}
		}
	}

	// Pre-add impact to event file to respect causality
	if (is_trigger==0 && next_impact!=NULL){
		grade_proxy->row = event_file->row;
		addEmptyEvent(event_file,next_impact,status);
	}

	//If a crosstalk event triggered, stop the process and once it was saved,
	//reprocess next_impact with this new event (which is now current)
	if(is_trigger==1){
		processGradedEvent(grade_proxy,sample_length,next_impact,det,event_file,0,save_crosstalk,-99,status);
	}
}
