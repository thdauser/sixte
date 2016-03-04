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

	//printf("grade1: %i grade2 %i => %i\n",*grade1,*grade2,makeGrading(*grade1,*grade2));

	if ( (*grade1<0) || (*grade2<0) ){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Grading of impacts failed in calcGradingTimes");
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
void impactsToEvents(AdvDet *det,PixImpFile *piximpactfile,TesEventFile* event_file,int save_crosstalk,int* const status){

	const double sample_length = 1./(det->SampleFreq);

	PixImpact impact;

	GradeProxy grade_proxys[det->npix];
	for (int ii=0;ii<det->npix;ii++){
		grade_proxys[ii].impact = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_VOID_STATUS(grade_proxys[ii].impact,*status);
		grade_proxys[ii].impact->energy=0.;
		grade_proxys[ii].crosstalk = NULL;
		grade_proxys[ii].n_active_crosstalk=0;
		grade_proxys[ii].times = NULL;
		grade_proxys[ii].row = event_file->row;
		grade_proxys[ii].crosstalk_energy=0.;
		grade_proxys[ii].nb_crosstalk_influence=0;
	}

	int id = -1;
	// Iterate over impacts
	while (getNextImpactFromPixImpFile(piximpactfile,&impact,status)){
		id = impact.pixID;
		// thermal crosstalk
		if (det->pix[id].thermal_cross_talk !=NULL){
			applyMatrixCrossTalk(det->pix[id].thermal_cross_talk,grade_proxys,sample_length,&impact,det,event_file,save_crosstalk,status);
		}
		// electrical crosstalk
		if (det->pix[id].electrical_cross_talk !=NULL){
			applyMatrixCrossTalk(det->pix[id].electrical_cross_talk,grade_proxys,sample_length,&impact,det,event_file,save_crosstalk,status);
		}

		// intermod crosstalk
		if (det->pix[id].intermodulation_cross_talk !=NULL){
			applyIntermodCrossTalk(det->pix[id].intermodulation_cross_talk,grade_proxys,sample_length,&impact,det,event_file,save_crosstalk,status);
		}

		// Process impact
		processGradedEvent(&(grade_proxys[id]),sample_length,&impact,det,event_file,0,status);
	}

	// now we need to clean the remaining events
	for (int ii=0;ii<det->npix;ii++){
		if (grade_proxys[ii].times != NULL){
			processGradedEvent(&(grade_proxys[ii]),sample_length,NULL,det,event_file,0,status);
			free(grade_proxys[ii].times);
		}
		free(grade_proxys[ii].impact);
		free(grade_proxys[ii].crosstalk);
	}
}


/** get the weight for the intermodulation crosstalk */
static double get_intermod_weight(IntermodulationCrossTalk* cross_talk,GradeProxy* grade_proxys,
		PixImpact* impact,int* const status){

	return -1.0;
}

/** Apply matrix cross talk: create new events on concerned pixels if corresponding energy is above the detection threshold, affect previous event otherwise */
void applyIntermodCrossTalk(IntermodulationCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){

	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;

	// Iterate over affected pixels
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){

		int num_comb = cross_talk->num_pixel_combinations[ii];

		for (int jj=0;jj<num_comb;jj++){
			double crosstalk_weight = get_intermod_weight(cross_talk,grade_proxys,impact,status);
			if (crosstalk_weight > 0){
				crosstalk_impact.energy = impact->energy*crosstalk_weight; // need to do this differently !! (depends on ampl1 and ampl2 and dt)
				crosstalk_impact.pixID = cross_talk->cross_talk_pixels[ii][jj]->pindex;
				crosstalk_impact.time = impact->time;
				crosstalk_impact.ph_id = -impact->ph_id;
				crosstalk_impact.src_id = impact->src_id;
				processCrosstalkEvent(&(grade_proxys[crosstalk_impact.pixID]),sample_length,&crosstalk_impact,det,event_file,save_crosstalk,status);
			}
		}
	}
}


/** Apply matrix cross talk: create new events on concerned pixels if corresponding energy is above the detection threshold, affect previous event otherwise */
void applyMatrixCrossTalk(MatrixCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){

	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;

	// Iterate over affected pixels
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){
		crosstalk_impact.energy = impact->energy*cross_talk->cross_talk_weights[ii];
		crosstalk_impact.pixID = cross_talk->cross_talk_pixels[ii]->pindex;
		crosstalk_impact.time = impact->time;
		crosstalk_impact.ph_id = -impact->ph_id;
		crosstalk_impact.src_id = impact->src_id;
		processCrosstalkEvent(&(grade_proxys[crosstalk_impact.pixID]),sample_length,&crosstalk_impact,det,event_file,save_crosstalk,status);
	}
}


/** Processes a crosstalk event using addCrosstalkEvent or processGradedEvent depending on whether it is above threshold or not */
void processCrosstalkEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){
	// the crosstalk event is above threshold and should be treated as a normal event
	if (impact->energy>=det->threshold_event_lo_keV) {
		processGradedEvent(grade_proxy,sample_length,impact,det,event_file,1,status);
	// the crosstalk event is below threshold and might influence another event
	} else {
		addCrosstalkEvent(grade_proxy,sample_length,impact,det,event_file,save_crosstalk,status);
	}
}

/** Adds a below threshold crosstalk event to the grading proxy */
void addCrosstalkEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){
	// save crosstalk event for debugging purposes ?
	if (save_crosstalk) addRMFImpact(event_file,impact,-2,-2,CROSSTALK,0,0.,status);

	assert(impact->energy<det->threshold_event_lo_keV);
	if (grade_proxy->crosstalk==NULL){
		grade_proxy->crosstalk = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_VOID_STATUS(grade_proxy->crosstalk,*status);
	}

	// No crosstalk event available to possibly pileup with yet, copy info and go
	if(grade_proxy->n_active_crosstalk==0){
		copyPixImpact(grade_proxy->crosstalk,impact);
		grade_proxy->n_active_crosstalk=1;
	}
	else{
		// Should there be an assert(impact->time-grade_proxy->crosstalk->time>=0) ? I am worried about rounding errors but am not sure if I should be

		// If the two crosstalks happened in the same time bin, their combined energy could reach the threshold -> add energy and try to trigger
		if(impact->time-grade_proxy->crosstalk->time <= sample_length) { // TODO: code pileup length at pixel level (or use timedep table ?)
			grade_proxy->crosstalk->energy+=impact->energy;
			grade_proxy->n_active_crosstalk++; // Update number of crosstalks
			if (grade_proxy->crosstalk->energy >= det->threshold_event_lo_keV) {
				// Save crosstalk information
				double influence = grade_proxy->crosstalk->energy;
				int nb_influence = grade_proxy->n_active_crosstalk;

				// Process now triggered event
				grade_proxy->n_active_crosstalk=0;// Crosstalk has been dealt with
				grade_proxy->crosstalk->time = impact->time;// to conserve causality
				processGradedEvent(grade_proxy,sample_length,grade_proxy->crosstalk,det,event_file,1,status);

				// This new event comes from crosstalk and we should know it
				grade_proxy->nb_crosstalk_influence = nb_influence;
				grade_proxy->crosstalk_energy = influence;
			}
		// If they did not happen in the same time bin, the previous one was indeed not detected and can
	    // safely be treated for possible influence on the previous triggered event
		} else {
			if(grade_proxy->times!=NULL){ //if there indeed was a previous impact
				int has_affected = computeCrosstalkInfluence(det,grade_proxy->impact,grade_proxy->crosstalk,&(grade_proxy->crosstalk_energy));
				if (has_affected){
					grade_proxy->nb_crosstalk_influence+=grade_proxy->n_active_crosstalk;
				}
			}
			// Now that the previous crosstalk is dealt with, only keep info of the new one
			copyPixImpact(grade_proxy->crosstalk,impact);
			grade_proxy->n_active_crosstalk=1;
		}
	}

}

/** Processes a graded event : update grading proxy and save previous event */
void processGradedEvent(GradeProxy* grade_proxy,const double sample_length,PixImpact* next_impact,AdvDet* det,TesEventFile* event_file,int is_crosstalk,int* const status){
	long grade1, grade2;
	int grading;
	int grading_index;
	struct RMF* rmf;
	long channel;
	PixImpact* impact_to_save=grade_proxy->impact;

	// If this is the first event on the pixel
	if(grade_proxy->times==NULL){
		assert(next_impact!=NULL); // this function can only be called with NULL next_impact if there are saved impacts
		grade_proxy->times = (gradingTimeStruct*) malloc (sizeof(gradingTimeStruct));
		CHECK_NULL_VOID(grade_proxy->times,*status,"malloc failed");

		grade_proxy->times->previous = -1.0;
		grade_proxy->times->current = next_impact->time;
		grade_proxy->times->next = -1.0;
		copyPixImpact(grade_proxy->impact,next_impact);
		// If by chance, there is a non-treated below threshold crosstalk event that happened just before
		// and should pileup, treat this case. This in principle can only happen on the first triggered event
		if (grade_proxy->n_active_crosstalk > 0){
			if (next_impact->time-grade_proxy->crosstalk->time<sample_length){ // TODO: code pileup length at pixel level (or use timedep table ?)
				grade_proxy->impact->energy+=grade_proxy->crosstalk->energy;
				grade_proxy->crosstalk_energy+=grade_proxy->crosstalk->energy;
				grade_proxy->nb_crosstalk_influence+=grade_proxy->n_active_crosstalk;
				grade_proxy->n_active_crosstalk = 0; // crosstalk is dealt with
			}
		}
	} else {
		// Add assert that at this stage, any active crosstalk is between the two triggered impacts ?

		// Update grade proxy
		if (next_impact!=NULL){
			grade_proxy->times->next = next_impact->time;
		} else {
			grade_proxy->times->next = -1.0;
		}

		// Calculate grades
		calcGradingTimes(sample_length,*(grade_proxy->times),&grade1,&grade2,status);

		// treat pileup case
		// TODO : add pileup limit as a pixel parameter (for the moment, equal to one sample)
		if(grade1==0){
			assert(next_impact!=NULL); // if next_impact is NULL, grade1 should be DEFAULTGOODSAMPLE
			next_impact->energy+=grade_proxy->impact->energy;
			// If there is a non triggered crosstalk between the two pileup events
			if (grade_proxy->n_active_crosstalk>0){
				next_impact->energy+=grade_proxy->crosstalk->energy;
				grade_proxy->crosstalk_energy+=grade_proxy->crosstalk->energy;
				grade_proxy->nb_crosstalk_influence+=grade_proxy->n_active_crosstalk;
				grade_proxy->n_active_crosstalk=0; // crosstalk is dealt with
			}
			impact_to_save->energy=0.;
			grade1=-1;
			grading=PILEUP;
			// do not change grading times in this case as this is a pileup (this event will no be saved in the end)
		} else {
			// Check whether there is some non-triggered crosstalk to take into account
			if (grade_proxy->n_active_crosstalk>0){
				int has_affected = computeCrosstalkInfluence(det,impact_to_save,grade_proxy->crosstalk,&(grade_proxy->crosstalk_energy));
				if (has_affected){
					grade_proxy->nb_crosstalk_influence+=grade_proxy->n_active_crosstalk;
				}
				grade_proxy->n_active_crosstalk=0;
			}

			grading_index = makeGrading(grade1,grade2,&(det->pix[impact_to_save->pixID]));

			if (grading_index>=0){

				rmf = det->pix[impact_to_save->pixID].grades[grading_index].rmf;
				grading = det->pix[impact_to_save->pixID].grades[grading_index].value;

				// Determine the measured detector channel (PI channel) according
				// to the RMF.
				// The channel is obtained from the RMF using the corresponding
				// HEAdas routine which is based on drawing a random number.
				returnRMFChannel(rmf, grade_proxy->impact->energy, &channel); //use total energy here to take pileup into account
				//printf("%f\n",grade_proxy->impact->energy);

				// Check if the photon is really measured. If the PI channel
				// returned by the HEAdas RMF function is '-1', the photon is not
				// detected. This should not happen as the rmf is supposedly
				// normalized
				if (channel<rmf->FirstChannel) {
					// flag as invalid (seemed better than discarding)
					if (!is_crosstalk){
						char msg[MAXMSG];
						sprintf(msg,"Impact found as not detected due to failure during RMF based energy allocation for energy %g keV and phi_id %ld",
								grade_proxy->impact->energy,grade_proxy->impact->ph_id);
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
			} else {
				impact_to_save->energy=0.;
				grading=INVGRADE;
			}

			// Move ahead if it is worth it
			if (next_impact!=NULL){
				grade_proxy->times->previous = grade_proxy->times->current;
				grade_proxy->times->current = grade_proxy->times->next;
			}
		}
		// Add processed event to event file
		updateSignal(event_file,grade_proxy->row,impact_to_save->energy,grade1,grade2,grading,grade_proxy->nb_crosstalk_influence,grade_proxy->crosstalk_energy,status);
		CHECK_STATUS_VOID(*status);

		// Finish moving ahead if worth it
		if (next_impact!=NULL){
			if (grading!=PILEUP){
				grade_proxy->nb_crosstalk_influence = 0;
				grade_proxy->crosstalk_energy=0.;
			}
			copyPixImpact(grade_proxy->impact,next_impact);
		}
	}

	// Pre-add impact to event file to respect causality
	if (next_impact!=NULL){
		grade_proxy->row = event_file->row;
		addEmptyEvent(event_file,next_impact,status);
	}
}

