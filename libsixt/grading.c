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


   Copyright 2016 Philippe Peille, IRAP; Thomas Dauser, ECAP
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


static imodProxy* newImodProxy(int* status){
	imodProxy* ipro = (imodProxy*) malloc ( sizeof(imodProxy)  );
	CHECK_MALLOC_RET_STATUS(ipro,NULL,*status);
	ipro->pixelHits = NULL;
	ipro->num_pix = 0;
	return ipro;
}
static void freeImodProxy(imodProxy* ipro){
	if (ipro != NULL){
		free(ipro->pixelHits);
	}
	free(ipro);
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
			CHECK_STATUS_VOID(*status);
		}
		// electrical crosstalk
		if (det->pix[id].electrical_cross_talk !=NULL){
			applyMatrixCrossTalk(det->pix[id].electrical_cross_talk,grade_proxys,sample_length,&impact,det,event_file,save_crosstalk,status);
			CHECK_STATUS_VOID(*status);
		}

		// intermod crosstalk
		if (det->pix[id].intermodulation_cross_talk !=NULL){
//			applyIntermodCrossTalk(det->pix[id].intermodulation_cross_talk,grade_proxys,sample_length,
//					&impact,det,event_file,save_crosstalk,status);
			CHECK_STATUS_VOID(*status);
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


static double get_ampli(double ener){
	return ener / 14.0 * 0.25;  // 14 keV is an amplitute of 0.25
}

// reverse function of get_ampli
static double get_imod_xt_energy(double ampli){
	return ampli * 14.0 / 0.25;
}

/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
static int binary_search(double val, double* arr, int n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
static int inv_binary_search(double val, double* arr, int n){
	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] > val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}

/** get the weight for the intermodulation crosstalk */
/**static double get_intermod_weight(ImodTab* cross_talk,	PixImpact* imp1,
		PixImpact* imp2, double dt, int* const status){

	double cross_talk_weight = 0.0;

	// calculate the amplitude for the given energy
	double ampl[]  = {get_ampli(imp1->energy),get_ampli(imp2->energy)};
	double d_ampl[2];
	int ind_ampl[] = {-1,-1};

	// get the amplitude bins in the weightening table
	for (int ii=0; ii<2 ; ii++){

		ind_ampl[ii] = inv_binary_search(ampl[ii],cross_talk->ampl,cross_talk->n_ampl);
		if (ind_ampl[ii] < 0 ){
			headas_chat(3," *** warning: intermodulation cross talk amplitude %.3e not tabulated, skipping this event\n",ampl[ii]);
			return 0.0;
		}
		assert ( (ind_ampl[ii]) >= 0 && (ind_ampl[ii] < cross_talk->n_ampl) );

		// be careful: amplitude is defined desceding
		d_ampl[ii] = (ampl[ii]-cross_talk->ampl[ind_ampl[ii]]) / ( cross_talk->ampl[ind_ampl[ii]+1] - cross_talk->ampl[ind_ampl[ii]] );
	}

	// get the time bin in the weightening table
	int ind_dt = binary_search(dt,cross_talk->dt,cross_talk->n_dt);

	assert ( (ind_dt) >= 0 && (ind_dt < cross_talk->n_dt) );
	double d_dt =  (dt - cross_talk->dt[ind_dt]) / ( cross_talk->dt[ind_dt+1] - cross_talk->dt[ind_dt] );

	// printf(" --> %i (%.3e) %i (%.3e) %i (%.3e) \n",ind_ampl[0],d_ampl[0],ind_ampl[1],d_ampl[1],ind_dt,d_dt);

	// interpolate in 3d: start with interpolation in dt
	double c00 = cross_talk->matrix[ ind_dt ][ind_ampl[0]][ind_ampl[1]] * (1-d_dt) +
				 cross_talk->matrix[ind_dt+1][ind_ampl[0]][ind_ampl[1]] * d_dt;
	double c01 = cross_talk->matrix[ ind_dt ][ind_ampl[0]][ind_ampl[1]+1] * (1-d_dt) +
				 cross_talk->matrix[ind_dt+1][ind_ampl[0]][ind_ampl[1]+1] * d_dt;
	double c10 = cross_talk->matrix[ ind_dt ][ind_ampl[0]+1][ind_ampl[1]] * (1-d_dt) +
				 cross_talk->matrix[ind_dt+1][ind_ampl[0]+1][ind_ampl[1]] * d_dt;
	double c11 = cross_talk->matrix[ ind_dt ][ind_ampl[0]+1][ind_ampl[1]+1] * (1-d_dt) +
				 cross_talk->matrix[ind_dt+1][ind_ampl[0]+1][ind_ampl[1]+1] * d_dt;

	// now along the first amplitude ampl[0]
	double c0 = c00 * (1-d_ampl[0]) + c10 * d_ampl[0];
	double c1 = c01 * (1-d_ampl[0]) + c11 * d_ampl[0];

	// ... and finally the last step
	cross_talk_weight  = c0 * (1-d_ampl[1]) + c1 * d_ampl[1];

	// printf("  |-> weights: %.3e %.3e %.3e %.3e | %.3e %.3e | %.3e \n",c00,c01,c10,c11,c0,c1,cross_talk_weight);

	// check that the outcome has a reasonable value
	if (cross_talk_weight < 0 ){
		*status=EXIT_FAILURE;
		printf(" *** error: something went wrong when interpolating the intermodulation cross talk table ( weight %.3e < 0) ... \n",
				cross_talk_weight);
		return 0.0;
	}

	return cross_talk_weight;
} */


/** check if a pixel is already in the intermodulation proxy
 * (means that it produced intermodulation crosstalk before and has to be skipped */
/**static int isPixelInProxy(imodProxy* ipro, long pindex){

	// make sure that the proxy is not NULL
	assert(ipro!=NULL);

	for (int ii=0; ii<ipro->num_pix; ii++){
		if (ipro->pixelHits[ii] == pindex){
			return 1;
		}
	}
	return 0;
} */

/**static void doImodCrossTalk_iter(IntermodulationCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact,long pindex_parent, AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){


	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;

	// Iterate over affected pixels
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){

		int ind_pix = cross_talk->cross_talk_pixels[ii][0]->pindex;
		int ind_pix_xt = cross_talk->cross_talk_pixels[ii][1]->pindex;

		// make sure that in the parent crosstalk pixel no event is created
//		if (pindex_parent == ind_pix_xt ){
//			//isPixelInProxy(imod_proxy,ind_pix)){
//			printf(" ** pixel %i is already in proxy, skipping ... \n", ind_pix_xt+1);
//			continue;
//		}

		// is there already an impact in the pixel
		// double dt = -1.0;
		if (grade_proxys[ind_pix].times==NULL){
			continue;
		}

		double dt = (impact->time - grade_proxys[ind_pix].times->current);

		// use the last bin of the tabulated time table for tmax (todo: should this value be fixed?)
		double tmax = cross_talk->cross_talk_info[ii]->dt[cross_talk->cross_talk_info[ii]->n_dt-1];

		// imod crosstalk only if 0 <= dt <= tmax (currently > 0 to avoid secondary imod events!)
		if ( (dt > 0) && (dt <= tmax) ){

			// check that the time is indeed greater than 0

			// get the amplitude of the cross talk event
			double crosstalk_ampli = get_intermod_weight(cross_talk->cross_talk_info[ii],
					impact, grade_proxys[ind_pix].impact, dt , status);


			if (crosstalk_ampli > 0){
				crosstalk_impact.energy = get_imod_xt_energy(crosstalk_ampli);
				crosstalk_impact.pixID = ind_pix_xt;
				crosstalk_impact.time = impact->time;       // will have the time of the latest impact
				crosstalk_impact.ph_id = -impact->ph_id;
				crosstalk_impact.src_id = -1;   // todo: what should we set here? could be two different sources

				// we have produced intermodulation cross talk
				headas_chat(7,"t=%.4e: Impact in %i and previous event in %i spaced dt=%.3e produce XT-signal with E=%.3e (%.3e) in %i \n",
						crosstalk_impact.time,impact->pixID+1,ind_pix+1,dt,
						get_imod_xt_energy(crosstalk_ampli),crosstalk_ampli,ind_pix_xt+1);

				processCrosstalkEvent(&(grade_proxys[crosstalk_impact.pixID]),sample_length,&crosstalk_impact,det,event_file,save_crosstalk,status);

				// now that we created another event, we need to see if it produces crosstalk as well
//				doImodCrossTalk_iter(det->pix[crosstalk_impact.pixID].intermodulation_cross_talk,grade_proxys,
//						sample_length,crosstalk_impact,impact->pixID,det,event_file,save_crosstalk,status);

			}
		}
	}

} */


//static void calc_imod_xt_influence(GradeProxy* grade_proxys,PixImpact* signal, PixImpact* perturber, int* status){
//}

/** Apply matrix cross talk: create new events on concerned pixels if corresponding
    energy is above the detection threshold, affect previous event otherwise */
/*void applyIntermodCrossTalk(GradeProxy* grade_proxys,PixImpact* impact, AdvDet* det,
		const double sample_length,TesEventFile* event_file,
		int save_crosstalk,int* const status){


	// find out the channel of pixel the photon hit
	Channel* active_chan = det->pix[impact->pixID]->channel;

	// loop over all pixels in the channel to see if there is a hit
	for (int ii=0; ii < active_chan->num_pixels; ii++){
		if (grade_proxys[active_chan->pixels[ii]->pindex].times!=NULL){

			double active_ind = active_chan->pixels[ii]->pindex;
			int dt = (impact->time - grade_proxys[active_ind].times->current);

			assert(dt>=0);

			// see if we need to calculate the crosstalk influence of
			// the perturber on the signal
			if (dt < IMOD_XT_UPPER_TAU*sample_length){
				calc_imod_xt_influence(grade_proxys[active_ind],grade_proxys[active_ind].impact,
						impact,status);

				// if signals are very close the signel might also influence the
 				//   perturber pulse
				//  Note: (1) we require here that the above impact is written to the proxy already
				//         (2) only the absolute value is of importance at the current implementation
				//             (actually it would be dt<0 in this case)
				if (dt < IMOD_XT_LOWER_TAU*sample_length){
					int perturb_ind = impact->pixID;
					calc_imod_xt_influence(grade_proxys[perturb_ind],impact,
							grade_proxys[perturb_ind].impact,status);
				}

				// process events now (not before!!)

			}

		}
	}


}  */


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

