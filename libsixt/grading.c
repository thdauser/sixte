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
		grade_proxys[ii].xtalk_proxy = NULL;
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
			applyMatrixEnerdepCrossTalk(det->pix[id].electrical_cross_talk,grade_proxys,sample_length,&impact,det,event_file,save_crosstalk,status);
			CHECK_STATUS_VOID(*status);
		}

		// intermod crosstalk
		if (det->crosstalk_intermod_table!=NULL){
			applyIntermodCrossTalk(grade_proxys,&impact,det,
					sample_length,event_file,save_crosstalk,status);
			CHECK_STATUS_VOID(*status);
		}

		// Process impact
		processGradedEvent(&(grade_proxys[id]),sample_length,&impact,det,event_file,0,status);
		CHECK_STATUS_VOID(*status);
	}

	// now we need to clean the remaining events
	for (int ii=0;ii<det->npix;ii++){
		if (grade_proxys[ii].times != NULL){
			int is_crosstalk=0;
			if (grade_proxys[ii].impact->ph_id<0){
				is_crosstalk=1;
			}
			processGradedEvent(&(grade_proxys[ii]),sample_length,NULL,det,event_file,is_crosstalk,status);
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

/** get the weight (i.e. additional energy here) for the intermodulation crosstalk */
static double get_intermod_weight(ImodTab* cross_talk, double df, double dt,
		double ampli_perturber, double ampli_signal, int* const status){

	double cross_talk_weight = 0.0;

	// calculate the amplitude for the given energy
	double ampl[]  = {ampli_perturber, ampli_signal};
	double d_ampl[2];
	int ind_ampl[] = {-1,-1};

	// get the amplitude bins in the weightening table
	for (int ii=0; ii<2 ; ii++){

		ind_ampl[ii] = binary_search(ampl[ii],cross_talk->ampl,cross_talk->n_ampl);
		if (ind_ampl[ii] < 0 ){
			headas_chat(3," *** warning: intermodulation cross talk amplitude %.3e not tabulated, skipping this event\n",ampl[ii]);
			return 0.0;
		}
		assert ( (ind_ampl[ii]) >= 0 && (ind_ampl[ii] < cross_talk->n_ampl) );

		d_ampl[ii] = (ampl[ii]-cross_talk->ampl[ind_ampl[ii]]) / ( cross_talk->ampl[ind_ampl[ii]+1] - cross_talk->ampl[ind_ampl[ii]] );
	}

	// get the time bin in the weightening table
	int ind_dt = binary_search(dt,cross_talk->dt,cross_talk->n_dt);

	assert ( (ind_dt) >= 0 && (ind_dt < cross_talk->n_dt) );
	double d_dt =  (dt - cross_talk->dt[ind_dt]) / ( cross_talk->dt[ind_dt+1] - cross_talk->dt[ind_dt] );

	// get the frequency bin in the weightening table
	int ind_df = binary_search(df,cross_talk->freq,cross_talk->n_freq);

	assert ( (ind_dt) >= 0 && (ind_dt < cross_talk->n_freq) );
	double d_df =  (df - cross_talk->freq[ind_df]) / ( cross_talk->freq[ind_df+1] - cross_talk->freq[ind_df] );

	int iarr[] = {ind_df, ind_dt,ind_ampl[0],ind_ampl[1]};
	double fac[] = {d_df, d_dt, d_ampl[0], d_ampl[1]};

	// get the crosstalk energy (in eV) and convert to keV
	double cross_energy = interp_lin_ndim(cross_talk->matrix,iarr,fac,4)*1e-3;

	// printf("  |-> weights: %.3e %.3e %.3e %.3e | %.3e %.3e | %.3e \n",c00,c01,c10,c11,c0,c1,cross_talk_weight);

	// check that the outcome has a reasonable value (less than 14 keV)
	if (fabs(cross_talk_weight) > 14  ){
		*status=EXIT_FAILURE;
		printf(" *** error: something went wrong when interpolating the intermodulation cross talk table ( weight = fabs(%.3e) > 28 keV) ... \n",
				cross_talk_weight);
		return 0.0;
	}

	return cross_energy;
}

// Conversion from Energy to Amplitude
// (Memo Roland, 10.06.2016 : 14keV = 0.25phi)
static double conv_ener2ampli(double ener){
	return ener / 14.0 * 0.5;
}


// taken from Roland's Memo (V2, 10.06.2016)
static double get_imod_df(double f_sig, double f_per){

	if (f_sig > f_per){
		return (f_sig-f_per)*(1.0e6/f_per) ; // correspond to 1MHz
	} else {
		return (f_sig-f_per)*(5.0e6/f_per) ; // correspond to 5MHz
	}

}

static void calc_imod_xt_influence(AdvDet* det,PixImpact* signal, PixImpact* perturber,
		PixImpact* crosstalk_impact, double dt, int* status){

	double df = get_imod_df(det->pix[signal->pixID].freq,det->pix[perturber->pixID].freq);
	double ampli_signal = conv_ener2ampli((signal->energy));
	double ampli_perturber = conv_ener2ampli((perturber->energy));


//	printf("      [ using df=%.1f dt=%e ampli_pert=%e ampli_sig=%e ]\n",
//			df,dt,ampli_perturber,ampli_signal);
	double energy_weight = get_intermod_weight(det->crosstalk_intermod_table, df, dt,
			ampli_perturber, ampli_signal, status);

    crosstalk_impact->energy = energy_weight; //  plus time-dep!!!
    crosstalk_impact->pixID = signal->pixID;             //  event is created in the singal-Pixel
    crosstalk_impact->time = perturber->time;            //  the event has the time of the perturber
    crosstalk_impact->ph_id = -perturber->ph_id;         //  negative ID to show it's a crosstalk event
    crosstalk_impact->src_id = perturber->src_id;


	headas_chat(7," -> energy %e eV ; dt=%e ; df=%.1f kHz \n",crosstalk_impact->energy*1e3,dt,df*1e-3);

}

/** Apply matrix cross talk: create new events on concerned pixels if corresponding
    energy is above the detection threshold, affect previous event otherwise */
void applyIntermodCrossTalk(GradeProxy* grade_proxys,PixImpact* impact, AdvDet* det,
		const double sample_length,TesEventFile* event_file,
		int save_crosstalk,int* const status){

    // This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
    PixImpact crosstalk_impact;
    crosstalk_impact.detposition.x = 0.;
    crosstalk_impact.detposition.y = 0.;
    crosstalk_impact.pixposition.x = 0.;
    crosstalk_impact.pixposition.y = 0.;

    PixImpact crosstalk_impact_ind;
    crosstalk_impact_ind.detposition.x = 0.;
    crosstalk_impact_ind.detposition.y = 0.;
    crosstalk_impact_ind.pixposition.x = 0.;
    crosstalk_impact_ind.pixposition.y = 0.;


	// find out the channel of pixel the photon hit
	Channel* active_chan = det->pix[impact->pixID].channel;

	// loop over all pixels in the channel to see if there is a hit
	for (int ii=0; ii < active_chan->num_pixels; ii++){
		int active_ind = active_chan->pixels[ii]->pindex;

		// do not count crosstalk with ourself!
		if (active_ind==impact->pixID) continue;

		// check if there has been an event in this pixel
		if (grade_proxys[active_ind].times!=NULL){

			double dt = (impact->time - grade_proxys[active_ind].times->current);

			// events should always be in time order (although crosstalk can create events earlier than the current time)
			assert(dt>=0);

			// see if we need to calculate the crosstalk influence of
			// the perturber on the signal
			if ( dt <= det->crosstalk_intermod_table->dt_max){

				printf(" *** [dt=%.2f sec] trigger crosstalk event in pix=%ld (%.2f MHz), with perturber %ld (%.2f MHz)\n",
						dt, grade_proxys[active_ind].impact->pixID,
						det->pix[active_ind].freq*1e-6,
						impact->pixID,det->pix[impact->pixID].freq*1e-6);


				calc_imod_xt_influence(det,grade_proxys[active_ind].impact,
						impact,&crosstalk_impact,dt,status);
				CHECK_STATUS_VOID(*status);

				// if signals are very close the signel might also influence the
 				//   perturber pulse
				//  Note:  (1) we require here that the above impact is written to the proxy already
				//         (2) dt -> -dt  in this case, as the order is reversed
				if ( -dt >= det->crosstalk_intermod_table->dt_min){
					dt = -dt;
					calc_imod_xt_influence(det,impact, grade_proxys[active_ind].impact,
							&crosstalk_impact_ind,dt,status);
					CHECK_STATUS_VOID(*status);


					 processCrosstalkEvent(&(grade_proxys[crosstalk_impact_ind.pixID]),sample_length,&crosstalk_impact_ind,
							 det,event_file,save_crosstalk,status);
					 CHECK_STATUS_VOID(*status);
				}

				// process events now (not before!!)
				processCrosstalkEvent(&(grade_proxys[crosstalk_impact.pixID]),sample_length,&crosstalk_impact,
						 det,event_file,save_crosstalk,status);
				CHECK_STATUS_VOID(*status);



				printf("  \n");

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

/** Same as applyMatricCrosstalk , but now the weights are energy dependent */
void applyMatrixEnerdepCrossTalk(MatrixEnerdepCrossTalk* cross_talk,GradeProxy* grade_proxys,const double sample_length,
		PixImpact* impact,AdvDet* det,TesEventFile* event_file,int save_crosstalk,int* const status){

	// This could be done in impactsToEvents, but seems less readable (it is just an information proxy, will be reused at each iteration)
	PixImpact crosstalk_impact;
	crosstalk_impact.detposition.x = 0.;
	crosstalk_impact.detposition.y = 0.;
	crosstalk_impact.pixposition.x = 0.;
	crosstalk_impact.pixposition.y = 0.;

	// make sure we have the same (length) energy array here
	assert(cross_talk->n_ener == det->crosstalk_elec_carrier_olap->n_ener_p);
	assert(cross_talk->n_ener == det->crosstalk_elec_common_imp->n_ener_p);
	double * ener_p = det->crosstalk_elec_carrier_olap->ener_p; // every electrical crosstalk energy vector has to be the same

	// Iterate over affected pixels
	for (int ii=0;ii<cross_talk->num_cross_talk_pixels;ii++){

		if ((impact->energy <= ener_p[0]) || (impact->energy >= ener_p[cross_talk->n_ener-1])) {
			printf(" *** warning : impact event energy %g outside the tabulated values for electrical crosstalk [%g,%g]}n",
					impact->energy,ener_p[0],ener_p[cross_talk->n_ener-1]);
			printf("     ---> skipping this event!\n");
		} else {

			// now determine the energy bin
			int ind = binary_search(impact->energy,det->crosstalk_elec_carrier_olap->ener_p,cross_talk->n_ener);

			double fac = (impact->energy - ener_p[ind]) / (ener_p[ind+1] - ener_p[ind]);

			crosstalk_impact.energy =
					    (1-fac)*cross_talk->cross_talk_weights[ii][ind] +
						(fac)*cross_talk->cross_talk_weights[ii][ind+1];
			crosstalk_impact.pixID = cross_talk->cross_talk_pixels[ii]->pindex;
			crosstalk_impact.time = impact->time;
			crosstalk_impact.ph_id = -impact->ph_id;
			crosstalk_impact.src_id = impact->src_id;

//			printf(" [%ld] crosstalk energy %.3f keV (for impact %.1f keV, fs=%.0f, fp=%.0f) \n",
//					crosstalk_impact.ph_id,crosstalk_impact.energy, impact->energy,
//					det->pix[impact->pixID].freq*1e-3,det->pix[crosstalk_impact.pixID].freq*1e-3);


			processCrosstalkEvent(&(grade_proxys[crosstalk_impact.pixID]),sample_length,&crosstalk_impact,
					det,event_file,save_crosstalk,status);
		}
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
	if(grade_proxy->xtalk_proxy==NULL){
		grade_proxy->xtalk_proxy = newCrosstalkProxy(status);
		CHECK_STATUS_VOID(*status);
	}

	// No crosstalk event available to possibly pileup with yet, copy info and go
	if(grade_proxy->xtalk_proxy->n_active_crosstalk==0){
		addCrosstalk2Proxy(grade_proxy->xtalk_proxy,impact,status);
	}
	else{
		// Should there be an assert(impact->time-grade_proxy->crosstalk->time>=0) ? I am worried about rounding errors but am not sure if I should be

		// If the two crosstalks happened in the same time bin, their combined energy could reach the threshold -> add energy and try to trigger
		PixImpact* previous_xtalk = grade_proxy->xtalk_proxy->xtalk_impacts[(grade_proxy->xtalk_proxy->current_crosstalk_index - 1)% grade_proxy->xtalk_proxy->xtalk_proxy_size];
		if(impact->time-previous_xtalk->time <= sample_length) { // TODO: code pileup length at pixel level (or use timedep table ?)
			previous_xtalk->energy+=impact->energy;
			previous_xtalk->nb_pileup++;// Update number of pileups
			if (previous_xtalk->energy >= det->threshold_event_lo_keV) {
				// Save crosstalk information
				double influence = previous_xtalk->energy;
				int nb_influence = previous_xtalk->nb_pileup;

				// Process now triggered event
				grade_proxy->xtalk_proxy->n_active_crosstalk--;// Crosstalk has been dealt with
				grade_proxy->xtalk_proxy->current_crosstalk_index--;
				previous_xtalk->time = impact->time;// to conserve causality
				processGradedEvent(grade_proxy,sample_length,previous_xtalk,det,event_file,1,status);

				// This new event comes from crosstalk and we should know it
				grade_proxy->nb_crosstalk_influence = nb_influence;
				grade_proxy->crosstalk_energy = influence;
			}
		// If they did not happen in the same time bin, treat crosstalk from events farther away than maximum backward influence (they will not
		// affect other events) and add new one
		} else {
			computeAllCrosstalkInfluence(det,grade_proxy->impact,grade_proxy->xtalk_proxy,&(grade_proxy->crosstalk_energy),
					&(grade_proxy->nb_crosstalk_influence),impact->time,grade_proxy->times==NULL,0);
			addCrosstalk2Proxy(grade_proxy->xtalk_proxy,impact,status);
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
	} else {
		// Add assert that at this stage, any active crosstalk is between the two triggered impacts ? NOT TRUE ANY MORE

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
			impact_to_save->energy=0.;
			grade1=-1;
			grading=PILEUP;
			// do not change grading times in this case as this is a pileup (this event will no be saved in the end)
		} else {

//			printf(" energy : %e (%e) [%i]\n",impact_to_save->energy,grade_proxy->impact->energy,is_crosstalk);

			// Get crosstalk influence from crosstalks closer to current time than max backward influence (others should have been dealt with)
			computeAllCrosstalkInfluence(det,impact_to_save,grade_proxy->xtalk_proxy,&(grade_proxy->crosstalk_energy),
					&(grade_proxy->nb_crosstalk_influence),grade_proxy->times->next,0,1);

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
		updateSignal(event_file,grade_proxy->row,impact_to_save->energy,grade1,grade2,grading,
				grade_proxy->nb_crosstalk_influence,grade_proxy->crosstalk_energy,status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("updating singal energy in the event file failed");
			return;
		}

		// Finish moving ahead if worth it
		if (next_impact!=NULL){
			grade_proxy->nb_crosstalk_influence = 0;
			grade_proxy->crosstalk_energy=0.;
			copyPixImpact(grade_proxy->impact,next_impact);
		}
	}

	// Pre-add impact to event file to respect causality
	if (next_impact!=NULL){
		grade_proxy->row = event_file->row;
		addEmptyEvent(event_file,next_impact,status);
	}
}
