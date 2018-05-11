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

#include "crosstalk.h"

/** Adds a cross talk pixel to the matrix */
static void add_xt_pixel(MatrixCrossTalk* matrix,AdvPix* pixel,int xt_index, double xt_weigth,int* const status){
	CHECK_STATUS_VOID(*status);

	// Allocate matrix if necessary
	if(matrix==NULL){
		matrix = newMatrixCrossTalk(status);
		CHECK_MALLOC_VOID_STATUS(matrix,*status);
	}

	// Increase matrix size
	matrix->cross_talk_pixels = realloc(matrix->cross_talk_pixels,(matrix->num_cross_talk_pixels+1)*sizeof(*(matrix->cross_talk_pixels)));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_pixels,*status);
	matrix->cross_talk_weights = realloc(matrix->cross_talk_weights,(matrix->num_cross_talk_pixels+1)*sizeof(*(matrix->cross_talk_weights)));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_weights,*status);
	matrix->cross_talk_index = realloc(matrix->cross_talk_index,(matrix->num_cross_talk_pixels+1)*sizeof(*(matrix->cross_talk_index)));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_index,*status);

	// Affect new values
	matrix->cross_talk_pixels[matrix->num_cross_talk_pixels] = pixel;
	matrix->cross_talk_weights[matrix->num_cross_talk_pixels] = xt_weigth;
	matrix->cross_talk_index[matrix->num_cross_talk_pixels]= xt_index;

	// Now, we can say that the matrix is effectively bigger
	matrix->num_cross_talk_pixels++;
}

/** Adds a cross talk pixel to the matrix */
static void add_xt_enerdep_pixel(MatrixEnerdepCrossTalk* matrix,AdvPix* pixel,double* xt_weigth,int n_ener,int* const status){
	CHECK_STATUS_VOID(*status);

	// Allocate matrix if necessary
	if (matrix==NULL){
		matrix = newMatrixEnerdepCrossTalk(pixel->ngrades, status);
		(matrix)->n_ener = n_ener;
		CHECK_MALLOC_VOID_STATUS(matrix,*status);
	}

	if (matrix->n_ener != n_ener){
		printf(" *** error : number of energy bins is %i, but should be %i \n",
				(matrix)->n_ener,n_ener);
		SIXT_ERROR(" adding energy-dependent crosstalk weight failed ");
		*status=EXIT_FAILURE;
		return;
	}

	// Increase matrix size
	matrix->cross_talk_pixels = realloc(matrix->cross_talk_pixels,(matrix->num_cross_talk_pixels+1)*sizeof(*(matrix->cross_talk_pixels)));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_pixels,*status);
	matrix->cross_talk_weights = realloc(matrix->cross_talk_weights,(matrix->num_cross_talk_pixels+1)*sizeof(*(matrix->cross_talk_weights)));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_weights,*status);

	// allocate space for the array (n_ener bins)
	matrix->cross_talk_weights[matrix->num_cross_talk_pixels] =(double*) malloc( n_ener *sizeof(double));
	CHECK_MALLOC_VOID_STATUS(matrix->cross_talk_weights[matrix->num_cross_talk_pixels],*status);

	// Affect new values
	matrix->cross_talk_pixels[matrix->num_cross_talk_pixels] = pixel;
	for (int ii=0; ii<n_ener; ii++){
		matrix->cross_talk_weights[matrix->num_cross_talk_pixels][ii] = xt_weigth[ii];
	}

	// Now, we can say that the matrix is effectively bigger
	matrix->num_cross_talk_pixels++;
}

// taken from Roland's Memo (V2, 10.06.2016)
double get_imod_df(double f_sig, double f_per, int* status){

	double f_high = 5.6e6; // changed to 5.6MHz
	double f_low  = 1.0e6;

	// integrate check df_nlin < f_high-fg_low

	if (f_sig > f_high || f_sig < f_low){
		*status=EXIT_FAILURE;
		printf(" *** error: expecting frequency %gMHz to be in the interval [%g,%g] MHz",
				f_sig,f_low,f_high);
		SIXT_ERROR("simulating non-linear crosstalk failed");
		return 0;
	}

	if (f_per > f_high || f_per < f_low){
		*status=EXIT_FAILURE;
		printf(" *** error: expecting frequency %gMHz to be in the interval [%g,%g] MHz",
				f_per,f_low,f_high);
		SIXT_ERROR("simulating non-linear crosstalk failed");
		return 0;
	}

	if (f_sig > f_per){
		return (f_sig-f_per)*(f_low/f_per) ;
	} else {
		return (f_sig-f_per)*(f_high/f_per) ;
	}

}

/** get the weight (i.e. additional energy here) for the intermodulation crosstalk */
double get_intermod_weight(AdvDet* det, int grading, double df, double dt,
		double ampli_perturber, double ampli_signal, int* const status){

	double cross_talk_weight = 0.0;

	// calculate the amplitude for the given energy
	double ampl[]  = {ampli_perturber, ampli_signal};
	double d_ampl[2];
	int ind_ampl[] = {-1,-1};

	// get the amplitude bins in the weightening table
	for (int ii=0; ii<2 ; ii++){

		ind_ampl[ii] = binary_search(ampl[ii],det->crosstalk_imod_table[grading].ampl,det->crosstalk_imod_table[grading].n_ampl);
		if (ind_ampl[ii] < 0 ){
			headas_chat(3," *** warning: intermodulation cross talk amplitude %.3e not tabulated, skipping this event\n",ampl[ii]);
			return 0.0;
		}
		assert ((ind_ampl[ii])>= 0 && (ind_ampl[ii]<det->crosstalk_imod_table[grading].n_ampl));

		d_ampl[ii] = (ampl[ii]-det->crosstalk_imod_table[grading].ampl[ind_ampl[ii]]) / ( det->crosstalk_imod_table[grading].ampl[ind_ampl[ii]+1] - det->crosstalk_imod_table[grading].ampl[ind_ampl[ii]] );
	}

	// get the time bin in the weighttable

	int ind_dt = binary_search(dt,det->crosstalk_imod_table[grading].dt,det->crosstalk_imod_table[grading].n_dt);
	assert ((ind_dt)>= 0 && (ind_dt < det->crosstalk_imod_table[grading].n_dt) );
	double d_dt =  (dt - det->crosstalk_imod_table[grading].dt[ind_dt]) / ( det->crosstalk_imod_table[grading].dt[ind_dt+1] - det->crosstalk_imod_table[grading].dt[ind_dt] );

	// get the frequency bin in the weightening table
	int ind_df = binary_search(df,det->crosstalk_imod_table[grading].freq,det->crosstalk_imod_table[grading].n_freq);
	if ( ind_df < 0 ){
		if ( df <= det->crosstalk_imod_table[grading].freq[0]){
			ind_df = 0;
		} else if ( df >= det->crosstalk_imod_table[grading].freq[det->crosstalk_imod_table[grading].n_freq-1]){
			ind_df = det->crosstalk_imod_table[grading].n_freq-1;
		}
		headas_chat(7, " **** WARNING: non-linear crosstalk frequency df=%e not available, using df=%e instead!\n",df,det->crosstalk_imod_table[grading].freq[ind_df]);
	}

	assert( (ind_df) >= 0 && (ind_df < det->crosstalk_imod_table[grading].n_freq) );
	double d_df =  (df - det->crosstalk_imod_table[grading].freq[ind_df]) / ( det->crosstalk_imod_table[grading].freq[ind_df+1] - det->crosstalk_imod_table[grading].freq[ind_df] );

	int iarr[] = {ind_df, ind_dt,ind_ampl[0],ind_ampl[1]};
	double fac[] = {d_df, d_dt, d_ampl[0], d_ampl[1]};

	// get the crosstalk energy (in eV) and convert to keV
	double cross_energy = interp_lin_ndim(det->crosstalk_imod_table[grading].matrix,iarr,fac,4)*1e-3;

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
// eMail Roland, 22.06.2016 : Take absolute value for conversion (as AC modulated signal)
// new Version (eMail Roland 15.11.2016): 12keV = 0.15phi
double conv_ener2ampli(double ener){
	return fabs(ener) / 12.0 * 0.15;
}

void calc_imod_xt_influence(AdvDet* det,PixImpact* signal, PixImpact* perturber, double* enweight, double dt, int grading, int* const status){

	/** Checking if perturber frequency is higher or lower than pixels*/
	double df = get_imod_df(det->pix[signal->pixID].freq,det->pix[perturber->pixID].freq,status);
	CHECK_STATUS_VOID(*status);
	double ampli_signal = conv_ener2ampli((signal->energy));
	double ampli_perturber = conv_ener2ampli((perturber->energy));
	double energy_weight = get_intermod_weight(det,grading, df, dt, ampli_perturber, ampli_signal, status);

	headas_chat(7," -> time: %g ; energy %e eV ; dt=%e ; df=%.1f kHz \n",signal->time,
			signal->energy*1e3,dt,df*1e-3);
    *enweight=energy_weight; //  plus time-dep!!!
}

void calc_prop_xt_influence(AdvDet* det,double signal_energy, double perturber_energy, double* enweight, double dt_in_frames, int grading){
	//Check if within the sample range for ener_p and ener_v
	if ((perturber_energy < det->crosstalk_TDM_prop[grading].ener_p[0]) ||
			(perturber_energy >= det->crosstalk_TDM_prop[grading].ener_p[det->crosstalk_TDM_prop[grading].n_ener_p-1])||
			(signal_energy < det->crosstalk_TDM_prop[grading].ener_v[0]) ||
			(signal_energy >= det->crosstalk_TDM_prop[grading].ener_p[det->crosstalk_TDM_prop[grading].n_ener_v-1])) {
		headas_chat(7, " *** warning : impact event energy %g or perturber energy %g is outside the tabulated values for proportional cross-talk [%g,%g]}n",
				signal_energy, perturber_energy, det->crosstalk_TDM_prop[grading].ener_p[0],det->crosstalk_TDM_prop[grading].ener_p[det->crosstalk_TDM_prop[grading].n_ener_p-1]);
		headas_chat(7, "     ---> skipping this event!\n");
	//Check if within the sample range in time
	} else if ((dt_in_frames<= det->crosstalk_TDM_prop[grading].samples[0]) ||
			(dt_in_frames >= det->crosstalk_TDM_prop[grading].samples[det->crosstalk_TDM_prop[grading].n_samples-1])){
		//The impact is not affected
		headas_chat(7, "  Event outside time frames\n");
		headas_chat(7, "     ---> skipping this event!\n");
	} else {
		//Finding the corresponding index for 0
		int ind_dt = binary_search(dt_in_frames,det->crosstalk_TDM_prop[grading].samples,det->crosstalk_TDM_prop[grading].n_samples);
		//printf("dt is %f and index is %i\n", dt_in_frames, ind_dt);
		assert ((ind_dt)>= 0 && (ind_dt < det->crosstalk_TDM_prop[grading].n_samples));
		float d_dt = (dt_in_frames-det->crosstalk_TDM_prop[grading].samples[ind_dt])/(det->crosstalk_TDM_prop[grading].samples[ind_dt+1]-det->crosstalk_TDM_prop[grading].samples[ind_dt]);

		//Perturber energy
		int ind_dep = binary_search(perturber_energy,det->crosstalk_TDM_prop[grading].ener_p,det->crosstalk_TDM_prop[grading].n_ener_p);
		assert ((ind_dep)>= 0 && (ind_dep < det->crosstalk_TDM_prop[grading].n_ener_p));
		//printf("energy is %f and index is %i\n", perturber_energy, ind_dep);
		float d_dep=(perturber_energy-det->crosstalk_TDM_prop[grading].ener_p[ind_dep])/(det->crosstalk_TDM_prop[grading].ener_p[ind_dep+1]-det->crosstalk_TDM_prop[grading].ener_p[ind_dep]);

		//Victim energy
		int ind_dev = binary_search(signal_energy,det->crosstalk_TDM_prop[grading].ener_v,det->crosstalk_TDM_prop[grading].n_ener_v);
		//printf("energy is %f and index is %i\n", signal_energy, ind_dev);
		assert ((ind_dev)>= 0 && (ind_dev < det->crosstalk_TDM_prop[grading].n_ener_v));
		float d_dev=(signal_energy-det->crosstalk_TDM_prop[grading].ener_v[ind_dev])/(det->crosstalk_TDM_prop[grading].ener_v[ind_dev+1]-det->crosstalk_TDM_prop[grading].ener_v[ind_dev]);

		int iarr[] = {ind_dt, ind_dep,ind_dev};
		double fac[] = {d_dt, d_dep, d_dev};

		// get the perturber energy (in keV) and rescale it
		*enweight = interp_lin_ndim(det->crosstalk_TDM_prop[grading].matrix,iarr,fac,3);
		//printf("In routine %f\n", *enweight);
	}
}

void calc_der_xt_influence(AdvDet* det,double signal_energy, double perturber_energy, double* enweight, double dt_in_frames, int grading){
	//Check if within the sample range for ener_p and ener_v (same)
	if ((perturber_energy < det->crosstalk_TDM_der[grading].ener_p[0]) ||
			(perturber_energy >= det->crosstalk_TDM_der[grading].ener_p[det->crosstalk_TDM_der[grading].n_ener_p-1]) ||
			(signal_energy < det->crosstalk_TDM_der[grading].ener_v[0]) ||
			(signal_energy >= det->crosstalk_TDM_der[grading].ener_p[det->crosstalk_TDM_der[grading].n_ener_v-1])) {
		headas_chat(7, " *** warning : energy outside the tabulated values for proportional crosstalk [%g,%g]}n",
				det->crosstalk_TDM_der[grading].ener_p[0],det->crosstalk_TDM_der[grading].ener_p[det->crosstalk_TDM_der[grading].n_ener_p-1]);
		headas_chat(7, "     ---> skipping this event!\n");
		//printf("Outside energy cross %f and impact %f\n", perturber_energy, signal_energy);
		//printf("Boundaries %f and %f\n", det->crosstalk_TDM_der[grading].ener_v[0],  det->crosstalk_TDM_der[grading].ener_p[det->crosstalk_TDM_der[grading].n_ener_v-1]);
	//Check if within the sample range in time
	} else if ((dt_in_frames<=det->crosstalk_TDM_der[grading].samples[0]) ||
			(dt_in_frames >= det->crosstalk_TDM_der[grading].samples[det->crosstalk_TDM_der[grading].n_samples-1])){
		//The impact is not affected
		headas_chat(7, "  Not influencing!\n");
		headas_chat(7, "     ---> skipping this event!\n");
		//printf("Outside time\n");
		//printf("Boundaries %f and %f\n", det->crosstalk_TDM_der[grading].samples[0],  det->crosstalk_TDM_der[grading].samples[det->crosstalk_TDM_der[grading].n_samples-1]);
	} else{
		//Finding the corresponding indices
		int ind_dt = binary_search(dt_in_frames,det->crosstalk_TDM_der[grading].samples,det->crosstalk_TDM_der[grading].n_samples);
		assert ((ind_dt)>= 0 && (ind_dt < det->crosstalk_TDM_der[grading].n_samples));
		//printf("dt is %f and index is %i\n", dt_in_frames, ind_dt);
		float d_dt = (dt_in_frames-det->crosstalk_TDM_der[grading].samples[ind_dt])/(det->crosstalk_TDM_der[grading].samples[ind_dt+1]-det->crosstalk_TDM_der[grading].samples[ind_dt]);

		int ind_dep = binary_search(perturber_energy,det->crosstalk_TDM_der[grading].ener_p,det->crosstalk_TDM_der[grading].n_ener_p);
		assert ((ind_dep)>= 0 && (ind_dep < det->crosstalk_TDM_der[grading].n_ener_p));
		//printf("energy is %f and index is %i\n", perturber_energy, ind_dep);
		float d_dep=(perturber_energy-det->crosstalk_TDM_der[grading].ener_p[ind_dep])/(det->crosstalk_TDM_der[grading].ener_p[ind_dep+1]-det->crosstalk_TDM_der[grading].ener_p[ind_dep]);

		int ind_dev = binary_search(signal_energy,det->crosstalk_TDM_der[grading].ener_v,det->crosstalk_TDM_der[grading].n_ener_v);
		assert ((ind_dev)>= 0 && (ind_dev < det->crosstalk_TDM_der[grading].n_ener_v));
		//printf("energy is %f and index is %i\n", signal_energy, ind_dev);
		float d_dev=(signal_energy-det->crosstalk_TDM_der[grading].ener_v[ind_dev])/(det->crosstalk_TDM_der[grading].ener_v[ind_dev+1]-det->crosstalk_TDM_der[grading].ener_v[ind_dev]);

		int iarr[] = {ind_dt, ind_dep,ind_dev};
		double fac[] = {d_dt, d_dep, d_dev};

		// get the crosstalk energy (in keV) and rescale it
		*enweight = interp_lin_ndim(det->crosstalk_TDM_der[grading].matrix,iarr,fac,3);
		//printf("In routine %f\n", *enweight);
	}
}


//Find which file we need to use for time dependency
CrosstalkTimedep* getTimeDep(AdvDet* det, CrosstalkProxy* xtalk_proxy, int ii, int grade, int* status){
	CrosstalkTimedep* buffer=NULL;
	PixImpact* crosstalk=xtalk_proxy->xtalk_impacts[ii];
	if(xtalk_proxy->type[ii]==-ELECCTK){
		buffer=&(det->crosstalk_elec_timedep[2*grade]);
	} else if(xtalk_proxy->type[ii]==ELECCTK){
		buffer=&(det->crosstalk_elec_timedep[2*grade+1]);
	} else if(xtalk_proxy->type[ii]==-THERCTK){
		buffer=&(det->crosstalk_ther_timedep[det->pix[crosstalk->pixID].thermal_cross_talk->cross_talk_index[crosstalk->weight_index]][2*grade]); //In case the timedep is different beteween grades
	} else if(xtalk_proxy->type[ii]==THERCTK){
		buffer=&(det->crosstalk_ther_timedep[det->pix[crosstalk->pixID].thermal_cross_talk->cross_talk_index[crosstalk->weight_index]][2*grade+1]);
	} else{
		printf("It seems as though there is an event of cross-talk type %i which has no timedependency given in XML file \n", xtalk_proxy->type[ii]);
		SIXT_ERROR("Wrong type");
		*status=EXIT_FAILURE;
	}
	return buffer;
}

/** Routine to empty a given cross-talk mechanism from proxy given their indices */
static void erasectk(CrosstalkProxy* xtalk_proxy, int* toerase, int erased_crosstalks, int* const status){
	assert(erased_crosstalks<=xtalk_proxy->n_active_crosstalk); //Not more to erase than current
	assert(erased_crosstalks>0); //Actually some to erase

	if (erased_crosstalks<xtalk_proxy->n_active_crosstalk){
		long int c=0;
		long int d=0;

		//Erasing amongst the actives
		for(int k=0;k<xtalk_proxy->n_active_crosstalk;k++){
			if(c<erased_crosstalks && k==toerase[c]){ //Useless to access indices higher than c itself
				c+=1;
			} else{
				copyPixImpact(xtalk_proxy->xtalk_impacts[d],xtalk_proxy->xtalk_impacts[k]); //Just shift the index
				xtalk_proxy->type[d]=xtalk_proxy->type[k];
				xtalk_proxy->is_saved[d]=xtalk_proxy->is_saved[k];
				d+=1;
			}
		}

		assert(c==erased_crosstalks); //Check if we have erased them all indeed

		xtalk_proxy->type=(int*) realloc(xtalk_proxy->type, (xtalk_proxy->xtalk_proxy_size-erased_crosstalks)*sizeof(int));
		xtalk_proxy->is_saved=(int*) realloc(xtalk_proxy->is_saved, (xtalk_proxy->xtalk_proxy_size-erased_crosstalks)*sizeof(int));
		xtalk_proxy->n_active_crosstalk-=erased_crosstalks;
	} else{
		//Freeing everything in the structure
		for(int k=0;k<xtalk_proxy->xtalk_proxy_size;k++){
			free(xtalk_proxy->xtalk_impacts[k]);
		}

		PixImpact** tmpimp = (PixImpact**) malloc(INITXTALKNB*sizeof(PixImpact*));
		CHECK_MALLOC_VOID_STATUS(tmpimp,*status);
		for (int ii=0;ii<INITXTALKNB;ii++){
			tmpimp[ii] = (PixImpact*) malloc(sizeof(PixImpact));
			CHECK_MALLOC_VOID_STATUS(tmpimp[ii],*status);
		}

		free(xtalk_proxy->type);
		free(xtalk_proxy->is_saved);
		free(xtalk_proxy->xtalk_impacts);
		xtalk_proxy->type=NULL;
		xtalk_proxy->is_saved=NULL;
		xtalk_proxy->xtalk_impacts=tmpimp;
		xtalk_proxy->n_active_crosstalk=0;
		xtalk_proxy->xtalk_proxy_size=INITXTALKNB;
	}
}

void get_imodtable_axis(int* nrows, double** val, char* extname, char* colname, fitsfile* fptr, int* status){
	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return;
	}

	// get the column id-number
	int colnum;
	if(fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

	// get the number of rows
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return;

	// allocate memory for the array
	*val=(double*)malloc(n*sizeof(double));
	CHECK_MALLOC_VOID_STATUS(*val,*status);

    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) n;
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);

	(*nrows) = (int) n;

	return;
}

// Loads thermal cross-talk for requested pixel
// Concretely, iterates over all the pixels to find neighbours
void load_thermal_cross_talk(AdvDet* det,int pixid,int* const status){
	double max_cross_talk_dist = det->xt_dist_thermal[det->xt_num_thermal-1];
	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;
	concerned_pixel->thermal_cross_talk=newMatrixCrossTalk(status);

	for (int i=0;i<det->npix;i++){
		// Cross-talk is not with yourself ;)
		if (i==pixid){
			continue;
		}
		current_pixel = &(det->pix[i]);

		// Initial quick distance check to avoid spending time on useless pixels
		if((fabs(current_pixel->sx-concerned_pixel->sx) > max_cross_talk_dist) ||
				(fabs(current_pixel->sy-concerned_pixel->sy) > max_cross_talk_dist)){
			continue;
		}

		// Get distance between two pixels
		//double pixel_distance = distance_two_pixels(current_pixel,concerned_pixel); //Deprecated way
		double pixel_distance = sqrt(pow(current_pixel->sx-concerned_pixel->sx,2)+pow(current_pixel->sy-concerned_pixel->sy,2));
		if (pixel_distance<0){ // distance should be positive
			*status = EXIT_FAILURE;
			printf("*** error: Distance between pixels %d and %d is negative\n",current_pixel->pindex,concerned_pixel->pindex);
			return;
		}

		// Iterate over cross-talk values and look for the first matching one
		for (int xt_index=0;xt_index<det->xt_num_thermal;xt_index++){
			if (pixel_distance<det->xt_dist_thermal[xt_index]){
				add_xt_pixel(concerned_pixel->thermal_cross_talk,current_pixel,xt_index,det->xt_weight_thermal[xt_index],status);
				CHECK_STATUS_VOID(*status);
				// If one cross talk was identified, go to next pixel (the future cross-talks should be lower order cases)
				break;
			}
		}
	}
}

void load_thermal_timedep(AdvDet* det, int i, int k, int* const status){
	CHECK_NULL_VOID(det->crosstalk_thermal_timedep_file[i],*status,"no file for the thermal time dependence specified");

	CrosstalkTimedep* tmp_lo=newCrossTalkTimedep(status);
	det->crosstalk_ther_timedep[i][2*k]=*tmp_lo;
	free(tmp_lo);
	CrosstalkTimedep* tmp_hi=newCrossTalkTimedep(status);
	det->crosstalk_ther_timedep[i][2*k+1]=*tmp_hi;
	free(tmp_hi);

	char* COLNAME_TIME_DELAY = "TIME_DELAY";
	char* COLNAME_WEIGHT = "WEIGHT";
	char* EXTNAME_LOW_GRAD;
	char* EXTNAME_HI_GRAD;
	// Concatenating the grade

	if (asprintf(&EXTNAME_LOW_GRAD,"PERTLO%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	if (asprintf(&EXTNAME_HI_GRAD,"PERTHI%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}


	fitsfile *fptr=NULL;

	do {
		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_thermal_timedep_file[i]);
		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the time dependence thermal table %s for grade %i \n",fullfilename, k);

		// Read time dep info for first file
		int n_time_low, n_weight_low;
		double* time_low;
		double* weight_low;
		get_imodtable_axis(&n_time_low,&time_low,EXTNAME_LOW_GRAD,COLNAME_TIME_DELAY,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_weight_low,&weight_low,EXTNAME_LOW_GRAD,COLNAME_WEIGHT,fptr,status);
		CHECK_STATUS_BREAK(*status);

		det->crosstalk_ther_timedep[i][2*k].length=n_time_low;
		det->crosstalk_ther_timedep[i][2*k].time=time_low;
		det->crosstalk_ther_timedep[i][2*k].weight=weight_low;

		// now we need to set the weight at t=0 (differs for tessim simulated LUT)
		int ind_t0_low = binary_search(0.0, det->crosstalk_ther_timedep[i][2*k].time, det->crosstalk_ther_timedep[i][2*k].length);

		double fac_low = (0.0 - det->crosstalk_ther_timedep[i][2*k].time[ind_t0_low]) /
				(det->crosstalk_ther_timedep[i][2*k].time[ind_t0_low+1] - det->crosstalk_ther_timedep[i][2*k].time[ind_t0_low]);

		det->crosstalk_ther_timedep[i][2*k].weight_t0 =
				(1-fac_low)*det->crosstalk_ther_timedep[i][2*k].weight[ind_t0_low] +
				(fac_low)  *det->crosstalk_ther_timedep[i][2*k].weight[ind_t0_low+1];

		// Read time dependence info for second the second case (same grading, different frequency difference)
		int n_time_hi, n_weight_hi;
		double* time_hi;
		double* weight_hi;
		get_imodtable_axis(&n_time_hi,&time_hi,EXTNAME_HI_GRAD,COLNAME_TIME_DELAY,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_weight_hi,&weight_hi,EXTNAME_HI_GRAD,COLNAME_WEIGHT,fptr,status);
		CHECK_STATUS_BREAK(*status);

		det->crosstalk_ther_timedep[i][2*k+1].length=n_time_hi;
		det->crosstalk_ther_timedep[i][2*k+1].time=time_hi;
		det->crosstalk_ther_timedep[i][2*k+1].weight=weight_hi;

		// now we need to set the weight at t=0 (differs for tessim simulated LUT)
		int ind_t0_hi = binary_search(0.0, det->crosstalk_ther_timedep[i][2*k+1].time, det->crosstalk_ther_timedep[i][2*k+1].length);

		double fac_hi = (0.0 - det->crosstalk_ther_timedep[i][2*k+1].time[ind_t0_hi]) /
				(det->crosstalk_ther_timedep[i][2*k+1].time[ind_t0_hi+1] - det->crosstalk_ther_timedep[i][2*k+1].time[ind_t0_hi]);

		det->crosstalk_ther_timedep[i][2*k+1].weight_t0 =
				(1-fac_hi)*det->crosstalk_ther_timedep[i][2*k+1].weight[ind_t0_hi] +
				(fac_hi)  *det->crosstalk_ther_timedep[i][2*k+1].weight[ind_t0_hi+1];

	}while(0);

	free(EXTNAME_LOW_GRAD);
	free(EXTNAME_HI_GRAD);

	if (fptr!=NULL) {fits_close_file(fptr,status);}
}

static void initElecTab(ElecTab** tab, int n_freq_s, int n_freq_p, int n_ener_p,
		double* freq_s, double* freq_p, double* ener_p, int* status){
	// make a short pointer
	ElecTab* t = (*tab);

	t->n_freq_s = n_freq_s;
	t->n_freq_p  = n_freq_p;
	t->n_ener_p = n_ener_p;

	t->freq_s = (double*) malloc(n_freq_s * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->freq_s,*status);
	for (int ii=0; ii<n_freq_s; ii++){
		t->freq_s[ii] = freq_s[ii];
	}

	t->freq_p = (double*) malloc(n_freq_p * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->freq_p,*status);
	for (int ii=0; ii<n_freq_p; ii++){
		t->freq_p[ii] = freq_p[ii];
	}

	t->ener_p = (double*) malloc(n_ener_p * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->ener_p,*status);
	for (int ii=0; ii<n_ener_p; ii++){
		t->ener_p[ii] = ener_p[ii];
	}

	// allocate the 3d matrix (n_freq_s x n_freq_p x n_ener_p)
	t->matrix = (double***) malloc (n_freq_s*sizeof(double**));
	CHECK_MALLOC_VOID_STATUS(t->matrix,*status);

	for (int ll=0; ll<n_freq_s; ll++){               // FREQ_S-LOOP
		t->matrix[ll] = (double**) malloc (n_freq_p*sizeof(double*));
		CHECK_MALLOC_VOID_STATUS(t->matrix[ll],*status);

		for (int ii=0; ii<n_freq_p; ii++){             // FREQ_P-LOOP
			t->matrix[ll][ii] = (double*) malloc (n_ener_p*sizeof(double));
			CHECK_MALLOC_VOID_STATUS(t->matrix[ll][ii],*status);

			for (int jj=0; jj<n_ener_p; jj++){      // ENER_P-LOOP
				t->matrix[ll][ii][jj] = 0.0;
			}
		}
	}
	return;
}

/** read one electrical crosstalk matrix from a fits file and store it in the structure (units given in keV) */
void read_elec_matrix(fitsfile* fptr,int n_freq_s, int n_freq_p, int n_ener_p, float scaling,
		ElecTab* tab, char* extname, int* status){

	// ampl and dt arrays should be initialized before
	assert(tab->freq_s!=NULL);
	assert(tab->freq_p!=NULL);
	assert(tab->ener_p!=NULL);

	// move to the extension containing the data
	int extver = 0;
	if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver ,status)){
		printf(" error: moving to extension %s failed \n",extname);
		return;
	}

	// get number of rows
	int naxis;
	if (fits_get_img_dim(fptr, &naxis, status) && (naxis!=4)){
		printf(" error: getting dimensions of electrical data array failed \n");
		return;
	}
	long naxes[naxis];


	if (fits_get_img_size(fptr, naxis, naxes, status) ){
		printf(" error: getting dimensions of electrical data array failed \n");
		return;
	}
	// check dimensions
	if (naxes[0]!=n_freq_s || naxes[1]!=n_freq_p || naxes[2]!=n_ener_p ){
		*status=EXIT_FAILURE;
		printf(" error: wrong dimensions of the intermodulation data array [%ld %ld %ld] \n",
				naxes[0],naxes[1],naxes[2]);
	}


	// actually read the table now
	int anynul=0;
	double* buf;
	buf = (double*)malloc(n_freq_s*n_freq_p*n_ener_p*sizeof(double));
	CHECK_MALLOC_VOID(buf);

	double nullval=0.0;
	long nelem = n_freq_s*n_freq_p*n_ener_p;
	long fpixel[naxis];
	for (int ii=0;ii<naxis;ii++){
		fpixel[ii]=1;
	}

	fits_read_pix(fptr, TDOUBLE, fpixel,nelem, &nullval,buf, &anynul, status);
	fits_report_error(stdout, *status);
	CHECK_STATUS_VOID(*status);

	headas_chat(7," start reading the electrical crosstalk matrix %s [%ld %ld %ld]",extname,n_freq_s,n_freq_p, n_ener_p);
	for (int ii=0; ii<n_freq_s ; ii++){  // freq_s loop
		for (int jj=0; jj<n_freq_p ; jj++){
			for (int kk=0; kk<n_ener_p ; kk++){
					int ind = 	( kk * n_freq_p + jj ) * n_freq_s + ii;
					assert(ind < nelem);
					tab->matrix[ii][jj][kk] = scaling*buf[ind]*1e-3;  // convert to keV !!!
			}
		}
	} // ------------------------  //  end (freq_s loop)
	headas_chat(7," ... done \n");

	free(buf);

	return;
}


/** load the electrical crosstalk tables */
void load_elec_table(AdvDet* det, int k ,int* status){

	// check if the table exists
	CHECK_NULL_VOID(det->crosstalk_elec_file,*status,"no file for the electrical crosstalk table specified");

	char* EXTNAME_FREQ_SIGNAL = "signal_frequency";
	char* EXTNAME_FREQ_PERTURBER = "perturber_frequency";
	char* EXTNAME_ENER_PERTURBER = "perturber_energy";
	char* COLNAME_FREQ_SIGNAL = "FREQ_S";
	char* COLNAME_FREQ_PERTURBER = "FREQ_P";
	char* COLNAME_ENER_PERTURBER = "EN_P";
	char* EXTNAME_FDM_CROSSTALK_GRAD;

	// Concatenating the grade
	if (asprintf(&EXTNAME_FDM_CROSSTALK_GRAD,"FDM_CROSSTALK_%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	fitsfile *fptr=NULL;

	do {
		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_elec_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the electrical crosstalk table %s for grade %i \n",fullfilename, k);
		/**headas_chat(5, "   ... reading the electrical crosstalk table %s \n",fullfilename);*/

		// read the extensions specifying the axes of the 3d matrix
		int n_freq_s, n_freq_p, n_ener_p;
		double* freq_s;
		double* freq_p;
		double* ener_p;
		get_imodtable_axis(&n_freq_s,&freq_s,EXTNAME_FREQ_SIGNAL,COLNAME_FREQ_SIGNAL,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_freq_p,&freq_p,EXTNAME_FREQ_PERTURBER,COLNAME_FREQ_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_p,&ener_p,EXTNAME_ENER_PERTURBER,COLNAME_ENER_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);
		// convert from eV to keV
		for (int ii=0; ii<n_ener_p; ii++){
			ener_p[ii] *=1e-3;
		}

		// initialize the elec crosstalk tables
		ElecTab* tmp_elec=&(det->crosstalk_elec[k]);
		initElecTab(&(tmp_elec), n_freq_s, n_freq_p, n_ener_p, freq_s, freq_p, ener_p, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing electrical crosstalk table in memory failed");
			break;
		}

		read_elec_matrix(fptr,n_freq_s,n_freq_p, n_ener_p, det->scaling, &(det->crosstalk_elec[k]),
				EXTNAME_FDM_CROSSTALK_GRAD,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading electrical crosstalk table %s  failed\n", fullfilename);
			break;
		}
		free(freq_s);
		free(freq_p);
		free(ener_p);
	} while(0); // END of Error handling loop

	free(EXTNAME_FDM_CROSSTALK_GRAD);
	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

void load_elec_timedep(AdvDet* det, int k, int* const status){
	CHECK_NULL_VOID(det->crosstalk_elec_timedep_file,*status,"no file for the elec time dependence specified");


	CrosstalkTimedep* tmp_lo=newCrossTalkTimedep(status);
	det->crosstalk_elec_timedep[2*k]=*tmp_lo;
	free(tmp_lo);
	CrosstalkTimedep* tmp_hi=newCrossTalkTimedep(status);
	det->crosstalk_elec_timedep[2*k+1]=*tmp_hi;
	free(tmp_hi);

	char* COLNAME_TIME_DELAY = "TIME_DELAY";
	char* COLNAME_WEIGHT = "WEIGHT";
	char* EXTNAME_LOW_GRAD;
	char* EXTNAME_HI_GRAD;
	// Concatenating the grade

	if (asprintf(&EXTNAME_LOW_GRAD,"PERTLO%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	if (asprintf(&EXTNAME_HI_GRAD,"PERTHI%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}


	fitsfile *fptr=NULL;

	do {
		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_elec_timedep_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the time dependence elecmal table %s for grade %i \n",fullfilename, k);

		// Read time dep info for first file
		int n_time_low, n_weight_low;
		double* time_low;
		double* weight_low;
		get_imodtable_axis(&n_time_low,&time_low,EXTNAME_LOW_GRAD,COLNAME_TIME_DELAY,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_weight_low,&weight_low,EXTNAME_LOW_GRAD,COLNAME_WEIGHT,fptr,status);
		CHECK_STATUS_BREAK(*status);

		det->crosstalk_elec_timedep[2*k].length=n_time_low;
		det->crosstalk_elec_timedep[2*k].time=time_low;
		det->crosstalk_elec_timedep[2*k].weight=weight_low;

		// now we need to set the weight at t=0 (differs for tessim simulated LUT)
		int ind_t0_low = binary_search(0.0, det->crosstalk_elec_timedep[2*k].time, det->crosstalk_elec_timedep[2*k].length);

		double fac_low = (0.0 - det->crosstalk_elec_timedep[2*k].time[ind_t0_low]) /
				(det->crosstalk_elec_timedep[2*k].time[ind_t0_low+1] - det->crosstalk_elec_timedep[2*k].time[ind_t0_low]);

		det->crosstalk_elec_timedep[2*k].weight_t0 =
				(1-fac_low)*det->crosstalk_elec_timedep[2*k].weight[ind_t0_low] +
				(fac_low)  *det->crosstalk_elec_timedep[2*k].weight[ind_t0_low+1];


		// Read time dependence info for second the second case (same grading, different frequency difference)
		int n_time_hi, n_weight_hi;
		double* time_hi;
		double* weight_hi;
		get_imodtable_axis(&n_time_hi,&time_hi,EXTNAME_HI_GRAD,COLNAME_TIME_DELAY,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_weight_hi,&weight_hi,EXTNAME_HI_GRAD,COLNAME_WEIGHT,fptr,status);
		CHECK_STATUS_BREAK(*status);

		det->crosstalk_elec_timedep[2*k+1].length=n_time_hi;
		det->crosstalk_elec_timedep[2*k+1].time=time_hi;
		det->crosstalk_elec_timedep[2*k+1].weight=weight_hi;

		// now we need to set the weight at t=0 (differs for tessim simulated LUT)
		int ind_t0_hi = binary_search(0.0, det->crosstalk_elec_timedep[2*k+1].time, det->crosstalk_elec_timedep[2*k+1].length);

		double fac_hi = (0.0 - det->crosstalk_elec_timedep[2*k+1].time[ind_t0_hi]) /
				(det->crosstalk_elec_timedep[2*k+1].time[ind_t0_hi+1] - det->crosstalk_elec_timedep[2*k+1].time[ind_t0_hi]);

		det->crosstalk_elec_timedep[2*k+1].weight_t0 =
				(1-fac_hi)*det->crosstalk_elec_timedep[2*k+1].weight[ind_t0_hi] +
				(fac_hi)  *det->crosstalk_elec_timedep[2*k+1].weight[ind_t0_hi+1];

	}while(0);

	free(EXTNAME_LOW_GRAD);
	free(EXTNAME_HI_GRAD);

	if (fptr!=NULL) {fits_close_file(fptr,status);}
}

/**static void print_elec_matrix(ElecTab* tab_cl,ElecTab* tab_ci){

	FILE * fp;
	fp = fopen ("elec_matrix_output_tab.txt", "w+");

	for (int kk=0; kk<tab_cl->n_ener_p ; kk++){
		for (int jj=0; jj<tab_cl->n_freq_p ; jj++){
				for (int ii=0; ii<tab_cl->n_freq_s ; ii++){
					fprintf(fp, "%i \t %i \t %i  \t %e \t %e \t %e \t %e \t %e \n",
							ii,jj,kk,tab_cl->freq_s[ii],tab_cl->freq_p[jj],tab_cl->ener_p[kk],
							tab_cl->matrix[ii][jj][kk],tab_ci->matrix[ii][jj][kk]);

			}
		}
	} // ------------------------  //  end (freq loop)

	fclose(fp);
}*/

// get the energy dependent weight (weight needs to be allocated before, values are added
// to the already given values in weight)
static void get_enerdep_weight(ElecTab* tab, double freq_s, double freq_p, double* weight, int n, int* status){

	assert(tab->n_ener_p==n);

	if ((freq_s > tab->freq_s[tab->n_freq_s-1] ) || (freq_s < tab->freq_s[0] )){
		printf(" *** error: signal pixel frequency %g outside tabulated range [%g,%g] \n",
				freq_s,tab->freq_s[0],tab->freq_s[tab->n_freq_s-1]);
		SIXT_ERROR("failed interpolating electrical crosstalk table");
		*status=EXIT_FAILURE;
		return;
	}
	if ((freq_p > tab->freq_p[tab->n_freq_p-1] ) || (freq_p < tab->freq_p[0] )){
		printf(" *** error: perturber pixel frequency %g outside tabulated range [%g,%g] \n",
				freq_p,tab->freq_p[0],tab->freq_p[tab->n_freq_p-1]);
		SIXT_ERROR("failed interpolating electrical crosstalk table");
		*status=EXIT_FAILURE;
		return;
	}

	// interpolate the arrays for the given frequencies
	int ind_freq_s = binary_search(freq_s,tab->freq_s,tab->n_freq_s);
	assert ( (ind_freq_s) >= 0 && (ind_freq_s < tab->n_freq_s-1) ); // should be valid after the above checks
	double d_freq_s =  (freq_s - tab->freq_s[ind_freq_s]) / (tab->freq_s[ind_freq_s+1] - tab->freq_s[ind_freq_s] );

	// interpolate the arrays for the given frequencies
	int ind_freq_p = binary_search(freq_p,tab->freq_p,tab->n_freq_p);
	assert ( (ind_freq_p) >= 0 && (ind_freq_p < tab->n_freq_p-1) );  // should be valid after the above checks
	double d_freq_p =  (freq_p - tab->freq_p[ind_freq_p]) / (tab->freq_p[ind_freq_p+1] - tab->freq_p[ind_freq_p] );

	// interpolate in 2d the weights
	for (int ii=0;ii<tab->n_ener_p;ii++){
		weight[ii] +=
				tab->matrix[ind_freq_s][ind_freq_p][ii]*(1-d_freq_s)*(1-d_freq_p) +
				tab->matrix[ind_freq_s+1][ind_freq_p][ii]*(d_freq_s)*(1-d_freq_p) +
				tab->matrix[ind_freq_s][ind_freq_p+1][ii]*(1-d_freq_s)*(d_freq_p) +
				tab->matrix[ind_freq_s+1][ind_freq_p+1][ii]*(d_freq_s)*(d_freq_p);
	}

}

// Loads electrical cross-talk for requested pixel
// Concretely, iterates over all the pixels of the channel
void load_electrical_cross_talk(AdvDet* det,int pixid,int* const status){
	CHECK_STATUS_VOID(*status);
	if (det->crosstalk_elec_file==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load electrical crosstalk with no corresponding file given at detector level");
		return;
	}
	if (det->crosstalk_elec_timedep_file==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load electrical crosstalk time dependence file with no corresponding file given at detector level");
		return;
	}

	// Load for each grade, every different table only once (on first iteration)
	if (det->crosstalk_elec == NULL){
		det->crosstalk_elec_timedep=(CrosstalkTimedep*)malloc(2*det->pix[0].ngrades*sizeof(CrosstalkTimedep));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_elec_timedep,*status);
		det->crosstalk_elec=(ElecTab*)malloc(det->pix[0].ngrades*sizeof(ElecTab));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_elec,*status);
		// load the tables for the electrical crosstalk
		for (int k=0; k<det->pix[pixid].ngrades; k++){
			load_elec_table(det, k, status);
			load_elec_timedep(det, k, status);
			CHECK_STATUS_VOID(*status);
		}
		assert(det->crosstalk_elec != NULL);
		assert(det->crosstalk_elec_timedep != NULL);
	}

	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;

	//Allocate memory
	concerned_pixel->electrical_cross_talk=newMatrixEnerdepCrossTalk(det->pix[pixid].ngrades, status);
	for (int k=0; k<det->pix[pixid].ngrades; k++){
		concerned_pixel->electrical_cross_talk[k].n_ener=det->crosstalk_elec[k].n_ener_p;
	}

	// Iterate over the channel
	for (int i=0;i<concerned_pixel->channel->num_pixels;i++){
		current_pixel = concerned_pixel->channel->pixels[i];
		// Cross-talk is not with yourself ;)
		if (current_pixel==concerned_pixel) continue;

		double freq_s = concerned_pixel->freq;
		double freq_p = current_pixel->freq;

		//Get the weights for each grading already
		for (int k=0; k<det->pix[pixid].ngrades; k++){

			// as the electrical crosstalk depends on the energy of the perturber, we need an array
			// of weights depending on the energy
			int n_weight=det->crosstalk_elec[k].n_ener_p;
			double weight[n_weight];

			for (int jj=0;jj<n_weight;jj++){
				weight[jj] = 0.0;
			}

			// Add Electrical Crosstalk
			get_enerdep_weight(&(det->crosstalk_elec[k]),freq_s,freq_p,weight,n_weight,status);

			// Add cross-talk pixel
			add_xt_enerdep_pixel(&(concerned_pixel->electrical_cross_talk[k]),current_pixel,weight,n_weight,status);
			CHECK_STATUS_VOID(*status);
		}
	}

}

static void initImodTab(ImodTab* tab, int n_ampl, int n_dt, int n_freq,
		double* ampl, double* dt, double* freq, int* status){

	// make a short pointer

	tab->n_ampl = n_ampl;
	tab->n_dt   = n_dt;
	tab->n_freq = n_freq;

	tab->ampl = (double*) malloc(n_ampl * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(tab->ampl,*status);
	for (int ii=0; ii<n_ampl; ii++){
		tab->ampl[ii] = ampl[ii];
	}


	tab->dt = (double*) malloc(n_dt * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(tab->dt,*status);

	for (int ii=0; ii<n_dt; ii++){
		tab->dt[ii] = dt[ii];
	}
	tab->dt_min=dt[0];
	tab->dt_max=dt[n_dt-1];

	tab->freq = (double*) malloc(n_freq * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(tab->freq,*status);

	for (int ii=0; ii<n_freq; ii++){
		tab->freq[ii] = freq[ii];
	}

	// allocate the 4d matrix (n_freq x n_dt x n_ampl x nampl)
	tab->matrix = (double****) malloc (n_freq*sizeof(double***));
	CHECK_MALLOC_VOID_STATUS(tab->matrix,*status);

	for (int ll=0; ll<n_freq; ll++){               // FREQ-LOOP
		tab->matrix[ll] = (double***) malloc (n_dt*sizeof(double**));
		CHECK_MALLOC_VOID_STATUS(tab->matrix[ll],*status);

		for (int ii=0; ii<n_dt; ii++){             // DT-LOOP
			tab->matrix[ll][ii] = (double**) malloc (n_ampl*sizeof(double*));
			CHECK_MALLOC_VOID_STATUS(tab->matrix[ll][ii],*status);

			for (int jj=0; jj<n_ampl; jj++){      // AMPL1-LOOP
				tab->matrix[ll][ii][jj] = (double*) malloc (n_ampl*sizeof(double));
				CHECK_MALLOC_VOID_STATUS(tab->matrix[ll][ii][jj],*status);

				for (int kk=0; kk<n_ampl; kk++){  // AMPL2-LOOP
					tab->matrix[ll][ii][jj][kk] = 0.0;
				}
			}
		}
	}
	return;
}

/** read one intermodulation matrix from a fits file and store it in the structure */
void read_intermod_matrix(fitsfile* fptr,int n_ampl, int n_dt, int n_freq,
		ImodTab* tab, char* extname, int* status){

	// ampl and dt arrays should be initialized before
	assert(tab->ampl!=NULL);
	assert(tab->dt!=NULL);
	assert(tab->freq!=NULL);

	// move to the extension containing the data
	int extver = 0;
//	if (fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status)){
	if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver ,status)){
		printf(" error: moving to extension %s failed \n",extname);
		return;
	}

	// get number of rows
//	int colnum;
//	if(fits_get_colnum(fptr, CASEINSEN,colname, &colnum, status)) return;
	int naxis;
	if (fits_get_img_dim(fptr, &naxis, status) && (naxis!=4)){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	long naxes[naxis];


	if (fits_get_img_size(fptr, naxis, naxes, status) ){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	// check dimensions
	if (naxes[0]!=n_freq || naxes[1]!=n_dt || naxes[2]!=n_ampl || naxes[3]!=n_ampl ){
		*status=EXIT_FAILURE;
		printf(" error: wrong dimensions of the intermodulation data array [%ld %ld %ld %ld] \n",
				naxes[0],naxes[1],naxes[2],naxes[3]);
	}

	// actually read the table now
	int anynul=0;
	double* buf;
	buf = (double*)malloc(n_freq*n_dt*n_ampl*n_ampl*sizeof(double));
	CHECK_MALLOC_VOID(buf);

	double nullval=0.0;
	long nelem = n_freq*n_dt*n_ampl*n_ampl;
	long fpixel[naxis];
	for (int ii=0;ii<naxis;ii++){
		fpixel[ii]=1;
	}

	fits_read_pix(fptr, TDOUBLE, fpixel,nelem, &nullval,buf, &anynul, status);
	fits_report_error(stdout, *status);
	CHECK_STATUS_VOID(*status);

	headas_chat(7," start reading the intermodulation crosstalk matrix [%i %i %i %i]",n_freq,n_dt,n_ampl,n_ampl);
	for (int ii=0; ii<n_freq ; ii++){  // freq loop
		for (int jj=0; jj<n_dt ; jj++){
			for (int kk=0; kk<n_ampl ; kk++){
				for (int ll=0; ll<n_ampl ; ll++){
					int ind = 	(( ll*n_ampl + kk ) * n_dt + jj ) * n_freq + ii;
					assert(ind < n_dt*n_ampl*n_ampl*n_freq);
					// amplitude is defined descending -> need to reverse it in output array
					tab->matrix[ii][jj][n_ampl-kk-1][n_ampl-ll-1] = buf[ind];
				}
			}
		}
	} // ------------------------  //  end (freq loop)
	headas_chat(7," ... done \n");
	free(buf);

	return;
}

/**static void print_imod_matrix(ImodTab* tab){

	FILE * fp;

	fp = fopen ("imod_matrix_output_tab.txt", "w+");

	for (int ll=0; ll<tab->n_ampl ; ll++){
		for (int kk=0; kk<tab->n_ampl ; kk++){
			for (int jj=0; jj<tab->n_dt ; jj++){
				for (int ii=0; ii<tab->n_freq ; ii++){  // freq loop
					fprintf(fp, "%i \t %i \t %i \t %i \t %e \t %e \t %e \t %e \t %e \n",
							ii,jj,kk,ll,tab->freq[ii],tab->dt[jj],tab->ampl[tab->n_ampl-kk-1],
							tab->ampl[tab->n_ampl-ll-1],
							tab->matrix[ii][jj][tab->n_ampl-kk-1][tab->n_ampl-ll-1]);

				}
			}
		}
	} // ------------------------  //  end (freq loop)

	fclose(fp);

}*/

void reverse_array_dbl(double* arr, int n){
	double t;
	for (int ii = 0; ii < n/2; ii++) {
		t           = arr[ii];
		arr[ii]     = arr[n-1-ii];
		arr[n-1-ii] = t;
	}
}

void storeEventtype(CrosstalkProxy* xtalk_proxy, int type, double df, int* status){
	int sgn;
	if (df<0.){ //Find sign of fs-fp
		sgn=-1;
	} else {
		sgn=1;
	}

	/** Allocation new space in the matrices*/
	xtalk_proxy->type=(int*)realloc(xtalk_proxy->type, (xtalk_proxy->n_active_crosstalk+1)*sizeof(int));
	xtalk_proxy->is_saved=(int*)realloc(xtalk_proxy->is_saved, (xtalk_proxy->n_active_crosstalk+1)*sizeof(int));


	/** Saving what kind of cross-talk mechanism, and the sign of the frequency difference */
	xtalk_proxy->type[xtalk_proxy->n_active_crosstalk]=type*sgn;
	xtalk_proxy->is_saved[xtalk_proxy->n_active_crosstalk]=0;

	if (abs(type)!=-IMODCTK && abs(type)!=-THERCTK && abs(type)!=-ELECCTK && abs(type)!=-PROPCTK1
			&& abs(type)!=-PROPCTK2 && abs(type)!=-DERCTK){
		printf("Found an incorrect cross-talk type %i \n", type);
		SIXT_ERROR("Wrong crosstalk type");
		*status=EXIT_FAILURE;
	}

}

/** load the intermodulation frequency table */
static void load_intermod_freq_table(AdvDet* det, int k, int* status){

	char* EXTNAME_AMPLITUDE = "amplitude";
	char* EXTNAME_TIMEOFFSET = "time_offset";
	char* EXTNAME_FREQOFFSET = "frequency_offset";
	char* COLNAME_AMPLITUDE = "AMPLITUDE";
	char* COLNAME_TIMEOFFSET = "DT_SEC";
	char* COLNAME_FREQOFFSET = "D_FREQ";
	char* EXTNAME_CROSSTALK_GRAD;

	if (asprintf(&EXTNAME_CROSSTALK_GRAD,"CROSSTALK_%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	fitsfile *fptr=NULL;

	do {

		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_intermod_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the intermodulation table %s\n",fullfilename);

		// read the extensions specifying the axes of the 3d matrix
		int n_ampl,n_dt,n_freq;
		double* ampl;
		double* dt;
		double* freq;
		get_imodtable_axis(&n_ampl,&ampl,EXTNAME_AMPLITUDE,COLNAME_AMPLITUDE,fptr,status);
		CHECK_STATUS_BREAK(*status);
		// amplitude is defined descending -> need to reverse it
		reverse_array_dbl(ampl,n_ampl);
		assert(ampl[0]<ampl[1]);

		get_imodtable_axis(&n_dt,&dt,EXTNAME_TIMEOFFSET,COLNAME_TIMEOFFSET,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_freq,&freq,EXTNAME_FREQOFFSET,COLNAME_FREQOFFSET,fptr,status);
		CHECK_STATUS_BREAK(*status);

		// initialize the intermodulation table
		initImodTab(&(det->crosstalk_imod_table[k]), n_ampl, n_dt, n_freq, ampl, dt, freq, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing intermodulation table in memory failed");
			break;
		}

		// set the minimal and maximal dt value (dt_min<0 and dt_max>0 required in the current implementation). NB: Hi and Low timedep files have the same dt (only wights differ).
		assert(det->crosstalk_imod_table[k].dt_min<0.0 && det->crosstalk_imod_table[k].dt_max>0.0);

		read_intermod_matrix(fptr,n_ampl,n_dt, n_freq, &(det->crosstalk_imod_table[k]),
				EXTNAME_CROSSTALK_GRAD,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading intermodulation table %s failed\n", fullfilename);
			break;
		}

//		print_imod_matrix(det->crosstalk_intermod_table);

	} while(0); // END of Error handling loop

	if (fptr!=NULL) {fits_close_file(fptr,status);}
	free(EXTNAME_CROSSTALK_GRAD);
	return;
}

/** Load intermodulation cross talk information into a single AdvPix*/
void load_intermod_cross_talk(AdvDet* det, int* status){

	CHECK_NULL_VOID(det->crosstalk_intermod_file,*status,"no file for the intermodulation table specified");
	/** Looping through all the grades*/
	if (det->crosstalk_imod_table==NULL){
		det->crosstalk_imod_table=(ImodTab*)malloc(det->pix[0].ngrades*sizeof(ImodTab));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_imod_table,*status);
		for(int k=0; k<det->pix->ngrades; k++){
			load_intermod_freq_table(det,k,status);
			CHECK_STATUS_VOID(*status);
		}
	}
}

static void initTDMTab(TDMTab** tab, int n_samples, int n_ener_p, int n_ener_v,
		double* samples, double* ener_p, double* ener_v, int* status){
	// make a short pointer
	TDMTab* t = (*tab);

	t->n_samples = n_samples;
	t->n_ener_p  = n_ener_p;
	t->n_ener_v = n_ener_v;

	t->samples = (double*) malloc(n_samples * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->samples,*status);
	for (int ii=0; ii<n_samples; ii++){
		t->samples[ii] = samples[ii];
	}

	t->ener_p = (double*) malloc(n_ener_p * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->ener_p,*status);
	for (int ii=0; ii<n_ener_p; ii++){
		t->ener_p[ii] = ener_p[ii];
	}

	t->ener_v = (double*) malloc(n_ener_v * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->ener_v,*status);
	for (int ii=0; ii<n_ener_v; ii++){
		t->ener_v[ii] = ener_v[ii];
	}

	// allocate the 3d matrix (n_samples x n_ener_p x n_ener_v)
	t->matrix = (double***) malloc (n_samples*sizeof(double**));
	CHECK_MALLOC_VOID_STATUS(t->matrix,*status);

	for (int ll=0; ll<n_samples; ll++){               // SAMPLE-LOOP
		t->matrix[ll] = (double**) malloc (n_ener_p*sizeof(double*));
		CHECK_MALLOC_VOID_STATUS(t->matrix[ll],*status);

		for (int ii=0; ii<n_ener_p; ii++){             // ENER_P-LOOP
			t->matrix[ll][ii] = (double*) malloc (n_ener_v*sizeof(double));
			CHECK_MALLOC_VOID_STATUS(t->matrix[ll][ii],*status);

			for (int jj=0; jj<n_ener_v; jj++){      // ENER_V-LOOP
				t->matrix[ll][ii][jj] = 0.0;
			}
		}
	}
	return;
}

/** read one electrical crosstalk matrix from a fits file and store it in the structure (units given in keV) */
void read_TDM_matrix(fitsfile* fptr,int n_samples, int n_ener_p, int n_ener_v, float scaling,
		TDMTab* tab, char* extname, int* status){

	// ampl and dt arrays should be initialized before
	assert(tab->samples!=NULL);
	assert(tab->ener_v!=NULL);
	assert(tab->ener_p!=NULL);

	// move to the extension containing the data
	int extver = 0;
	if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver ,status)){
		printf(" error: moving to extension %s failed \n",extname);
		return;
	}

	// get number of rows
	int naxis;
	if (fits_get_img_dim(fptr, &naxis, status) && (naxis!=4)){
		printf(" error: getting dimensions of proportional data array failed \n");
		return;
	}
	long naxes[naxis];


	if (fits_get_img_size(fptr, naxis, naxes, status) ){
		printf(" error: getting dimensions of proportional data array failed \n");
		return;
	}
	// check dimensions
	if (naxes[0]!=n_samples || naxes[1]!=n_ener_p || naxes[2]!=n_ener_v ){
		*status=EXIT_FAILURE;
		printf(" error: wrong dimensions of the proportional data array [%ld %ld %ld] \n",
				naxes[0],naxes[1],naxes[2]);
	}


	// actually read the table now
	int anynul=0;
	double* buf;
	buf = (double*)malloc(n_samples*n_ener_v*n_ener_p*sizeof(double));
	CHECK_MALLOC_VOID(buf);

	double nullval=0.0;
	long nelem = n_samples*n_ener_v*n_ener_p;
	long fpixel[naxis];
	for (int ii=0;ii<naxis;ii++){
		fpixel[ii]=1;
	}

	fits_read_pix(fptr, TDOUBLE, fpixel,nelem, &nullval,buf, &anynul, status);
	fits_report_error(stdout, *status);
	CHECK_STATUS_VOID(*status);

	headas_chat(7," start reading the proportional crosstalk matrix %s [%ld %ld %ld]",extname,n_samples,n_ener_p, n_ener_v);
	for (int ii=0; ii<n_samples ; ii++){  // freq_s loop
		for (int jj=0; jj<n_ener_p ; jj++){
			for (int kk=0; kk<n_ener_v ; kk++){
					int ind = 	( kk * n_ener_p + jj ) * n_samples + ii;
					assert(ind < nelem);
					tab->matrix[ii][jj][kk] = buf[ind]*scaling;  // already in keV !!!
			}
		}
	} // ------------------------  //  end (freq_s loop)
	headas_chat(7," ... done \n");

	free(buf);

	return;
}


void load_proportional_cross_talk(AdvDet* det,int pixid,int* const status){
	CHECK_STATUS_VOID(*status);
	if (det->TDM_prop_file==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load proportional crosstalk with no corresponding file given at detector level");
		return;
	}

	// Load table only once (on first iteration)
	if (det->crosstalk_TDM_prop == NULL){
		det->crosstalk_TDM_prop=(TDMTab*)malloc(det->pix[0].ngrades*sizeof(TDMTab));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_TDM_prop,*status);
		// load the tables for the proportional crosstalk
		for (int k=0; k<det->pix[pixid].ngrades; k++){
			load_prop_table(det, k, status);
			CHECK_STATUS_VOID(*status);
		}
		assert(det->crosstalk_TDM_prop != NULL);
	}

	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;

	//Allocate memory
	concerned_pixel->prop_cross_talk=newMatrixPropCrossTalk(status);
	int max_rows=concerned_pixel->channel->num_pixels; //rows in a given column, column start at 0 in structure (can differ)
	concerned_pixel->prop_cross_talk->cross_talk_pixels_1=(int*)realloc(concerned_pixel->prop_cross_talk->cross_talk_pixels_1,sizeof(int)); //only row N+1, column N
	concerned_pixel->prop_cross_talk->cross_talk_pixels_2=(int*)realloc(concerned_pixel->prop_cross_talk->cross_talk_pixels_2, (max_rows-1)*sizeof(int)); //All but itself in column N

	int circ_next_row=(concerned_pixel->row)%(det->max_rows)+1; //If pixel is at 'last' row, it will affect row 1

	//The last column may not have the same number of rows, if this is the case, since switches are sequential,
	//type 1 cross-talk cannot happen with row N+1 within the same column, we affect the next pixel to -1 since
	//it should not have a neighbor for type 1
	if (circ_next_row>max_rows){
		circ_next_row=-1;
	}
	//printf("Number of row %i, next row %i and max row %i\n", concerned_pixel->row, circ_next_row, max_rows);

	// Iterate over the channel to find all type 1/2 pixels
	for (int i=0;i<concerned_pixel->channel->num_pixels;i++){
		// Cross-talk is not with yourself ;)
		if (concerned_pixel->channel->pixels[i]->pindex==pixid) continue;

		current_pixel = &(det->pix[concerned_pixel->channel->pixels[i]->pindex]);

		int row = current_pixel->row;
		assert(current_pixel->channel->channel_id==concerned_pixel->channel->channel_id); //Check if we are in the same column indeed

		//Quick check if the pixel is concerned
		if (row==circ_next_row){
			//Type 1 for next row (N+1)
			concerned_pixel->prop_cross_talk->cross_talk_pixels_1[concerned_pixel->prop_cross_talk->type_1_pix]=current_pixel->pindex;
			concerned_pixel->prop_cross_talk->type_1_pix+=1;
		}

		//Type 2 for same column N
		concerned_pixel->prop_cross_talk->cross_talk_pixels_2[concerned_pixel->prop_cross_talk->type_2_pix]=current_pixel->pindex;
		concerned_pixel->prop_cross_talk->type_2_pix+=1;
	}

	assert(concerned_pixel->prop_cross_talk->type_1_pix<=1); //There should be 1 next row neighbor maximum
	assert(concerned_pixel->prop_cross_talk->type_2_pix==max_rows-1); //There should be max_rows-1 in column neighbors
	//printf("Number of type 1 pixels %i \n", concerned_pixel->prop_cross_talk->type_1_pix);
	//printf("Number of type 2 pixels %i \n", concerned_pixel->prop_cross_talk->type_2_pix);
}

/** load the proportional crosstalk tables */
void load_prop_table(AdvDet* det, int k ,int* status){

	// check if the table exists
	CHECK_NULL_VOID(det->TDM_prop_file,*status,"no file for the proportional crosstalk table specified");

	char* EXTNAME_SAMPLES = "OFFSETS";
	char* EXTNAME_ENER_VICTIM = "VENERGIES";
	char* EXTNAME_ENER_PERTURBER = "PENERGIES";
	char* COLNAME_SAMPLES = "OFFSET";
	char* COLNAME_ENER_VICTIM = "VENERGY";
	char* COLNAME_ENER_PERTURBER = "PENERGY";
	char* EXTNAME_PROP_CROSSTALK_GRAD;

	// Concatenating the grade
	if (asprintf(&EXTNAME_PROP_CROSSTALK_GRAD,"CROSSTALK_%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	fitsfile *fptr=NULL;

	do {
		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->TDM_prop_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the proportional crosstalk table %s for grade %i \n",fullfilename, k);

		// read the extensions specifying the axes of the 3d matrix
		int n_samples, n_ener_p, n_ener_v;
		double* samples;
		double* ener_v;
		double* ener_p;
		get_imodtable_axis(&n_samples,&samples,EXTNAME_SAMPLES,COLNAME_SAMPLES,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_v,&ener_v,EXTNAME_ENER_VICTIM,COLNAME_ENER_VICTIM,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_p,&ener_p,EXTNAME_ENER_PERTURBER,COLNAME_ENER_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);

		//Find the correct number of samples, which varies from grade to grade
		//For instance grade 1 has 1453 offset samples, grade 0 8181 etc... which are the samples + 1101 (the time before event)
		n_samples=(n_samples-det->pix[0].grades[0].gradelim_post)+det->pix[0].grades[k].gradelim_post;

		// initialize the elec crosstalk tables
		TDMTab* tmp_TDM=&(det->crosstalk_TDM_prop[k]);
		initTDMTab(&(tmp_TDM), n_samples, n_ener_p, n_ener_v, samples, ener_v, ener_p, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing proportional crosstalk table in memory failed");
			break;
		}

		read_TDM_matrix(fptr, n_samples, n_ener_p, n_ener_v,det->scaling,&(det->crosstalk_TDM_prop[k]),
				EXTNAME_PROP_CROSSTALK_GRAD,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading proportional crosstalk table %s  failed\n", fullfilename);
			break;
		}
		free(samples);
		free(ener_v);
		free(ener_p);
	} while(0); // END of Error handling loop

	free(EXTNAME_PROP_CROSSTALK_GRAD);
	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

void load_derivative_cross_talk(AdvDet* det,int pixid,int* const status){
	CHECK_STATUS_VOID(*status);
	if (det->TDM_der_file==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load derivative crosstalk with no corresponding file given at detector level");
		return;
	}

	// Load for each grade, every different table only once (on first iteration)
	if (det->crosstalk_TDM_der == NULL){
		det->crosstalk_TDM_der=(TDMTab*)malloc(det->pix[0].ngrades*sizeof(TDMTab));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_TDM_der,*status);
		// load the tables for the derivative crosstalk
		for (int k=0; k<det->pix[pixid].ngrades; k++){
			load_der_table(det, k, status);
			CHECK_STATUS_VOID(*status);
		}
		assert(det->crosstalk_TDM_der != NULL);
	}

	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;

	//Allocate memory
	concerned_pixel->der_cross_talk=newMatrixDerCrossTalk(status);
	int max_rows=concerned_pixel->channel->num_pixels; //row in a given column, column start at 0 in the structure (can differ)
	concerned_pixel->der_cross_talk->cross_talk_pixels=(int*)realloc(concerned_pixel->der_cross_talk->cross_talk_pixels,4*sizeof(int)); //for the moment only 4 considered TODO: Change once more accurate

	// Find the affected 'wiring' neighbors. This concerns row (N+1), (N+2), (N-1), (N-2)
	// in the same column. No issue as in type 1, since even columns with less rows affect the
	// next wiring neighbours.
	int circ_next_row=(concerned_pixel->row)%(max_rows)+1; //If pixel is at 'last' row, it will affect row 0
	int circ_next_next_row=(concerned_pixel->row+1)%(max_rows)+1;
	int circ_prev_row=(concerned_pixel->row+max_rows-2)%(max_rows)+1;
	int circ_prev_prev_row=(concerned_pixel->row+max_rows-3)%(max_rows)+1;
	//printf("Number of rows %i, row actual %i and other 4, %i, %i, %i, %i \n", max_rows, concerned_pixel->row, circ_prev_prev_row, circ_prev_row, circ_next_row, circ_next_next_row);


	// Iterate over the channel to find all type 1/2 pixels
	for (int i=0;i<concerned_pixel->channel->num_pixels;i++){
		// Cross-talk is not with yourself ;)
		if (concerned_pixel->channel->pixels[i]->pindex==pixid) continue;

		current_pixel = &(det->pix[concerned_pixel->channel->pixels[i]->pindex]);

		int row = current_pixel->row;
		assert(current_pixel->channel->channel_id==concerned_pixel->channel->channel_id); //Check if we are in the same column indeed

		//Quick check if the pixel is concerned
		if ((row==circ_next_row || row==circ_next_next_row || row==circ_prev_row || row==circ_prev_prev_row)){
			//Type 3 for next row (N+1), (N+2), (N-1), (N-2)
			concerned_pixel->der_cross_talk->cross_talk_pixels[concerned_pixel->der_cross_talk->num_cross_talk_pixels]=current_pixel->pindex;
			concerned_pixel->der_cross_talk->num_cross_talk_pixels+=1;
		}
	}
	//printf("Number of type 3 pixels %i \n", concerned_pixel->der_cross_talk->num_cross_talk_pixels);
	assert(concerned_pixel->der_cross_talk->num_cross_talk_pixels==4); //There should be 4 wiring neighbors
}

/** load the proportional crosstalk tables */
void load_der_table(AdvDet* det, int k ,int* status){

	// check if the table exists
	CHECK_NULL_VOID(det->TDM_der_file,*status,"no file for the derivative crosstalk table specified");

	char* EXTNAME_SAMPLES = "OFFSETS";
	char* EXTNAME_ENER_VICTIM = "VENERGIES";
	char* EXTNAME_ENER_PERTURBER = "PENERGIES";
	char* COLNAME_SAMPLES = "OFFSET";
	char* COLNAME_ENER_VICTIM = "VENERGY";
	char* COLNAME_ENER_PERTURBER = "PENERGY";
	char* EXTNAME_DER_CROSSTALK_GRAD;

	// Concatenating the grade
	if (asprintf(&EXTNAME_DER_CROSSTALK_GRAD,"CROSSTALK_%05ld",det->pix[0].grades[k].gradelim_post)== -1){
	  *status=EXIT_FAILURE;
	  printf(" *** error: allocating memory failed ");
	  return;
	}

	fitsfile *fptr=NULL;

	do {
		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->TDM_der_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the derivative crosstalk table %s for grade %i \n",fullfilename, k);
		/**headas_chat(5, "   ... reading the electrical crosstalk table %s \n",fullfilename);*/

		// read the extensions specifying the axes of the 3d matrix
		int n_samples, n_ener_p, n_ener_v;
		double* samples;
		double* ener_v;
		double* ener_p;
		get_imodtable_axis(&n_samples,&samples,EXTNAME_SAMPLES,COLNAME_SAMPLES,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_v,&ener_v,EXTNAME_ENER_VICTIM,COLNAME_ENER_VICTIM,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_p,&ener_p,EXTNAME_ENER_PERTURBER,COLNAME_ENER_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);

		//Find the correct number of samples, which varies from grade to grade
		//For instance grade 1 has 1453 offset samples, grade 0 8181 etc... which are the samples + 1101 (the time before event)
		n_samples=(n_samples-det->pix[0].grades[0].gradelim_post)+det->pix[0].grades[k].gradelim_post;

		// initialize the elec crosstalk tables
		TDMTab* tmp_TDM=&(det->crosstalk_TDM_der[k]);
		initTDMTab(&(tmp_TDM), n_samples, n_ener_p, n_ener_v, samples, ener_v, ener_p, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing derivative crosstalk table in memory failed");
			break;
		}

		read_TDM_matrix(fptr,n_samples, n_ener_p, n_ener_v,det->scaling,&(det->crosstalk_TDM_der[k]),
				EXTNAME_DER_CROSSTALK_GRAD,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading derivative crosstalk table %s  failed\n", fullfilename);
			break;
		}
		free(samples);
		free(ener_v);
		free(ener_p);
	} while(0); // END of Error handling loop

	free(EXTNAME_DER_CROSSTALK_GRAD);
	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

static int doCrosstalk(int id, AdvDet* det){
	/** crosstalk "ALL" now exculdes the IMOD crosstalk (does not show a large effect) **/
	if ( (id == det->crosstalk_id) || (det->crosstalk_id == CROSSTALK_ID_ALL )){
		return 1;
	} else {
		return 0;
	}
}

void init_crosstalk(AdvDet* det, int* const status){

	printf("\nCross-talk modes switched on: \n");

	// load the channel list
	if (det->readout_channels==NULL){
		if (det->tdm==0) {
			det->readout_channels = get_readout_channels(det,status);
		} else{
			det->readout_channels = get_readout_column(det,status);
		}
		CHECK_STATUS_VOID(*status);
	}

	if (det->pix[0].ngrades == 0){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load different grades for electrical cross-talk but got none");
	}

	// load thermal cross talk and associated time dependency
	if (det->xt_num_thermal>0 && doCrosstalk(CROSSTALK_ID_THERM,det)){
		printf(" - thermal crosstalk\n");
		det->crosstalk_ther_timedep=(CrosstalkTimedep**)malloc(det->xt_num_thermal*sizeof(CrosstalkTimedep*));
		CHECK_MALLOC_VOID_STATUS(det->crosstalk_ther_timedep,*status);
		for (int i=0; i<det->xt_num_thermal;i++){
			det->crosstalk_ther_timedep[i]=(CrosstalkTimedep*)malloc(2*det->pix[0].ngrades*sizeof(CrosstalkTimedep));
			CHECK_MALLOC_VOID_STATUS(det->crosstalk_ther_timedep[i],*status);
			for(int k=0; k<det->pix[0].ngrades; k++){
				load_thermal_timedep(det, i, k, status);
			}
		}

		for (int ii=0;ii<det->npix;ii++){
			load_thermal_cross_talk(det,ii,status);
			CHECK_STATUS_VOID(*status);
		}

	}

	if (det->tdm==0){
		// load electrical cross talk and associated time dependency
		if (det->crosstalk_elec_file!=NULL && doCrosstalk(CROSSTALK_ID_ELEC,det)){
			printf(" - electrical crosstalk (scaling=%.2e)\n", det->scaling);
			for (int ii=0;ii<det->npix;ii++){
				load_electrical_cross_talk(det,ii,status);
				CHECK_STATUS_VOID(*status);
			}
		}

		// load intermodulation cross talk and associated time dependency (NOT DONE FOR CROSSTALK==ALL)
		if (det->crosstalk_intermod_file!=NULL && doCrosstalk(CROSSTALK_ID_IMOD,det)){
			printf(" - intermodulation crosstalk\n");
			// loop through all pixels
			load_intermod_cross_talk(det,status);
			CHECK_STATUS_VOID(*status);
		}
	}

	if (det->tdm==1){
		// load proportional cross talk and associated time dependency (DONE ONLY WHEN TDM option is active)
		if (det->TDM_prop_file!=NULL && doCrosstalk(CROSSTALK_ID_TDM_PROP,det)){
			printf(" - proportional TDM crosstalk \n");
			// loop through all the pixels
			for (int ii=0;ii<det->npix;ii++){
				load_proportional_cross_talk(det, ii, status);
				CHECK_STATUS_VOID(*status);
			}

			/* addtional check: if tdm_prop1 or tdm_prop2 is given, only type1 or type2 will be simulated,
				and the other scaling factor set to 0 */
			char *buf;
			query_simput_parameter_string("doCrosstalk", &buf, status );
			if (strncmp(buf,"tdm_prop1",9)==0){
				det->prop_TDM_scaling_2=0.0;
				printf("    -> only simulating TYPE 1 proportional crosstalk (scaling=%.2e)\n", det->prop_TDM_scaling_1*det->scaling);
			} else if (strncmp(buf,"tdm_prop2",9)==0){
				det->prop_TDM_scaling_1=0.0;
				printf("    -> only simulating TYPE 2 proportional crosstalk (scaling=%.2e)\n", det->prop_TDM_scaling_2*det->scaling);
			}
			free(buf);
		}

		// load derivative cross talk and associated time dependency (DONE ONLY WHEN TDM option is active)
		if (det->TDM_der_file!=NULL && doCrosstalk(CROSSTALK_ID_TDM_DER,det)){
			printf(" - derivative TDM crosstalk \n");
			// loop through all the pixels
			for (int ii=0;ii<det->npix;ii++){
				load_derivative_cross_talk(det, ii, status);
				CHECK_STATUS_VOID(*status);
			}
		}
	}
}

static void add_pixel_to_readout(ReadoutChannels* read_chan, AdvPix* pix, int ic, int* status){

	// check if PixID already belongs to a Channel (duplicated naming not allowed)
	if (pix->channel != NULL){
		printf("*** error: Pixel with ID %i already belongs to a Channel (check channel_file)\n",pix->pindex);
		*status = EXIT_FAILURE;
		return;
	}

	if ( (ic <= 0) || ic > read_chan->num_channels ){
		printf("*** error: Channel with ID %i not a valid channel number\n",ic);
		*status = EXIT_FAILURE;
		return;
	}

	// add another pixel to the channel (ID starts at 1)
	// (note that we use here that realloc behaves like malloc for a NULL pointer)
	read_chan->channels[ic-1].pixels = (AdvPix**) realloc( read_chan->channels[ic-1].pixels,
					(read_chan->channels[ic-1].num_pixels+1) * sizeof(AdvPix*));
	CHECK_MALLOC_VOID_STATUS(read_chan->channels[ic-1].pixels,*status);
	read_chan->channels[ic-1].pixels[read_chan->channels[ic-1].num_pixels] = pix;
	read_chan->channels[ic-1].num_pixels++;

	return;
}


static int get_num_col(column_list* cl, int* status){

	int num_chans = -1;
	for(int ii=0; ii<cl->len; ii++){
		num_chans = MAX(num_chans,cl->col[ii]);
	}
	if (num_chans <= 0 ){
		CHECK_STATUS_RET(*status,0);
	}
	return num_chans;
}

static int get_num_chans(channel_list* cl, int* status){

	int num_chans = -1;
	for(int ii=0; ii<cl->len; ii++){
		num_chans = MAX(num_chans,cl->chan[ii]);
	}
	if (num_chans <= 0 ){
		CHECK_STATUS_RET(*status,0);
	}
	return num_chans;
}

static int get_num_rows(column_list* cl, int* status){

	int num_chans = -1;
	for(int ii=0; ii<cl->len; ii++){
		num_chans = MAX(num_chans,cl->row[ii]);
	}
	if (num_chans <= 0 ){
		CHECK_STATUS_RET(*status,0);
	}
	return num_chans;
}

ReadoutChannels* init_newReadout(int num_chan, int* status){

	ReadoutChannels* read_chan = (ReadoutChannels*) malloc(sizeof (ReadoutChannels) );
	CHECK_MALLOC_RET_NULL_STATUS(read_chan,*status);

	read_chan->num_channels = num_chan;

	read_chan->channels = (Channel*) malloc (num_chan*sizeof(Channel));
	CHECK_MALLOC_RET_NULL_STATUS(read_chan->channels,*status);

	read_chan->df_information_band = (double*) malloc (num_chan*sizeof(double));
	CHECK_MALLOC_RET_NULL_STATUS(read_chan->df_information_band,*status);

	for (int ii=0; ii<num_chan; ii++){
		read_chan->channels[ii].channel_id = ii+1; // channel IDs go from 1 to num_chan
		read_chan->channels[ii].num_pixels = 0;
		read_chan->channels[ii].pixels = NULL;
		read_chan->df_information_band[ii] = 0.0;
        read_chan->channels[ii].fdmsys = NULL;
	}

	return read_chan;
}

ReadoutChannels* get_readout_channels(AdvDet* det, int* status){

        // check if any channel file has been supplied
        if (det->channel_file == NULL){
                SIXT_ERROR("No channel file supplied");
                *status=EXIT_FAILURE;
                return NULL;
        }
	// get the complete file name
	char fullfilename[MAXFILENAME];
	strcpy(fullfilename,det->filepath);
	strcat(fullfilename,det->channel_file);

	// get the channel list
	channel_list* chans = load_channel_list(fullfilename, status);

	CHECK_STATUS_RET(*status,NULL);

	// check if the file agrees with the number of pixels
	if (chans->len != det->npix){
		printf("*** error: number of pixels from channel_file %s (%i) is not equal total number of pixels (%i)!\n",
				fullfilename,chans->len,det->npix);
		*status=EXIT_FAILURE;
		return NULL;
	}


	int num_chan = get_num_chans(chans,status);
	CHECK_STATUS_RET(*status,NULL);
	ReadoutChannels* read_chan = init_newReadout(num_chan, status);

	// sort pixels in the channel array
	AdvPix* pix=NULL;
	for (int ii=0; ii < chans->len; ii++ ){

		// check if PixID makes sense (PixID starts at 0)
		if ( (chans->pixid[ii] < 0) || (chans->pixid[ii] > det->npix-1) ){
			printf("*** error: Pixel-ID %i does not belong to the detector \n",chans->pixid[ii]);
			*status = EXIT_FAILURE;
			return NULL;
		}

		pix = &(det->pix[chans->pixid[ii]]);

		// set frequency of the pixel
		pix->freq = chans->freq[ii];
		pix->resfreq = pix->freq; // initialization, just to be safe
		if (pix->freq <= 0.0){
			printf("*** warning: assigning unrealistic frequency (%.3e) to pixel %i\n",pix->freq,pix->pindex);
		}

		// assign pixel to readout channel
		add_pixel_to_readout(read_chan, pix, chans->chan[ii], status);
		// make sure the pixel knows its channel (Channel ID starts at 1)
		pix->channel = &(read_chan->channels[chans->chan[ii]-1]);

//		printf("Readout Channel: Assigne PixID %i to Channel %i with Freq %.3e\n",
//				pix->pindex,pix->channel->channel_id,pix->freq);
	}

	// free the channel list
	free_channel_list(&chans);

	// set the information band for each channel now
    // set_df_information_band(read_chan,status);
	return read_chan;
}

/** Same as above but for TDM*/
ReadoutChannels* get_readout_column(AdvDet* det, int* status){

	// check if any channel file has been supplied
    if (det->channel_file == NULL){
    	SIXT_ERROR("No channel file supplied");
        *status=EXIT_FAILURE;
        return NULL;
    }

    // get the complete file name
	char fullfilename[MAXFILENAME];
	strcpy(fullfilename,det->filepath);
	strcat(fullfilename,det->channel_file);

	// get the channel list
	column_list* chans = load_column_list(fullfilename, status);

	CHECK_STATUS_RET(*status,NULL);

	// check if the file agrees with the number of pixels
	if (chans->len != det->npix){
		printf("*** error: number of pixels from channel_file %s (%i) is not equal total number of pixels (%i)!\n",
				fullfilename,chans->len,det->npix);
		*status=EXIT_FAILURE;
		return NULL;
	}

	int num_col = get_num_col(chans,status);
	det->max_rows=get_num_rows(chans,status);
	CHECK_STATUS_RET(*status,NULL);
	ReadoutChannels* read_chan = init_newReadout(num_col, status);

	// sort pixels in the channel array
	AdvPix* pix=NULL;
	for (int ii=0; ii < chans->len; ii++ ){

		// check if PixID makes sense (PixID starts at 0)
		if ( (chans->pixid[ii] < 0) || (chans->pixid[ii] > det->npix-1) ){
			printf("*** error: Pixel-ID %i does not belong to the detector \n",chans->pixid[ii]);
			*status = EXIT_FAILURE;
			return NULL;
		}

		pix = &(det->pix[chans->pixid[ii]]);

		// set row
		pix->row = chans->row[ii];
		if (pix->row<0){
			printf("*** warning: assigning unrealistic row (%i) to pixel %i\n",pix->row,pix->pindex);
		}

		// assign pixel to readout channel
		add_pixel_to_readout(read_chan, pix, chans->col[ii], status);
		// make sure the pixel knows its channel (Channel ID starts at 1)
		pix->channel = &(read_chan->channels[chans->col[ii]-1]);
	}
	// free the channel list
	free_column_list(&chans);
	return read_chan;
}

void free_channel_list(channel_list** chans){
	if (*(chans)!=NULL){
		free((*chans)->chan);
		free((*chans)->pixid);
		free((*chans)->freq);
		free(*chans);
	}
}

void free_column_list(column_list** chans){
	if (*(chans)!=NULL){
		free((*chans)->pixid);
		free((*chans)->row);
		free((*chans)->col);
		free(*chans);
	}
}

channel_list* load_channel_list(char* fname, int* status){

	channel_list* chans = (channel_list*)malloc(sizeof(channel_list));
	CHECK_MALLOC_RET_NULL(chans);

	// init parameters
	chans->len=0;
	chans->chan=NULL;
	chans->freq=NULL;
	chans->pixid=NULL;

	// open the file
	FILE *file=NULL;
	file=fopen(fname, "r");

	if (file == NULL){
		printf("*** error reading file %s \n",fname);
	 	 *status=EXIT_FAILURE;
	 	 return NULL;
	}

	int n=0;
	int* cha;
	int* pix;
	double* freq;
	while(!feof(file)){

	     pix = realloc(chans->pixid, (n+1)*sizeof(int));
	     cha = realloc(chans->chan,  (n+1)*sizeof(int));
	     freq = realloc(chans->freq,  (n+1)*sizeof(double));

	     if ((cha==NULL)||(pix==NULL)||(freq==NULL)){
	    	 	 free_channel_list(&chans);
	    	 	 SIXT_ERROR("error when reallocating arrays");
	    	 	 *status=EXIT_FAILURE;
	    	 	 return NULL;
	     } else {
	    	 chans->pixid = pix;
	    	 chans->chan  = cha;
	    	 chans->freq  = freq;
	     }
	     // check that always all three numbers are matched
	     if ((fscanf(file, "%i %i %lf\n",&(chans->pixid[n]),&(chans->chan[n]),&(chans->freq[n]))) != 3){
	    	 printf("*** formatting error in line %i of the channel file %s: check your input file\n",n+1,fname);
	    	 *status=EXIT_FAILURE;
	    	 return NULL;
	     }
	      // printf("reading channel list (line %i): %i %i %lf\n",n+1,chans->pixid[n],chans->chan[n],chans->freq[n]);
	     n++;
	}

	fclose(file);
    chans->len = n;
	return chans;
}

column_list* load_column_list(char* fname, int* status){

	column_list* chans = (column_list*)malloc(sizeof(column_list));
	CHECK_MALLOC_RET_NULL(chans);

	// init parameters
	chans->len=0;
	chans->col=NULL;
	chans->row=NULL;
	chans->pixid=NULL;

	// open the file
	FILE *file=NULL;
	file=fopen(fname, "r");

	if (file == NULL){
		printf("*** error reading file %s \n",fname);
	 	 *status=EXIT_FAILURE;
	 	 return NULL;
	}

	int n=0;
	int* col;
	int* pix;
	int* row;
	while(!feof(file)){

	     pix = realloc(chans->pixid, (n+1)*sizeof(int));
	     col = realloc(chans->col,  (n+1)*sizeof(int));
	     row = realloc(chans->row,  (n+1)*sizeof(double));

	     if ((col==NULL)||(pix==NULL)||(row==NULL)){
	    	 	 free_column_list(&chans);
	    	 	 SIXT_ERROR("error when reallocating arrays");
	    	 	 *status=EXIT_FAILURE;
	    	 	 return NULL;
	     } else {
	    	 chans->pixid = pix;
	    	 chans->col = col;
	    	 chans->row = row;
	     }
	     // check that always all three numbers are matched
	     if ((fscanf(file, "%i %i %i\n",&(chans->pixid[n]),&(chans->col[n]),&(chans->row[n]))) != 3){
	    	 printf("*** formatting error in line %i of the channel file %s: check your input file\n",n+1,fname);
	    	 *status=EXIT_FAILURE;
	    	 return NULL;
	     }
	      // printf("reading channel list (line %i): %i %i %i\n",n+1,chans->pixid[n],chans->col[n],chans->row[n]);
	     n++;
	}

	fclose(file);
    chans->len = n;
	return chans;
}

/** Compute influence of the crosstalk event on an impact using the timedependence table */
int computeCrosstalkInfluence(CrosstalkTimedep* buffer,PixImpact* impact,PixImpact* crosstalk, double energy, double* influence){

	double time_difference = crosstalk->time - impact->time;
	double energy_influence = 0.;

	if(buffer==NULL){
		SIXT_ERROR("No input time dependence matrix");
	}

	// if impact is close enough to have an influence
	if ((time_difference>buffer->time[0]) && (time_difference<buffer->time[buffer->length-1])){

		headas_chat(7,"TimeI:%.1e TimeC:%.2e CrosstalkE:%.2e ",impact->time,crosstalk->time,energy);
		// Binary search for to find interpolation interval
		int timedep_index=binary_search(time_difference,buffer->time,buffer->length);
		assert(timedep_index<buffer->length-1);

		// we need to weight to be one at t=0 (not the case for tessim LUT)
		double weight0 = buffer->weight[timedep_index]   / buffer->weight_t0;
		double weight1 = buffer->weight[timedep_index+1] / buffer->weight_t0;

		// influence previous impact
		energy_influence= energy*(weight0+
				(weight1-weight0) /(buffer->time[timedep_index+1]-buffer->time[timedep_index])
				*(time_difference-buffer->time[timedep_index]));
		impact->energy+=energy_influence;
		*influence+=energy_influence;
		headas_chat(7,", Influence:%.2e (fraction:%.2e)\n",*influence,*influence/energy);
		return 1;
	}
	return 0;
}


/** Cosntructor of CrosstalkProxy structure */
CrosstalkProxy* newCrosstalkProxy(int* const status){
	CrosstalkProxy* xtalk_proxy= (CrosstalkProxy*) malloc (sizeof (*xtalk_proxy));
	CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy,*status);

	xtalk_proxy->xtalk_impacts = (PixImpact**) malloc(INITXTALKNB*sizeof(PixImpact*));
	CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy->xtalk_impacts,*status);
	for (int ii=0;ii<INITXTALKNB;ii++){
		xtalk_proxy->xtalk_impacts[ii] = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy->xtalk_impacts[ii],*status);
	}
	xtalk_proxy->xtalk_proxy_size=INITXTALKNB;
	xtalk_proxy->n_active_crosstalk=0;
	xtalk_proxy->type=NULL;
	xtalk_proxy->is_saved=NULL;

	return xtalk_proxy;
}

/** Destructor of CrosstalkProxy structure */
void freeCrosstalkProxy(CrosstalkProxy** xtalk_proxy){
	if (*(xtalk_proxy)!=NULL){
		if ((*xtalk_proxy)->xtalk_impacts !=NULL){
			for (int ii=0;ii<((*xtalk_proxy)->xtalk_proxy_size);ii++){ //The other should be initialised but empty
				if((*xtalk_proxy)->xtalk_impacts[ii]!=NULL){
					free((*xtalk_proxy)->xtalk_impacts[ii]);
				}
			}
			free((*xtalk_proxy)->type);
			free((*xtalk_proxy)->is_saved);
			free((*xtalk_proxy)->xtalk_impacts);
			free(*xtalk_proxy);
		}
	}
}

/** Add crosstalk to proxy */
void addCrosstalk2Proxy(CrosstalkProxy* xtalk_proxy, float current_time, PixImpact* impact, int type, double df, int* const status){
	// If no more space, first check if the first events in the buffer are not too far from the current case (i.e. can we erase some)
	if (xtalk_proxy->n_active_crosstalk==xtalk_proxy->xtalk_proxy_size){
		int* toerase=NULL;
		int erased_crosstalks=0;
		float dt_current=0.;
		float dt_next=0.;
		for (int ii=0;ii<xtalk_proxy->xtalk_proxy_size;ii++){

			//This is the difference between the last impact time (!=simulation current time) and the cross-talk impact
			//If this is larger than DTMIN, the cross-talk event is too much in the past to affect, it can be erased
			dt_current=current_time-xtalk_proxy->xtalk_impacts[ii]->time;

			//This is the difference between the simulation current time and the cross-talk impact
			//If this is larger than DTMIN + DTMAX (i.e. cross-talk impact is too much in the past to influence current event) and if
			//-dt_current>DTMAX (i.e. the cross-talk impact is too much in the future to affect the last impact), we erase it.
			// /!\ THIS INTRODUCES A SMALL BIAS AS WE COULD HAVE A CROSS-TALK TRIGGER IN-BETWEEN THESE TWO TIMES
			// BUT 1) THIS HAPPENS EXTREMELY RARELY 2) EVEN IF THIS HAPPENS IT SHOULD NOT INFLUENCE THE CURRENT EVENT
			// AS TIME ELAPSED IS LARGE ENOUGH THIS APPROACH INCLUDES A POTENTIAL TRIGGER AT BORDER OF INFLUENCE
			dt_next=impact->time-xtalk_proxy->xtalk_impacts[ii]->time;

			if (dt_current>DTMIN){ //TODO Implement this as a variable not ad hoc
				toerase=(int*) realloc(toerase, (erased_crosstalks+1)*sizeof(int));
				toerase[erased_crosstalks]=ii;
				erased_crosstalks+=1;
			} else if (-dt_current>DTMAX&&dt_next>DTMIN+DTMAX){ //TODO Implement this as a variable not ad hoc
				//In this case, the last current event is clearly high-res, we erase all in between the last event and
				//the current time (minimal time of next impact).
				toerase=(int*) realloc(toerase, (erased_crosstalks+1)*sizeof(int));
				toerase[erased_crosstalks]=ii;
				erased_crosstalks+=1;
			}
		}

		if (erased_crosstalks>0){
			erasectk(xtalk_proxy, toerase, erased_crosstalks,status);
			free(toerase);
		}

	}
	// Check that there is room left for new crosstalk now
	if (xtalk_proxy->n_active_crosstalk==xtalk_proxy->xtalk_proxy_size){
		//Then allocate a new array of double size and copy old impacts
		PixImpact** piximp_list = (PixImpact**) malloc((xtalk_proxy->xtalk_proxy_size*2)*sizeof(PixImpact*));
		CHECK_MALLOC_VOID_STATUS(piximp_list,*status);
		for (int ii=0;ii<xtalk_proxy->xtalk_proxy_size;ii++){
			piximp_list[ii] = xtalk_proxy->xtalk_impacts[ii];
		}
		for (int ii=xtalk_proxy->xtalk_proxy_size;ii<2*xtalk_proxy->xtalk_proxy_size;ii++){
			piximp_list[ii] = (PixImpact*) malloc(sizeof(PixImpact));
			CHECK_MALLOC_VOID_STATUS(piximp_list[ii],*status);
		}
		free(xtalk_proxy->xtalk_impacts);
		xtalk_proxy->xtalk_impacts = piximp_list;
		xtalk_proxy->xtalk_proxy_size*=2;
	}

	// Copy crosstalk to proxy
	copyPixImpact(xtalk_proxy->xtalk_impacts[xtalk_proxy->n_active_crosstalk],impact);
	xtalk_proxy->xtalk_impacts[xtalk_proxy->n_active_crosstalk]->nb_pileup=0;
	xtalk_proxy->xtalk_impacts[xtalk_proxy->n_active_crosstalk]->weight_index=impact->weight_index;

	storeEventtype(xtalk_proxy, type, df, status);
	xtalk_proxy->n_active_crosstalk+=1;
}


/** Constructor of EventProxy structure */
EventProxy* newEventProxy(int* const status){
	EventProxy* proxy= (EventProxy*) malloc (sizeof (*proxy));
	CHECK_MALLOC_RET_NULL_STATUS(proxy,*status);

	proxy->impact = (PixImpact**) malloc(INITEVTPROXYNB*sizeof(PixImpact*));
	CHECK_MALLOC_RET_NULL_STATUS(proxy->impact,*status);
	for (int ii=0;ii<INITEVTPROXYNB;ii++){
		proxy->impact[ii] = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_RET_NULL_STATUS(proxy->impact[ii],*status);
	}

	proxy->ind_grading_previous=0;
	proxy->ind_grading_current=0;
	proxy->ind_grading_next=0;
	proxy->ind_event=0;
	proxy->nb_active=0;
	proxy->event_proxy_size=INITXTALKNB;

	return proxy;
}

/** Destructor of CrosstalkProxy structure */
void freeEventProxy(EventProxy** proxy){
	if (*(proxy)!=NULL){
		if ((*proxy)->impact !=NULL){
			for (int ii=0;ii<((*proxy)->event_proxy_size);ii++){
				free((*proxy)->impact[ii]);
			}
			free((*proxy)->impact);
			free(*proxy);
		}
	}
}

/** Computes energy deposited crosstalk energy on a given pixel */
void computeWeights(AdvDet* det, CrosstalkProxy* xtalk_proxy, PixImpact * impact, double* energies,
		int grade,int* const status){

	//Setting up proxy
	PixImpact* crosstalk=NULL;

	//Looping through all the current active ctk to compute weights wrt grading
	for (int ii=0;ii<xtalk_proxy->n_active_crosstalk;ii++){
		crosstalk = xtalk_proxy->xtalk_impacts[ii];

		//Getting the crosstalk type of the given impact
		//Intermod crosstalk (only if there is an impact indeed i.e. skip the first)
		//------------------
		if (abs(xtalk_proxy->type[ii])==-IMODCTK && impact->energy>0.) {
			//TODO This should be changed at some point, dE is assumed to be for dt=0, which is not the case...
			calc_imod_xt_influence(det,impact,crosstalk,&energies[ii],0,grade,status);
			CHECK_STATUS_VOID(*status);
			headas_chat(7,"  \n");

		//Thermal Crosstalk
		//------------------
		} else if (abs(xtalk_proxy->type[ii])==-THERCTK) {
			energies[ii] = crosstalk->energy*det->pix[crosstalk->pixID].thermal_cross_talk->cross_talk_weights[crosstalk->weight_index];

		//Electrical Crosstalk
		//----------------------
		} else if (abs(xtalk_proxy->type[ii])==-ELECCTK) {
			// make sure we have the same (length) energy array here
			assert(det->pix[impact->pixID].electrical_cross_talk[grade].n_ener == det->crosstalk_elec[grade].n_ener_p);
			double * ener_p = det->crosstalk_elec[grade].ener_p; // every electrical crosstalk energy vector has to be the same

			if ((crosstalk->energy < ener_p[0]) || (crosstalk->energy >= ener_p[det->pix[impact->pixID].electrical_cross_talk[grade].n_ener-1])) {
				headas_chat(7, " *** warning : impact event energy %g on pix %li is outside the tabulated values for electrical crosstalk [%g,%g]}n",
						crosstalk->energy, crosstalk->pixID, ener_p[0],ener_p[det->pix[crosstalk->pixID].electrical_cross_talk[grade].n_ener-1]);
				headas_chat(7, "     ---> skipping this event!\n");
			} else {
				// now determine the energy bin
				int ind = binary_search(crosstalk->energy ,det->crosstalk_elec[grade].ener_p,det->pix[impact->pixID].electrical_cross_talk[grade].n_ener);
				double fac = (crosstalk->energy - ener_p[ind]) / (ener_p[ind+1] - ener_p[ind]);
				energies[ii] =(1-fac)*det->pix[crosstalk->pixID].electrical_cross_talk[grade].cross_talk_weights[crosstalk->weight_index][ind]
					+(fac)*det->pix[crosstalk->pixID].electrical_cross_talk[grade].cross_talk_weights[crosstalk->weight_index][ind+1];
			}

		//Proportional Crosstalk
		//------------------
		} else if (abs(xtalk_proxy->type[ii])==-PROPCTK1 || abs(xtalk_proxy->type[ii])==-PROPCTK2){
			double dt_in_frames=0.;
			double energy_fictional_victim=0.;
			double crosstalk_effect=0.;
			calc_prop_xt_influence(det,energy_fictional_victim,crosstalk->energy,&crosstalk_effect, dt_in_frames, grade);
			if (abs(xtalk_proxy->type[ii])==-PROPCTK1){
				energies[ii]=crosstalk_effect*det->prop_TDM_scaling_1/1.e-2; //Scaled at 1% of amplitude;
			} else if (abs(xtalk_proxy->type[ii])==-PROPCTK2){
				energies[ii]=crosstalk_effect*det->prop_TDM_scaling_2/1.e-2; //Scaled at 1% of amplitude
			}

		//Derivative Crosstalk
		//------------------
		} else if (abs(xtalk_proxy->type[ii])==-DERCTK){
			double dt_in_frames=0.;
			double energy_fictional_victim=0.;
			double crosstalk_effect=0.;
			calc_der_xt_influence(det,energy_fictional_victim,crosstalk->energy,&crosstalk_effect, dt_in_frames, grade);
			energies[ii]=crosstalk_effect*det->der_TDM_scaling/1.e-2; //Scaled at 1% of amplitude;
		}
	}
}

/** Computes the effect of the deposited ctk energy computed on current event*/
void computeTimeDependency(AdvDet* det, CrosstalkProxy* xtalk_proxy,PixImpact * impact, double* energies, double* xtalk_energy,int* nb_influences,
		 TesEventFile* event_file, int save_crosstalk, int grade, double sample_length, int* const status){

	double energy_influence = 0.;
	int has_influenced= 0;
	int erased_crosstalks=0;
	int* toerase=NULL;
	PixImpact* crosstalk=NULL;
	PixImpact* to_save=(PixImpact*)malloc(sizeof(PixImpact));
	CHECK_MALLOC_VOID_STATUS(to_save,*status);

	for (int ii=0;ii<xtalk_proxy->n_active_crosstalk;ii++){
		//Getting the time dependency file
		crosstalk = xtalk_proxy->xtalk_impacts[ii];
		double dt = (crosstalk->time - impact->time);

		//If needed, we save first occurrence (to respect causality), modify the saving index (to avoid multiple saves)
		//and do it only if energy is non 0 (should never happen but failsafe)...
		if (save_crosstalk==1 && xtalk_proxy->is_saved[ii]==0 && energies[ii]!=0){
			copyPixImpact(to_save,crosstalk);
			to_save->pixID=impact->pixID; //We put back the ID of the pixel that receives the ctk (not the perturber)
			to_save->energy=energies[ii]; //We put the corresponding energy deposited (not full energy effect on pulse)
			addRMFImpact(event_file,to_save,-2,-2,xtalk_proxy->type[ii],0,energies[ii],status);
			xtalk_proxy->is_saved[ii]=1; //The event is now saveds
		}

		// If a previous crosstalk has no chance of influencing another events, we can clear it immediately and go on
		if (impact->time-crosstalk->time > DTMIN){
			toerase=(int*)realloc(toerase, (erased_crosstalks+1)*sizeof(int));
			toerase[erased_crosstalks]=ii;
			erased_crosstalks+=1;
			continue;
		}

		//Thermal or electrical cross-talk (only a LUT multiplication)
		if (abs(xtalk_proxy->type[ii])==-ELECCTK || abs(xtalk_proxy->type[ii])==-THERCTK){
			//Getting time dependency
			CrosstalkTimedep* buffer=NULL;
			buffer=getTimeDep(det, xtalk_proxy, ii, grade, status);

			// Compute influence if asked or if the crosstalk is going to be erased (otherwise multiple counting)
			// or if specifically asked to compute all the influences
			energy_influence=0.;
			has_influenced=computeCrosstalkInfluence(buffer,impact, crosstalk,energies[ii],&energy_influence);
			if (has_influenced){
				*nb_influences+=crosstalk->nb_pileup+1;
				*xtalk_energy+=energy_influence;
			}
		//Non-linear cross-talk
		} else if ((abs(xtalk_proxy->type[ii])==-IMODCTK)){
			if (dt <= det->crosstalk_imod_table[grade].dt_max && dt >= det->crosstalk_imod_table[grade].dt_min){
				energy_influence=0.;
				calc_imod_xt_influence(det,impact,crosstalk,&energy_influence,dt,grade,status);
				if (energy_influence!=0.){
					*nb_influences+=crosstalk->nb_pileup+1;
					*xtalk_energy+=energy_influence;
				}
			}
		//Proportional cross-talk
		} else if ((abs(xtalk_proxy->type[ii])==-PROPCTK1)|| (abs(xtalk_proxy->type[ii])==-PROPCTK2)){
			double dt_in_frames=dt/sample_length;
			energy_influence=0.;
			calc_prop_xt_influence(det,impact->energy,crosstalk->energy, &energy_influence, dt_in_frames, grade);
			if (energy_influence!=0.){
				if ((abs(xtalk_proxy->type[ii])==-PROPCTK1) && (det->prop_TDM_scaling_1>1e-9)){
					*nb_influences+=crosstalk->nb_pileup+1;
					*xtalk_energy+=energy_influence*det->prop_TDM_scaling_1/1.e-2; //Scaled at 1% of amplitude;
				} else if ((abs(xtalk_proxy->type[ii])==-PROPCTK2) && (det->prop_TDM_scaling_2>1e-9)){
					*nb_influences+=crosstalk->nb_pileup+1;
					*xtalk_energy+=energy_influence*det->prop_TDM_scaling_2/1.e-2; //Scaled at 1%
				}
			}
		//Derivative cross-talk
		} else if (abs(xtalk_proxy->type[ii])==-DERCTK){
			double dt_in_frames=dt/sample_length;
			energy_influence=0.;
			calc_der_xt_influence(det,impact->energy,crosstalk->energy,&energy_influence, dt_in_frames, grade);
			if (energy_influence!=0.){
				*nb_influences+=crosstalk->nb_pileup+1;
				*xtalk_energy+=energy_influence*det->der_TDM_scaling/1.e-2; //Scaled at 1% of amplitude;
			}
		}
	}

	if (erased_crosstalks>0){
		erasectk(xtalk_proxy, toerase, erased_crosstalks, status);
		free(toerase);
	}

	free(to_save);
}

/**Checks for pile-up or triggering ctk event*/
void checkTrigger(AdvDet* det, PixImpact* impact, CrosstalkProxy* xtalk_proxy, GradeProxy* grade_proxy, double* energies, TesEventFile* event_file,
		double next_time, double sample_length,int grade, int* is_trigger, int save_crosstalk,int* const status){

	int ii=1; //Starting index for pile-up
	double cumul_ener=0; //Cumulative energy

	while(ii<xtalk_proxy->n_active_crosstalk){
		PixImpact* previous_xtalk=(PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_VOID_STATUS(previous_xtalk,*status);
		copyPixImpact(previous_xtalk, xtalk_proxy->xtalk_impacts[ii-1]); //Copying it in case we need to save it (will be erased in proxy otherwise)
		int* toerase=NULL;
		int end=0; //Boolean in case of multiple pile-ups
		cumul_ener=energies[ii-1];

		//Check it the ctk event is above threshold.
		//Only consider those between current-dtmin and next+dtmax others will not influence
		//(Condition made in case of recursive call due to cross-talk trigger)
		if((previous_xtalk->time<next_time-DTMIN || next_time==-1.0) && energies[ii-1]>det->threshold_event_lo_keV){
			toerase=(int*)realloc(toerase, sizeof(int));
			toerase[0]=ii-1;

			//Does it happen to pile up? (if more than one event in proxy)
			//If so, add all the energies corresponding to the event erase it from proxy then process it!
			if(xtalk_proxy->n_active_crosstalk>1 &&
					fabs(xtalk_proxy->xtalk_impacts[ii]->time-previous_xtalk->time) <= sample_length){
				previous_xtalk->nb_pileup++;// Update number of pileups
				cumul_ener+=energies[ii]; //Cumul energy
				toerase=(int*)realloc(toerase, (previous_xtalk->nb_pileup+1)*sizeof(int)); //In case we have to erase the event
				toerase[previous_xtalk->nb_pileup]=ii;

				while(end==0){ //Counting how many are actually piling up
					if(ii+previous_xtalk->nb_pileup<xtalk_proxy->n_active_crosstalk){ //If there actually is an event afterwards
						PixImpact* next_xtalk = xtalk_proxy->xtalk_impacts[ii+previous_xtalk->nb_pileup];
						if(next_xtalk->time-previous_xtalk->time <= sample_length){
							previous_xtalk->nb_pileup++;// Update number of pileups
							toerase=(int*)realloc(toerase, (previous_xtalk->nb_pileup+1)*sizeof(int));
							toerase[previous_xtalk->nb_pileup]=ii+previous_xtalk->nb_pileup-1;
							cumul_ener+=energies[ii+previous_xtalk->nb_pileup-1];
						} else{
							end=1;
						}
					} else{
						end=1;
					}
				}
			}

			//Failsafe for causality, this condition is always true normally. Process new event
			if(xtalk_proxy->xtalk_impacts[ii-1+previous_xtalk->nb_pileup]->time>impact->time){
				//Copying information
				previous_xtalk->time = xtalk_proxy->xtalk_impacts[ii-1+previous_xtalk->nb_pileup]->time;// to conserve causality
				previous_xtalk->energy=cumul_ener;
				previous_xtalk->pixID=impact->pixID; //As otherwise we have the pertuber index

				//Reset grade proxy to 0 as we will reprocess the event
				grade_proxy->nb_crosstalk_influence=0.;
				grade_proxy->crosstalk_energy=0.;

				//Erasing the event from proxy
				erasectk(xtalk_proxy, toerase, previous_xtalk->nb_pileup+1,status); // Erase the previous event(s), which now form(s) a pulse
				free(toerase);
				end=0;

				// Process now triggered event and exit the function
				*is_trigger=1;
				processGradedEvent(grade_proxy,sample_length,previous_xtalk,det,event_file,1,save_crosstalk,grade,status);
				free(previous_xtalk);
				return;
			}
			ii+=previous_xtalk->nb_pileup+1;
			free(toerase);
			end=0;

		//If the impact is not above threshold by itself check if a possible pile-up may trigger (if more than one event)
		} else if(xtalk_proxy->n_active_crosstalk>1 && (previous_xtalk->time<next_time-DTMIN  || next_time==-1.0)
				&& fabs(xtalk_proxy->xtalk_impacts[ii]->time-previous_xtalk->time) <= sample_length) { // TODO: code pileup length at pixel level (or use timedep table)
			previous_xtalk->nb_pileup++;// Update number of pileups
			cumul_ener+=energies[ii+previous_xtalk->nb_pileup-1]; //Cumul energy
			toerase=(int*)realloc(toerase, (previous_xtalk->nb_pileup+1)*sizeof(int)); //In case we have to erase the event
			toerase[previous_xtalk->nb_pileup-1]=ii-1; //If it triggers this event will need to be erased
			toerase[previous_xtalk->nb_pileup]=ii;

			while(end==0){ //Counting how many are actually piling up
				if(ii+previous_xtalk->nb_pileup<xtalk_proxy->n_active_crosstalk){ //If there actually is an event afterwards
					PixImpact* next_xtalk = xtalk_proxy->xtalk_impacts[ii+previous_xtalk->nb_pileup];
					if(next_xtalk->time-previous_xtalk->time <= sample_length){
						previous_xtalk->nb_pileup++;// Update number of pileups
						toerase=(int*)realloc(toerase, (previous_xtalk->nb_pileup+1)*sizeof(int));
						toerase[previous_xtalk->nb_pileup]=ii+previous_xtalk->nb_pileup-1;
						cumul_ener+=energies[ii+previous_xtalk->nb_pileup-1];
					} else{
						end=1;
					}
				} else{
					end=1;
				}
			}

			//Now check if this piled up event does trigger
			if (cumul_ener >= det->threshold_event_lo_keV) {

				//Failsafe for causality, this condition is always true normally. Process new event
				if(xtalk_proxy->xtalk_impacts[ii-1+previous_xtalk->nb_pileup]->time>impact->time){
					//Copying information
					previous_xtalk->time = xtalk_proxy->xtalk_impacts[ii-1+previous_xtalk->nb_pileup]->time;// to conserve causality
					previous_xtalk->energy=cumul_ener;
					previous_xtalk->pixID=impact->pixID; //As otherwise we have the pertuber index /!\ HERE BECAUSE WE ERASE!

					//Reset grade proxy to 0 as we will reprocess the event
					grade_proxy->nb_crosstalk_influence=0.;
					grade_proxy->crosstalk_energy=0.;

					//Erasing the event
					erasectk(xtalk_proxy, toerase, previous_xtalk->nb_pileup+1,status); // Erase the previous event(s), which now form(s) a pulse
					free(toerase);
					end=0;
					// Process now triggered event and exit the function
					*is_trigger=1;
					processGradedEvent(grade_proxy,sample_length,previous_xtalk,det,event_file,1,save_crosstalk,grade,status);
					free(previous_xtalk);
					return;
				}
			}
			free(toerase);
			end=0;
			ii+=previous_xtalk->nb_pileup+1;

		//If neither pile-up nor trigger occurs, just go on
		} else{
			ii++;
		}
		free(previous_xtalk);
	}
}

/** Compute crosstalk influence */
void computeAllCrosstalkInfluence(AdvDet* det,PixImpact * impact,CrosstalkProxy* xtalk_proxy, GradeProxy* grade_proxy, TesEventFile* event_file,double* xtalk_energy,int* nb_influences,
		double next_time, double sample_length, int* is_trigger, int save_crosstalk, int grade, int* const status){

	// If there is no crosstalk in the proxy get out without any influence
	if (xtalk_proxy==NULL || xtalk_proxy->n_active_crosstalk==0){
		*xtalk_energy=0.;
		*nb_influences=0;
		return;
	}

	// Initialise the deposited dE on the victim pixel
	//TODO: This could also be done on the addCrosstalk2Proxy level, but more flexible here, although slower
	double energies[xtalk_proxy->n_active_crosstalk];
	for(int jj=0;jj<xtalk_proxy->n_active_crosstalk;jj++){
		energies[jj]=0;
	}

	//Computing deposited energy of active ctk events on given victim event
	computeWeights(det, xtalk_proxy, impact, energies, grade, status);

	//Now we check whether two of the events can trigger
	//If pile-up occurs in the pixel between the two current event (should be rare!)
	checkTrigger(det,impact,xtalk_proxy,grade_proxy,energies,event_file,next_time,sample_length,grade,is_trigger,save_crosstalk,status);

	//If no outside trigger occurs, then we compute the values of the affected energy via the time dependency and save
	//(Not before as we did not know if there was a trigger)
	if (grade_proxy->is_first==1 && *is_trigger==0){
		computeTimeDependency(det,xtalk_proxy,impact, energies, xtalk_energy,nb_influences,
				 event_file, save_crosstalk, grade, sample_length, status);
	}
}
