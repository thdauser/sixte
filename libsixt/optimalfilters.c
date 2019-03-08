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

#include "optimalfilters.h"

/** OptimalFilterCollection constructor. Returns a pointer to an empty OptimalFilterCollection data
    structure. */
OptimalFilterCollection* newOptimalFilterCollection(int* const status){
  OptimalFilterCollection* opt_filter_collection= malloc(sizeof*opt_filter_collection);
  if (NULL==opt_filter_collection) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for OptimalFilterCollection failed");
    return(opt_filter_collection);
  }

  opt_filter_collection->nfilters = 0;
  opt_filter_collection->ph_a = 0;
  opt_filter_collection->ph_b = 0;
  opt_filter_collection->optimal_filters = NULL;

  return(opt_filter_collection);
}

/** Allocates memory for an OptimalFilterCollection structure. */
void allocateOptimalFilterCollection(OptimalFilterCollection* opt_filter_collection,int nfilters,int* const status){
	opt_filter_collection->nfilters = nfilters;
	opt_filter_collection->optimal_filters = malloc(nfilters*sizeof(*(opt_filter_collection->optimal_filters)));
	if (NULL==opt_filter_collection->optimal_filters) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for OptimalFilterCollection->optimal_filters failed");
	}
//	for (int i = 0 ; i < nfilters ; i++){
//		opt_filter_collection->optimal_filters[i] = newOptimalFilterCollection(status);
//	}
}

/** OptimalFilterCollection Destructor. */
void freeOptimalFilterCollection(OptimalFilterCollection** const opt_filter_collection){
  if (NULL!=*opt_filter_collection) {
	  for (int i = 0 ; i < (*opt_filter_collection)->nfilters ; i++){
		  freeOptimalFilter(&(*opt_filter_collection)->optimal_filters[i]);
	  }
	  (*opt_filter_collection)->nfilters=0;
	  free(*opt_filter_collection);
  }
}

/** OptimalFilter constructor. Returns a pointer to an empty OptimalFilter data structure. */
OptimalFilter* newOptimalFilter(int* const status){
	OptimalFilter* opt_filter = malloc(sizeof*opt_filter);
	if (NULL==opt_filter) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for OptimalFilter failed");
		return(opt_filter);
	}

	opt_filter->filter_duration=0;
	opt_filter->filter=NULL;

	return(opt_filter);
}

/** Allocates memory for an OptimalFilter structure. */
void allocateOptimalFilter(OptimalFilter* opt_filter,int filter_duration,int* const status){
	opt_filter->filter_duration=filter_duration;
	opt_filter->filter=malloc(filter_duration*sizeof(double));
	if (NULL==opt_filter->filter){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for OptimalFilter->filter failed");
	}
}

/** OptimalFilter destructor. */
void freeOptimalFilter(OptimalFilter* opt_filter){
	if (NULL!=opt_filter) {
		free(opt_filter->filter);
		opt_filter->filter_duration=0;
	}
}


/** Create and retireve an OptimalFilterCollection from a file. */
OptimalFilterCollection*  getOptimalFilterCollection(const char* const filename, int nfilters, int* const status){

	//Create OptimalFilterCollection structure
	OptimalFilterCollection* opt_filter_collection = newOptimalFilterCollection(status);
	CHECK_STATUS_RET(*status, opt_filter_collection);

	//Allocate OptimalFilterCollection structure
	allocateOptimalFilterCollection(opt_filter_collection,nfilters,status);
	CHECK_STATUS_RET(*status, opt_filter_collection);



	//Open FITS file in READONLY mode
	fitsfile* fptr = NULL;
	if (fits_open_file(&fptr, filename, READONLY, status)) return(opt_filter_collection);

	//Move to the first hdu
	int hdu_type;
	if (fits_movabs_hdu(fptr,2, &hdu_type, status)) return(opt_filter_collection);

	//Iterate over the filters and populate the structure
	char column_name[12];
	int column_number=0;
	int anynul=0;
	for (int i = 0 ; i < nfilters ; i++){
	    int filter_duration=(int)pow(2,i+2);
		sprintf(column_name, "OPTFILT%04d",filter_duration);

	    //Get column number
	    fits_get_colnum(fptr, CASEINSEN,column_name, &column_number, status);
		CHECK_STATUS_RET(*status, opt_filter_collection);

		//Allocate corresponding OptimalFilter
		allocateOptimalFilter(&(opt_filter_collection->optimal_filters[i]),filter_duration,status);
		CHECK_STATUS_RET(*status, opt_filter_collection);

		//Actually read column and save it into the structure
		fits_read_col(fptr, TDOUBLE, column_number, 1,1,filter_duration,NULL, opt_filter_collection->optimal_filters[i].filter, &anynul, status);
		CHECK_STATUS_RET(*status, opt_filter_collection);
	}
	//Get the pulse height to keV conversion parameters
	sprintf(column_name, "PH_A");
	fits_get_colnum(fptr, CASEINSEN,column_name, &column_number, status);
	CHECK_STATUS_RET(*status, opt_filter_collection);
	fits_read_col(fptr, TDOUBLE, column_number, 1,1,1,NULL, &(opt_filter_collection->ph_a), &anynul, status);
	CHECK_STATUS_RET(*status, opt_filter_collection);
	sprintf(column_name, "PH_B");
	fits_get_colnum(fptr, CASEINSEN,column_name, &column_number, status);
	CHECK_STATUS_RET(*status, opt_filter_collection);
	fits_read_col(fptr, TDOUBLE, column_number, 1,1,1,NULL, &(opt_filter_collection->ph_b), &anynul, status);
	CHECK_STATUS_RET(*status, opt_filter_collection);

	fits_close_file(fptr, status);
	CHECK_STATUS_RET(*status, opt_filter_collection);

	return(opt_filter_collection);
}



/** Derivate a stream */
void derivate_stream(double * data_stream,double * derivated_stream,unsigned long stream_length){
	for (unsigned long i=0 ; i < stream_length-1 ; i++) {
		derivated_stream[i] = data_stream[i+1]-data_stream[i];
	}
}

/** Remove a pulse from a data stream */
void subtractPulse(double * data_stream,unsigned long pulse_time,double * pulse_template,double factor,int pulse_length,unsigned long stream_length){
	for (unsigned long i=pulse_time; (i < stream_length && i < (pulse_time+pulse_length)); i++){
		data_stream[i] -= factor*pulse_template[i-pulse_time];
	}
}

/** Filter the pulse and return the energy */
double filterPulse(double * data_stream,int pulse_time,double * filter,int filter_length){
	double energy=0.;
	for (int i=0; i < filter_length ; i++){
		energy = energy + data_stream[pulse_time+i]*filter[i];
	}
	return(energy);
}

/** Trigger on the pulses and updates the event_list accordingly */
int triggerEvents(TesRecord* record,TesEventList* event_list,double * derivated_pulse,int derivated_pulse_length,
		double threshold,double pulse_template_height,double saturation_value,int derivate_exclusion,
		int normal_exclusion,int* const status){
	//Derivate the data stream
	double *derivated_stream = malloc((record->trigger_size - 1)*sizeof(*derivated_stream));
	if (NULL==derivated_stream){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for derivated_stream failed");
		CHECK_STATUS_RET(*status,1);
	}
	derivate_stream(record->adc_double,derivated_stream,record->trigger_size);

	//Iterate over the data points and trigger
	unsigned long ii = 1;
	int pulse_detected=0;
	double pulse_height=0;
	unsigned long pulse_time = 0;
	int grade1 = 0;
	int distance_to_wrong_reconstruction = normal_exclusion;
	unsigned long ii_end = record->trigger_size - 1;
	char current_value_is_saturated=0;
	char next_value_is_saturated=0;
	//If we start already above half of the threshold, the start of the stream is probably not event clean
	if (derivated_stream[0] > 0.5*threshold){
		distance_to_wrong_reconstruction = 0;
	}
	while(ii < ii_end){
		//Handle incomplete record
		if (record->adc_double[ii+1]==0){
			break;
		}

		//Detect when we are in the tail of a non-suppressed pulse
		if (derivated_stream[ii]<-0.5*threshold){
			distance_to_wrong_reconstruction=0;
		} else{
			distance_to_wrong_reconstruction++;
		}

		//Detect saturated pulses
		current_value_is_saturated = (record->adc_double[ii] >= saturation_value || record->adc_double[ii+1] >= saturation_value);
		if (current_value_is_saturated && !pulse_detected){
			if (event_list->index>0 && event_list->grades1[event_list->index-1] != -3){
				event_list->grades1[event_list->index-1] = -2;
			}
		} else {

			//Triggering part
			if(!pulse_detected){
				if(derivated_stream[ii]>threshold){
					pulse_height = derivated_stream[ii];
					pulse_time = ii-1;
					pulse_detected=1;
				}
				//We triggered, look for pulse height
			} else {
				next_value_is_saturated = (ii<ii_end-1 && (record->adc_double[ii+1] >= saturation_value || record->adc_double[ii+2] >= saturation_value));
				if(derivated_stream[ii] > pulse_height && !next_value_is_saturated){
					pulse_height = derivated_stream[ii];
					//We have reached the top of the pulse -> check for pileup, remove pulse from stream and add event to list
				} else {
					//It the pulse peak is saturated, estimate pulse height from last correct value
					if(current_value_is_saturated){
						pulse_height = (derivated_stream[pulse_time+1]-derivated_stream[pulse_time])/derivated_pulse[0]*pulse_template_height;
					} else if (ii<ii_end-1 && next_value_is_saturated){
						pulse_height = (pulse_height-derivated_stream[pulse_time])/derivated_pulse[ii-pulse_time-1]*pulse_template_height;
					}

					//Test if we have a pileup by checking whether the pulse subtraction would be too big
					if((derivated_stream[pulse_time+1] - pulse_height/pulse_template_height*derivated_pulse[0]) < -0.75*threshold){
						//This is a pileup, re-evaluate pulse height
						pulse_height = derivated_stream[pulse_time+1]/derivated_pulse[0]*pulse_template_height;
					}

					//Detect wrong reconstructions
					if (distance_to_wrong_reconstruction < derivate_exclusion || (pulse_time>0 && record->adc_double[pulse_time-1] >= saturation_value) || record->adc_double[pulse_time] >= saturation_value || record->adc_double[pulse_time+1] >= saturation_value){
						grade1 = -3;
					} else if (distance_to_wrong_reconstruction < normal_exclusion){
						grade1 = -2;
					} else {
						grade1 = -1;
					}

					//Remove pulse from stream for future detections
					subtractPulse(derivated_stream,pulse_time+1,derivated_pulse,pulse_height/pulse_template_height,derivated_pulse_length,record->trigger_size - 1);

					//If after subtraction, there is some residual signal, flag the event as misreconstructed
					if(derivated_stream[pulse_time] > 0.75*threshold || derivated_stream[pulse_time]<-0.75*threshold){
						grade1 = -3;
					}

					//Add Event to list
					addEventToList(event_list,pulse_time+1,pulse_height,grade1,status);
					CHECK_STATUS_RET(*status,ii+1);

					//Reboot search to pulse time
					ii = pulse_time+1;
					pulse_detected = 0;
				}
			}
		}
		ii++;
	}
	free(derivated_stream);
	return(ii+1);
}

/** Computes the energy of the detected pulses and save the result in the event list */
void computeEnergy(TesRecord* record,TesEventList* event_list,OptimalFilterCollection* opt_filter_collection,double * pulse_template,
		int pulse_length,double calfac,const char identify,int* const status){
	//Add dummy event to correctly compute the grades
	addEventToList(event_list,record->trigger_size-1,0.,0,status);
	CHECK_STATUS_VOID(*status);

	//Nominal optimal_filter index
	int nominal_filter_index = floor(log(pulse_length)/log(2)) - 2;

	//Iterate over pulses
	int pulse_time;
	int filter_index;
	int current_pulse_duration;
	double energy;
	long ph_id;
	double impact_time;
	double pulse_real_time;
	double pulse_distance;
	char pulse_identified=0;
	if(identify){
		popPhID(record->phid_list,&ph_id,&impact_time,status);
	}
	for (int i=0; i < event_list->index -1 ; i++){
		pulse_time = event_list->event_indexes[i];
		pulse_distance = event_list->event_indexes[i+1] - pulse_time;
		if(pulse_distance>=pulse_length && event_list->grades1[i] > -2){
			current_pulse_duration = opt_filter_collection->optimal_filters[nominal_filter_index].filter_duration;
			//Compute energy
			event_list->energies[i] = filterPulse(record->adc_double,pulse_time,opt_filter_collection->optimal_filters[nominal_filter_index].filter,
					current_pulse_duration)/calfac;
			event_list->grades1[i] = current_pulse_duration;
		} else if(pulse_distance>=4 && event_list->grades1[i] > -2) {
			//Get available optimal filter index
			filter_index = floor(log(pulse_distance)/log(2)) - 2;
			current_pulse_duration = opt_filter_collection->optimal_filters[filter_index].filter_duration;
			//Compute energy
			energy = filterPulse(record->adc_double,pulse_time,opt_filter_collection->optimal_filters[filter_index].filter,
					current_pulse_duration);
			//Remove pulse from data
			subtractPulse(record->adc_double,pulse_time,pulse_template,energy,
					pulse_length,record->trigger_size);
			//Update list
			event_list->energies[i] = energy/calfac;
		} else{
			if (i>0 && event_list->grades1[i-1]==-3 && (pulse_time - event_list->event_indexes[i-1])==1){
				current_pulse_duration = -3;
			} else {
				current_pulse_duration = event_list->grades1[i];
			}
			//Compute energy
			energy = (opt_filter_collection->ph_a*event_list->pulse_heights[i]*0.5 + opt_filter_collection->ph_b);
			event_list->energies[i] = energy/calfac;
			//Remove pulse from data
			subtractPulse(record->adc_double,pulse_time,pulse_template,energy,
					pulse_length,record->trigger_size);
		}
		event_list->grades1[i] = current_pulse_duration;
		if(i>0){
			event_list->grades2[i] = pulse_time - event_list->event_indexes[i-1];
		} else {
			event_list->grades2[i] = pulse_length;
		}

		if(identify){
			pulse_real_time = record->time + pulse_time*record->delta_t;
			//Get next impact until it has a chance to match the pulse
			while(impact_time<pulse_real_time-record->delta_t){
				if (!popPhID(record->phid_list,&ph_id,&impact_time,status)){
					break;
				}
			}
			//This impact does match the reconstructed pulse -> add its ph_id
			if(fabs((impact_time-pulse_real_time)) < record->delta_t){
				event_list->ph_ids[i] = ph_id;
				pulse_identified=1;
			} else { //this is a False detection...
				pulse_identified=0;
				event_list->ph_ids[i] = 0;
				//if this was not found as misreconstructed, flag as real false detection
				if (event_list->grades1[i] !=-3){
					event_list->grades1[i] = -4;
				//otherwise flag as correctable false detection
				} else {
					event_list->grades1[i]=-5;
				}
			}
			//Test if this is not a pileup. This test is only valid because we are not simulating
			//the pulse phase wrt the sampling process yet.
			//Will need to refine this once we know the behavior of the pulse reconstruction in this case.
			//This will for instance be more complicated in cases like this:
			//   |  |  |
			//   V  V  V
			//  |    |    |
			//The pulse processing might either associate the first two or the last two (TBC) while for now it is always
			//the first two. Test on energy?
			while (pulse_identified && popPhID(record->phid_list,&ph_id,&impact_time,status) && fabs(impact_time-pulse_real_time) < record->delta_t){
				event_list->ph_ids[i] = -ph_id;
			}
		}
	}
	event_list->index--;
}



/** Wrapper around the whole pulse reconstruction */
void reconstructRecord(TesRecord* record,TesEventList* event_list,ReconstructInit* reconstruct_init,const char identify,int* const status){
	//printf("You lazy guy. Put that in reconstruct init\n");
	int final_length = triggerEvents(record,event_list,reconstruct_init->derivated_template,reconstruct_init->pulse_length-1,
			reconstruct_init->threshold,reconstruct_init->pulse_template_height,reconstruct_init->saturation_value,
			reconstruct_init->derivate_exclusion,reconstruct_init->normal_exclusion,status);
	record->trigger_size = final_length;
	allocateWholeTesEventList(event_list,identify,status);
	computeEnergy(record,event_list,reconstruct_init->opt_filter_collection,reconstruct_init->pulse_template,
			reconstruct_init->pulse_length,reconstruct_init->calfac,identify,status);
}


/** Constructor. Returns a pointer to an empty ReconstructInit data
    structure. */
ReconstructInit* newReconstructInit(int* const status){
	ReconstructInit* reconstruct_init = malloc(sizeof(*reconstruct_init));
	if (NULL==reconstruct_init) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for ReconstructInit failed");
		return(reconstruct_init);
	}

	// Initialize pointers with NULL.
	reconstruct_init->opt_filter_collection =NULL;
	reconstruct_init->pulse_template        =NULL;
	reconstruct_init->derivated_template    =NULL;

	// Initialize values.
	reconstruct_init->pulse_template_height =0;
	reconstruct_init->pulse_length=0;
	reconstruct_init->calfac=0;
	reconstruct_init->threshold=0;

	return(reconstruct_init);
}

/** Destructor. */
void freeReconstructInit(ReconstructInit* reconstruct_init){
	if(NULL!=reconstruct_init){
		freeOptimalFilterCollection(&(reconstruct_init->opt_filter_collection));
		free(reconstruct_init->pulse_template);
		free(reconstruct_init->derivated_template);
	}
	free(reconstruct_init);
	reconstruct_init = NULL;
}

/** Initializes the different variables necessary for the reconstruction */
void initializeReconstruction(ReconstructInit* reconstruct_init,char* const optimal_filter_file,int pulse_length,
		char* const pulse_template_file,double threshold,double calfac,int normal_exclusion,int derivate_exclusion,
		double saturation_value,int* const status){
	// Load OptimalFilterCollection structure
	reconstruct_init->opt_filter_collection = getOptimalFilterCollection(optimal_filter_file, floor(log(pulse_length)/log(2)) - 1,status);
	CHECK_STATUS_VOID(*status);

	//Get pulse template
	//TODO : do that properly from XML file. Add test on sampling rate
	TESProfiles * tes_profiles = newTESProfiles(status);
	readTESProfiles(pulse_template_file,
			"Pulse1",
			tes_profiles,
			status);
	CHECK_STATUS_VOID(*status);

	reconstruct_init->pulse_template = malloc((pulse_length)*sizeof(*(reconstruct_init->pulse_template)));
	if (NULL==reconstruct_init->pulse_template) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for pulse_template failed");
	}
	CHECK_STATUS_VOID(*status);
	memcpy(reconstruct_init->pulse_template,tes_profiles->profiles[0].adc_value[0],pulse_length*sizeof(double));
	destroyTESProfiles(tes_profiles);

	//Derivate the pulse template
	reconstruct_init->derivated_template= malloc((pulse_length-1)*sizeof(*(reconstruct_init->derivated_template)));
	if (NULL==reconstruct_init->derivated_template) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for derivated_pulse failed");
	}
	CHECK_STATUS_VOID(*status);
	derivate_stream(reconstruct_init->pulse_template,reconstruct_init->derivated_template,pulse_length);

	//Get pulse height
	for (int i=0; i< pulse_length-1;i++){
		if(reconstruct_init->derivated_template[i]<reconstruct_init->pulse_template_height){
			break;
		}
		reconstruct_init->pulse_template_height=reconstruct_init->derivated_template[i];
	}

	reconstruct_init->pulse_length=pulse_length;
	reconstruct_init->threshold=threshold;
	reconstruct_init->calfac=calfac;
	reconstruct_init->normal_exclusion=normal_exclusion;
	reconstruct_init->derivate_exclusion=derivate_exclusion;
	reconstruct_init->saturation_value=saturation_value;
}
