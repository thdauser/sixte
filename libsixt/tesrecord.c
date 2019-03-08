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

#include "tesrecord.h"

/** Constructor. Returns a pointer to an empty RecordStruct data
    structure. */
TesRecord* newTesRecord(int* const status){
  CHECK_STATUS_RET(*status,NULL);
  TesRecord* record=malloc(sizeof(TesRecord));
  if (NULL==record) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesRecord failed");
    return(NULL); // yes, we could return record here, but best be obvious
  }
  // Initialize pointers with NULL.
  record->adc_array    = NULL;
  record->adc_double   = NULL;
  record->phid_list    = NULL;
  record->error_double = NULL;

  // Initialize values.
  record->trigger_size =0;
  record->time         =0;
  record->pixid        =0;

  return(record);
}

/** Allocate memory for a RecordStruct data */
void allocateTesRecord(TesRecord * record,unsigned long triggerSize,double delta_t,unsigned char wait_list,int* const status){
	record->trigger_size = triggerSize;
	record->delta_t = delta_t;

	//Allocate adc_array
	record->adc_array = malloc(triggerSize*sizeof(*(record->adc_array)));
	CHECK_NULL_VOID(record->adc_array,*status,"allocateTesRecord: adc_array: memory allocation failed");
	//Allocate adc_double
	record->adc_double = malloc(triggerSize*sizeof(*(record->adc_double)));
	CHECK_NULL_VOID(record->adc_array,*status,"allocateTesRecord: adc_double: memory allocation failed");

	//Allocate error_double
	record->error_double = malloc(triggerSize*sizeof(*(record->error_double)));
	CHECK_NULL_VOID(record->adc_array,*status,"allocateTesRecord: error_double: memory allocation failed");

	//Allocate phid_list
	record->phid_list = newAllocatedPhIDList(MAXIMPACTNUMBER,wait_list,status);
	CHECK_STATUS_VOID(*status);

}

// useful routine to generate a new tes record in one step
TesRecord *createTesRecord(unsigned long triggerSize,double delta_t,unsigned char wait_list,int* const status){
  TesRecord *record=newTesRecord(status);
  CHECK_STATUS_RET(*status,NULL);
  allocateTesRecord(record,triggerSize,delta_t,wait_list,status);
  CHECK_STATUS_RET(*status,NULL);
  return record;
}


/** Change the record length of a tes record **/
void resizeTesRecord(TesRecord *record,unsigned long triggerSize, int* const status) {
  CHECK_STATUS_VOID(*status);

  // reallocate the buffers
  // note: since a tes record has no concept of how full it is,
  // this can destroy data!
  uint16_t* new_adc=realloc(record->adc_array,triggerSize*sizeof(uint16_t));
  CHECK_NULL_VOID(new_adc,*status,"Cannot reallocate ADC array");
  double *new_doub=realloc(record->adc_double,triggerSize*sizeof(double));
  CHECK_NULL_VOID(new_doub,*status,"Cannot reallocate ADC double array");
  double *new_err=realloc(record->error_double,triggerSize*sizeof(double));
  CHECK_NULL_VOID(new_err,*status,"Cannot reallocate ADC double array");

  record->adc_array=new_adc;
  record->adc_double=new_doub;
  record->error_double=new_err;
  record->trigger_size=triggerSize;
}


/** Destructor of the RecordStruct data structure. */
void freeTesRecord(TesRecord** const record){
	if (*record!=NULL){
		free((*record)->adc_double);
		free((*record)->adc_array);
		free((*record)->error_double);
		freePhIDList((*record)->phid_list);
		free(*record);
		*record = NULL;
	}
}

/** Constructor and allocater. Returns a pointer to an allocated PHIDList data
    structure. */
PhIDList* newAllocatedPhIDList(int size,unsigned char wait_list,int* const status){
	PhIDList* phid_list=malloc(sizeof(*phid_list));
	if (NULL==phid_list) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for PhIDList failed");
		return(phid_list);
	}

	// Allocate phid_list.
	phid_list->phid_array = malloc(size*sizeof(*(phid_list->phid_array)));
	if (NULL==phid_list->phid_array) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for PhIDList failed");
		return(phid_list);
	}

	//Allocate time array if this is a wait_list
	if(wait_list){
		phid_list->times = malloc(size*sizeof(*(phid_list->times)));
		if (NULL==phid_list->times) {
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for PhIDList failed");
			return(phid_list);
		}
	}

	// Initialize values.
	phid_list->size  =size;
	phid_list->index =0;
	phid_list->wait_list = wait_list;
	phid_list->n_elements=0;

	return(phid_list);
}

/** Destructor of the RecordStruct data structure. */
void freePhIDList(PhIDList * list){
	if(NULL!=list){
		free(list->phid_array);
		if(list->wait_list){
			free(list->times);
		}
	}
	free(list);
	list=NULL;
}

/** Append ph_id to list */
void appendPhID(PhIDList * const list,long ph_id,double time,int* const status){
	if(list->wait_list){
		if(list->n_elements >= list->size){
			*status=EXIT_FAILURE;
			SIXT_ERROR("Number of impacts in record greater than the maximum allocated number -> abort.\nCheck your count rate or MAXIMPACTNUMBER in testriggerfile.h");
			CHECK_STATUS_VOID(*status);
		} else {
			list->phid_array[list->index % list->size] = ph_id;
			list->times[list->index % list->size] = time;
			list->index++;
			list->n_elements++;
		}
	} else {
		if(list->index >= list->size){
			*status=EXIT_FAILURE;
			SIXT_ERROR("Number of impacts in record greater than the maximum allocated number -> abort.\nCheck your count rate or MAXIMPACTNUMBER in testriggerfile.h");
			CHECK_STATUS_VOID(*status);
		} else {
			list->phid_array[list->index] = ph_id;
			list->index++;
		}
	}
}

/** Get first element in wait list */
int popPhID(PhIDList * const list,long* ph_id,double* time,int* const status){
	if(!(list->wait_list)){
		*status=EXIT_FAILURE;
		SIXT_ERROR("Pop function only defined for wait list");
		CHECK_STATUS_RET(*status,0);
	}
	if(list->n_elements==0){
		return(0);
	}
	*ph_id = list->phid_array[(list->index-list->n_elements) % list->size];
	*time = list->times[(list->index-list->n_elements) % list->size];
	list->n_elements--;
	return(1);
}
