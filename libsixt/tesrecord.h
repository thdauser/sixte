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

#ifndef TESRECORD_H
#define TESRECORD_H 1

#include "sixt.h"

////////////////////////////////////////////////////////////////////////
// Constants.
////////////////////////////////////////////////////////////////////////

/** Maximum number of impacts per trigger. */
#define MAXIMPACTNUMBER (50000)

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

typedef struct{
	/** Array containing the phIDs */
	long* phid_array;

	/** Array containing the corresponding impact times (only relevant in case of wait list) */
	double* times;

	/** Boolean to state if this should be a wait list */
	unsigned char wait_list;

	/** Number of elements in the wait list */
	int n_elements;

	/** Current index in the list */
	int index;

	/** Size of the list */
	int size;
}PhIDList;

typedef struct{
	/** Number of ADC values in the record */
	unsigned long trigger_size;

	/** Start time of the record */
	double time;

	/** Time difference between two samples */
	double delta_t;

	/** Buffer to read a record of ADC values */
	uint16_t* adc_array;

	/** Double version of the record */
	double* adc_double;

	/** Error signal */
	double* error_double;
    
    /** EXTEND of the record */
	long extend;

	/** PIXID of the record */
	long pixid;

	/** Array of the PH_ID in the record */
	PhIDList* phid_list;

}TesRecord;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////

/** Constructor. Returns a pointer to an empty RecordStruct data
    structure. */
TesRecord* newTesRecord(int* const status);

/** Allocates memory for a RecordStruct data */
void allocateTesRecord(TesRecord * record,unsigned long triggerSize,double delta_t,unsigned char wait_list,int* const status);

TesRecord *createTesRecord(unsigned long triggerSize,double delta_t,unsigned char wait_list,int* const status);

/** resizes the  TesRecord **/
void resizeTesRecord(TesRecord *record,unsigned long triggerSize, int* const status);

/** Destructor of the RecordStruct data structure. */
void freeTesRecord(TesRecord** const record);


/** Constructor and allocater. Returns a pointer to an allocated PHIDList data
    structure. */
PhIDList* newAllocatedPhIDList(int size,unsigned char wait_list,int* const status);

/** Destructor of the RecordStruct data structure. */
void freePhIDList(PhIDList * list);

/** Append ph_id to list */
void appendPhID(PhIDList * const list,long ph_id,double time,int* const status);

/** Get first element in wait list */
int popPhID(PhIDList * const list,long* ph_id,double* time,int* const status);


#endif /* TESRECORD_H */
