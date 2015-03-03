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

#ifndef TESEVENTLIST_H
#define TESEVENTLIST_H 1

#include "sixt.h"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

typedef struct {
	/** Current size of the list */
	int size;

	/** Current size of the energy/grade lists */
	int size_energy;

	/** Current end index of the list */
	int index;

	/** Index arrival time of the photons inside a record */
	int * event_indexes;

	/** Pulse height of the photons */
	double * pulse_heights;

	/** Energy of the photons */
	double * energies;

	/** Grade 1: length of the filter used during the reconstruction */
	int * grades1;

	/** Grade 2: distance in samples to the previous pulse */
	int * grades2;

	/** PH_ID of the reconstructed photons */
	long * ph_ids;

} TesEventList;

typedef struct {
	/** Pointer to the FITS file. */
	fitsfile* fptr;

	/** Number of the current row in the FITS file. The numbering
	starts at 1 for the first line. If row is equal to 0, no row
	has been read or written so far. */
	long row;

	/** Column numbers for time, energy, grade1, grade2, and pixID columns */
	int timeCol,energyCol,grade1Col,grade2Col,pixIDCol,phIDCol;

} TesEventFile;

////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////

/** TesEventList constructor. Returns a pointer to an empty TesEventList data
    structure. */
TesEventList* newTesEventList(int* const status);

/** TesEventList Destructor. */
void freeTesEventList(TesEventList* event_list);

/** Allocates memory for a TesEventList structure for the triggering stage:
 *  only event_index, pulse_height and grades1 */
void allocateTesEventListTrigger(TesEventList* event_list,int size,int* const status);

/** Allocates memory for the energy and grade2 arrays according to
 *  current size if necessary. Only allocate ph_id if specified */
void allocateWholeTesEventList(TesEventList* event_list,unsigned char allocate_ph,int* const status);

/** Appends the index and pulse_height lists to the list */
void addEventToList(TesEventList* event_list,int index,double pulse_height,int grade1,int* const status);


/** TesEventFile constructor. Returns a pointer to an empty TesEventFile data
    structure. */
TesEventFile* newTesEventFile(int* const status);

/** Create and open a new TesEventFile. */
TesEventFile* opennewTesEventFile(const char* const filename,
				  char* const telescop,
				  char* const instrume,
				  char* const filter,
				  char* const ancrfile,
				  char* const respfile,
				  const double mjdref,
				  const double timezero,
				  const double tstart,
				  const double tstop,
				  const char clobber,
				  int* const status);

/** TesEventFile Destructor. */
void freeTesEventFile(TesEventFile* file, int* const status);

/** Adds the data contained in the event list to the given file */
void saveEventListToFile(TesEventFile* file,TesEventList * event_list,
		double start_time,double delta_t,long pixID,int* const status);

#endif /* TESEVENTLIST_H */
