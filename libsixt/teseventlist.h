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
#include "pixelimpact.h"

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

	//int * event_indexes;	
	double * event_indexes;	//SIRENA

	/** Pulse height of the photons */
	double * pulse_heights;

	/** Average of the first 4 samples of the derivative of the event (pulse) */
	double * avgs_4samplesDerivative;  //SIRENA
	
	/** Low resolution energy estimator (4 samples-long filter) */
	double * Es_lowres;  //SIRENA

	/** Offset relative to the central point of the parabola */
	double * phis;  //SIRENA

	/** Number of samples shifted to find the maximum of the parabola */
	int * lagsShifts;  //SIRENA

	/** Baseline calculated just previously to the pulse (in general)(see 'getB') */
	double * bsln;  //SIRENA

        /** Rms of the baseline calculated just previously to the pulse (in general)(see 'getB') */
	double * rmsbsln;  //SIRENA

	/** Pulse grade */
	int * grading;  //SIRENA

	/** Energy of the photons */
	double * energies;

	/** Grade 1: length of the filter used during the reconstruction */
	int * grades1;

	/** Grade 2: distance in samples to the previous pulse */
	int * grades2;

	/** PH_ID of the reconstructed photons */
	long * ph_ids;

	/** PIX_ID of the reconstructed photons */
	long * pix_ids;

	/** Tstart of the reconstructed photons (in time) */
	double * tstarts;

	/** Tend of the reconstructed photons (in time) */
	double * tends;

	/** Rise time of the reconstructed photons (in time) */
	double * risetimes;

	/** Fall time of the reconstructed photons (in time) */
	double * falltimes;

} TesEventList;

typedef struct {
	/** Pointer to the FITS file. */
	fitsfile* fptr;

	/** Number of the current row in the FITS file. The numbering
	starts at 1 for the first line. If row is equal to 0, no row
	has been read or written so far. */
	long row;

	/** Total number of rows */
	long nrows;

	/** Column numbers for time, energy, grade1, grade2, pixID, RA and DEC columns */
	int timeCol,energyCol,avg_4samplesDerivativeCol,E_lowresCol,grade1Col,grade2Col,phiCol,lagsShiftCol,bslnCol,rmsbslnCol,pixIDCol,riseCol,fallCol,phIDCol,raCol,decCol,detxCol,detyCol,gradingCol,srcIDCol,nxtCol,extCol; //SIRENA

} TesEventFile;

////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////

/** TesEventList constructor. Returns a pointer to an empty TesEventList data
    structure. */
TesEventList* newTesEventList(int* const status);
TesEventList* newTesEventListSIRENA(int* const status);

/** TesEventList Destructor. */
void freeTesEventList(TesEventList* event_list);
void freeTesEventListSIRENA(TesEventList* event_list);

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
TesEventFile* newTesEventFileSIRENA(int* const status);

/** Create and open a new TesEventFile. */
TesEventFile* opennewTesEventFile(const char* const filename,
				  SixtStdKeywords* keywords,
				  const char clobber,
				  int* const status);

TesEventFile* opennewTesEventFileSIRENA(const char* const filename,
				  SixtStdKeywords* keywords,
			          const char* const sirenaVersion,
				  const char clobber,
				  int* const status);

/** Opens a TES event file with the given mode */
TesEventFile* openTesEventFile(const char* const filename,const int mode, int* const status);

/** TesEventFile Destructor. */
void freeTesEventFile(TesEventFile* file, int* const status);

/** Adds the data contained in the event list to the given file */
void saveEventListToFile(TesEventFile* file,TesEventList * event_list,
		double start_time,double delta_t,long pixID,int* const status);
void saveEventListToFileSIRENA(TesEventFile* file,TesEventList * event_list,
		double start_time,double delta_t,long pixID,int* const status);

/** Updates the RA, DEC and DETX/Y columns with the given coordinates */
void updateRaDecDetXY(TesEventFile* file,double ra, double dec, float detx,float dety,int* const status);

/** Add event as reconstructed with the RMF method */
void addRMFImpact(TesEventFile* file,PixImpact * impact,int grade1,int grade2,int grading,int n_xt,double e_xt,int* const status);

/** Adds an event whose signal as not been evaluated yet (necessity in order to keep causality in event file) */
void addEmptyEvent(TesEventFile* file,PixImpact* impact, int* const status);

/** Update signal and grading columns of an event */
//void updateSignal(TesEventFile* file,long row,double energy,double avg_4samplesDerivative,long grade1,long grade2,int grading,int n_xt,double e_xt,int* const status);
void updateSignal(TesEventFile* file,long row,double energy,long grade1,long grade2,int grading,int n_xt,double e_xt,int* const status);

#endif /* TESEVENTLIST_H */
