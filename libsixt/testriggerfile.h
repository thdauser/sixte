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

#ifndef TESTRIGGERFILE_H
#define TESTRIGGERFILE_H 1

#include "sixt.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"
#include "tesrecord.h"
#include "optimalfilters.h"
#include "tesinitialization.h"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////

typedef struct{
	/** Pointer to the FITS file. */
	fitsfile* fptr;

	/** Total number of rows in the FITS file. */
	long nrows;

	/** Number of the current row in the FITS file. The numbering
      starts at 1 for the first line. If row is equal to 0, no row
      has been read or written so far. */
	long row;

	/** Number of ADC values per record */
	int trigger_size;

	/** Time interval between two ADC values in the file */
	double delta_t;

	/** Column numbers for time, trigger, impact, and pixID columns */
	int timeCol,trigCol,ph_idCol,pixIDCol;

}TesTriggerFile;

////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty TesTriggerFile data
    structure. */
TesTriggerFile* newTesTriggerFile(int triggerSize,int* const status);

/** Destructor. */
void freeTesTriggerFile(TesTriggerFile** const file, int* const status);

/** Create and open a new TesTriggerFile. */
TesTriggerFile* opennewTesTriggerFile(const char* const filename,
				  SixtStdKeywords * keywords,
				  char* const xmlfile,
				  char* const impactlist,
				  int triggerSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status);

/** Create and open a new TesTriggerFile. */
TesTriggerFile* openexistingTesTriggerFile(const char* const filename,SixtStdKeywords* keywords,int* const status);

/** Save pixels, NES/NET and monoen keywords to the given FITS file */
void saveTriggerKeywords(fitsfile* fptr,int firstpix,int lastpix,int numberpix,float monoen,
		int* const numberSimulated,int* const numberTrigger,int* const status);

/** Writes the ADC curves in the TES trigger format*/
void triggerWithImpact(TESDataStream* const stream,TESGeneralParameters * par,
		TESInitStruct* init,float monoen,const char write_file,const char reconstruct,
		ReconstructInit* reconstruct_init,char* const tes_event_file,int event_list_size,
		const char identify,int* const status);


/** Populates a TesRecord structure with the next record */
int getNextRecord(TesTriggerFile* const file,TesRecord* record,int* const status);

/** Writes a record to a file */
void writeRecord(TesTriggerFile* outputFile,TesRecord* record,int* const status);

#endif /* TESTRIGGERFILE_H */
