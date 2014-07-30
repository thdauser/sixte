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


   Copyright 2014 Philippe Peille, IRAP
*/

#ifndef TESTRIGGERFILE_H
#define TESTRIGGERFILE_H 1

#include "sixt.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"

////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;
  
  /** Total number of rows in the FITS file. */
  long nrows;

  /** Number of the current row in the FITS file. The numbering
      starts at 1 for the first line. If row is equal to 0, no row
      has been read or written so far. */
  long row;

  /** Column numbers for time and trigger columns */
  int timeCol,trigCol;

  /** ID of the pixel corresponding to this file */
  int pixID;

} TesTriggerFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty TesTriggerFile data
    structure. */
TesTriggerFile* newTesTriggerFile(int* const status);

/** Destructor. */
void freeTesTriggerFile(TesTriggerFile** const file, int* const status);

/** Create and open a new TesTriggerFile. */
TesTriggerFile*  opennewTesTriggerFile(const char* const filename,
				  char* const telescop,
				  char* const instrume,
				  char* const filter,
				  char* const ancrfile,
				  char* const respfile,
				  char* const xmlfile,
				  char* const impactlist,
				  const double mjdref,
				  const double timezero,
				  const double tstart,
				  const double tstop,
				  const int pixID,
				  int triggerSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status);

/** Writes the ADC curves in the TES trigger format*/
void writeTriggerFileWithImpact(TESDataStream* const stream,
				char* const tesTriggerFilename,char* const telescop,
				char* const instrume,char* const filter,
				char* const ancrfile,char* const respfile,
				char* const xmlfile,char* const impactlist,
				const double mjdref,const double timezero,
				double tstart,double tstop,const int triggerSize,
				const int preBufferSize,const double sampleFreq,
				const char clobber,const int pixlow,const int Npix,
				int* const status);


#endif /* TESTRIGGERFILE_H */
