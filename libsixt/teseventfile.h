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

#ifndef TESEVENTFILE_H
#define TESEVENTFILE_H 1

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

  /** Column numbers for time and event columns */
  int timeCol,evtCol;

  /** ID of the pixel corresponding to this file */
  int pixID;

} TesEventFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty TesEventFile data
    structure. */
TesEventFile* newTesEventFile(int* const status);

/** Destructor. */
void freeTesEventFile(TesEventFile** const file, int* const status);

/** Create and open a new TesEventFile. */
TesEventFile*  opennewTesEventFile(const char* const filename,
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
				  int eventSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status);

/** Writes the ADC curves in the TES event format*/
void writeEvents2FITS(TesEventFile** outputFiles,TESDataStream* stream,PixImpFile* impfile,
		      int pixlow,int Npix,double tstart,double tstartTES,
		      double tstop,double sampleFreq,int eventSize,int preBufferSize,
		      float monoen,int* const status);


#endif /* TESEVENTFILE_H */
