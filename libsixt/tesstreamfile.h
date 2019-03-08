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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef TESSTREAMFILE_H
#define TESSTREAMFILE_H 1

#include "sixt.h"
#include "tesdatastream.h"

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

  /** Number of pixels to handle. */
  int Npix;

  /** Array containing the column numbers for each pixel. */
  int *columnNumber;

  /** Array containing the extension number for each pixel. */
  int *extNumber;

  /** Column number for time. */
  int timeColumnNumber;

} TesStreamFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty TesStreamFile data
    structure. */
TesStreamFile* newTesStreamFile(int* const status);

/** Destructor. */
void freeTesStreamFile(TesStreamFile** const file, int* const status);

/** Open an existing TesStreamFile. */
TesStreamFile* openTesStreamFile(const char* const filename,int pixlow,int pixhigh,
				 const int mode, int* const status);

/** Generate a TESDataStream from tesstreamfile. */
TesStreamFile* opennewTesStreamFile(const char* const filename,
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
				    const char clobber,
				    int* const status);

/** Generate a TESDataStream from tesstreamfile. */
TESDataStream* generateTESDataStreamFromFile(TESDataStream* stream,TesStreamFile* file,double tstart_stream,
					     double tstart,double tstop, double sampleFreq,int* const status);


#endif /* TESSTREAMFILE_H */
