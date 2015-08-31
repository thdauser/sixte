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


   Copyright 2014 Thorsten Brand, FAU
*/

#ifndef PIXIMPFILE_H
#define PIXIMPFILE_H 1

#include "sixt.h"
#include "pixelimpact.h"
#include "point.h"

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

  /** Column numbers in the FITS binary table. */
  int ctime, cenergy, cx, cy, cph_id, csrc_id, cpix_id, cu, cv,cgrade1,cgrade2,ctotalenergy;

} PixImpFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty PixImpFile data
    structure. */
PixImpFile* newPixImpFile(int* const status);

/** Destructor. */
void freePixImpFile(PixImpFile** const file, int* const status);

/** Open an existing PixImpFile. */
PixImpFile* openPixImpFile(const char* const filename,
			   const int mode, int* const status);

/** Create and open a new PixImpFile. The new file is generated
    according to the specified template. */
PixImpFile* openNewPixImpFile(const char* const filename,
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
			      const char clobber,
			      int* const status);

/** Return the next pixel impact from the file. Increments the internal row
    counter by 1 (e.g. if 'row==0' at the beginning of the function
    call, the first row from the FITS table is read and the counter is
    increased to 'row==1'). */
int getNextImpactFromPixImpFile(PixImpFile* const file, 
			   PixImpact* const impact,
			   int* const status);

/** Append a new entry to the PixImpFile. */
void addImpact2PixImpFile(PixImpFile* const ilf, 
			  PixImpact* const impact, 
			  int* const status);

/** Get the tstart and tstop times from the file. */
void getPixImpFileTimeValues(PixImpFile* const ilf,
			  double *mjdref,
			  double *timezero,
			  double *tstart,
			  double *tstop,
			  int* const status);

/** update the grading in the pixImpactFile **/
void updateGradingPixImp(PixImpFile* const ilf,
		  int row, long grade1, long grade2,double totalenergy,
		  int* const status);

#endif /* PIXIMPFILE_H */
