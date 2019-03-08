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


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef IMPACTFILE_H
#define IMPACTFILE_H 1

#include "sixt.h"
#include "impact.h"
#include "point.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


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
  int ctime, cenergy, cx, cy, cph_id, csrc_id;

} ImpactFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty ImpactFile data
    structure. */
ImpactFile* newImpactFile(int* const status);

/** Destructor. */
void freeImpactFile(ImpactFile** const file, int* const status);

/** Open an existing ImpactFile. */
ImpactFile* openImpactFile(const char* const filename,
			   const int mode, int* const status);

/** Create and open a new ImpactFile. The new file is generated
    according to the specified template. */
ImpactFile* openNewImpactFile(const char* const filename,
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

/** Return the next impact from the file. Increments the internal row
    counter by 1 (e.g. if 'row==0' at the beginning of the function
    call, the first row from the FITS table is read and the counter is
    increased to 'row==1'). */
void getNextImpactFromFile(ImpactFile* const file,
			   Impact* const impact,
			   int* const status);

/** Append a new entry to the ImpactFile. */
void addImpact2File(ImpactFile* const ilf,
		    Impact* const impact,
		    int* const status);


#endif /* IMPACTFILE_H */
