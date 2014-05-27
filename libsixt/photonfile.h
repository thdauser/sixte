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
*/

#ifndef PHOTONFILE_H
#define PHOTONFILE_H 1

#include "sixt.h"
#include "photon.h"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;
  
  /** Total number of rows in the FITS file. */
  long nrows;

  /** Number of the current row in the FITS file. Numbering starts at
      1 for the first line. If row is equal to 0, no row has been
      read or written so far. */
  long row;

  /** Column numbers in the FITS binary table. */
  int ctime, cenergy, cra, cdec, cph_id, csrc_id;

} PhotonFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory and initializes pointers and
    variables. */
PhotonFile* newPhotonFile(int* const status);

/** Close open PhotonFile (FITS file). */
void freePhotonFile(PhotonFile** const plf, int* const status);

/** Open an existing photon list FITS file and return the
    corresponding PhotonFile object. The access_mode parameter can
    be either READONLY or READWRITE. */
PhotonFile* openPhotonFile(const char* const filename, 
			   const int access_mode,
			   int* const status);


/** Create new photon list FITS file according to a given template and
    return the PhotonFile object for the open file. If the clobber
    parameter has a value different from 0, any existing files will be
    overwritten. */
PhotonFile* openNewPhotonFile(const char* const filename, 
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

/** Read the next Photon from the PhotonFile. This routine
    increases the internal counter of the PhotonFile data
    structure. The return value is the error status. */
int PhotonFile_getNextRow(PhotonFile* const plf, Photon* const ph);

/** Read a specific row from the PhotonFile. This routine does NOT
    increase the internal counter of the PhotonFile data
    structure. The return value of the function is the error
    status. */
int PhotonFile_getRow(PhotonFile* const plf, 
		      Photon* const ph, const long row);

/** Append a new photon to the to PhotonFile. */
int addPhoton2File(PhotonFile* const plf, Photon* const ph);


#endif /* PHOTONFILE_H */
