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

#ifndef LADEVENTFILE_H 
#define LADEVENTFILE_H 1

#include "sixt.h"
#include "ladevent.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event file for the GenDet generic detector model. */
typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;

  /** Total number of rows in the file. */
  long nrows;

  /** Current row in the file. */
  long row;

  /** Column numbers. */
  int ctime, csignal, cpanel, cmodule, celement, canode, cph_id, csrc_id;

} LADEventFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty LADEventFile
    data structure. */
LADEventFile* newLADEventFile(int* const status);

/** Destructor. */
void freeLADEventFile(LADEventFile** const file, int* const status);

/** Create and open a new LADEventFile. The new file is generated
    according to the specified template. */
LADEventFile* openNewLADEventFile(const char* const filename,
				  char* const ancrfile,
				  char* const respfile,
				  const double mjdref,
				  const double timezero,
				  const double tstart,
				  const double tstop,
				  const char clobber,
				  int* const status);

/** Open an existing LADEventFile. */
LADEventFile* openLADEventFile(const char* const filename,
			       const int mode, 
			       int* const status);

/** Append a new event to the event file. */
void addLADEvent2File(LADEventFile* const file, 
		      LADEvent* const event, 
		      int* const status);

/** Read the LADEvent at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getLADEventFromFile(const LADEventFile* const file,
			 const int row, 
			 LADEvent* const event,
			 int* const status);

/** Update the LADEvent at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateLADEventInFile(const LADEventFile* const file,
			  const int row, LADEvent* const event,
			  int* const status);


#endif /* LADEVENTFILE_H */
