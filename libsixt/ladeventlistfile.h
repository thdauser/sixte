#ifndef LADEVENTLISTFILE_H 
#define LADEVENTLISTFILE_H 1

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

} LADEventListFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty LADEventListFile
    data structure. */
LADEventListFile* newLADEventListFile(int* const status);

/** Destructor. */
void freeLADEventListFile(LADEventListFile** const file, 
			  int* const status);

/** Create and open a new LADEventListFile. The new file is generated
    according to the specified template. */
LADEventListFile* openNewLADEventListFile(const char* const filename,
					  char* const ancrfile,
					  char* const respfile,
					  const double mjdref,
					  const double timezero,
					  const double tstart,
					  const double tstop,
					  const char clobber,
					  int* const status);

/** Open an existing LADEventListFile. */
LADEventListFile* openLADEventListFile(const char* const filename,
				       const int mode, 
				       int* const status);

/** Append a new event to the event file. */
void addLADEvent2File(LADEventListFile* const file, 
		      LADEvent* const event, 
		      int* const status);

/** Read the LADEvent at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getLADEventFromFile(const LADEventListFile* const file,
			 const int row, 
			 LADEvent* const event,
			 int* const status);

/** Update the LADEvent at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateLADEventInFile(const LADEventListFile* const file,
			  const int row, LADEvent* const event,
			  int* const status);


#endif /* LADEVENTLISTFILE_H */
