#ifndef LADSIGNALFILE_H 
#define LADSIGNALFILE_H 1

#include "sixt.h"
#include "ladsignal.h"


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

} LADSignalFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty LADSignalFile
    data structure. */
LADSignalFile* newLADSignalFile(int* const status);

/** Destructor. */
void freeLADSignalFile(LADSignalFile** const file, int* const status);

/** Create and open a new LADSignalFile. The new file is generated
    according to the specified template. */
LADSignalFile* openNewLADSignalFile(const char* const filename,
				    const char clobber,
				    int* const status);

/** Open an existing LADSignalFile. */
LADSignalFile* openLADSignalFile(const char* const filename,
				 const int mode, 
				 int* const status);

/** Append a new event to the event file. */
void addLADSignal2File(LADSignalFile* const file, 
		       LADSignal* const event, 
		       int* const status);

/** Read the LADSignal at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getLADSignalFromFile(const LADSignalFile* const file,
			  const int row, 
			  LADSignal* const event,
			  int* const status);

/** Update the LADSignal at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateLADSignalInFile(const LADSignalFile* const file,
			   const int row, LADSignal* const event,
			   int* const status);


#endif /* LADSIGNALFILE_H */
