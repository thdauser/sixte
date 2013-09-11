#ifndef PATTERNFILE_H 
#define PATTERNFILE_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "pattern.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event pattern file for the GenDet generic detector model. */
typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;

  /** Total number of rows in the file. */
  long nrows;

  /** Column numbers. */
  int ctime, cframe, cpi, csignal, crawx, crawy, cra, cdec, 
    cph_id, csrc_id, cnpixels, ctype, cpileup, csignals, cpis;

} PatternFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty PatternFile data
    structure. */
PatternFile* newPatternFile(int* const status);

/** Destructor. */
void destroyPatternFile(PatternFile** const file, int* const status);

/** Create and open a new PatternFile. */
PatternFile* openNewPatternFile(const char* const filename,
				char* const telescop,
				char* const instrume,
				char* const filter,
				char* const ancrfile,
				char* const respfile,
				const double mjdref,
				const double timezero,
				const double tstart,
				const double tstop,
				const int nxdim,
				const int nydim,
				const char clobber,
				int* const status);

/** Open an existing PatternFile. */
PatternFile* openPatternFile(const char* const filename,
			     const int mode, int* const status);

/** Append a new pattern to the event file. */
void addPattern2File(PatternFile* const file, 
		     Pattern* const pattern, 
		     int* const status);

/** Read the Pattern at the specified row from the file. The
    numbering for the rows starts at 1 for the first line. */
void getPatternFromFile(const PatternFile* const file,
			const int row, Pattern* const pattern,
			int* const status);

/** Update the Pattern at the specified row in the file. The
    numbering for the rows starts at 1 for the first line. */
void updatePatternInFile(const PatternFile* const file,
			 const int row, Pattern* const pattern,
			 int* const status);

/** Fill the PatternFile with data from the given EventFile. The
    routine is used to transfer events without pattern identification
    to a PatternFile. This is needed, e.g., if no split events are
    simulated. */
void copyEvents2PatternFile(const EventFile* const elf,
			    PatternFile* const plf,
			    int* const status);


#endif /* PATTERNFILE_H */
