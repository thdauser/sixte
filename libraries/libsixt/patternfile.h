#ifndef PATTERNFILE_H 
#define PATTERNFILE_H 1

#include "sixt.h"
#include "pattern.h"
#include "eventlistfile.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event pattern file for the GenDet generic detector model. */
typedef struct {
  
  EventListFile* eventlistfile;

  /** Column numbers. */
  int cpat_type, cpileup, cphas;

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
				int* const status);

/** Open an existing EventListFile. */
PatternFile* openPatternFile(const char* const filename,
			     const int mode, int* const status);

/** Append a new pattern to the event file. */
void addPattern2File(PatternFile* const file, 
		     Pattern* const pattern, 
		     int* const status);


#endif /* PATTERNFILE_H */
