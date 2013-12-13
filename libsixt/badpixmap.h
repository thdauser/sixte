#ifndef BADPIXMAP_H 
#define BADPIXMAP_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Bad pixel map. */
typedef struct {

  /** Dimensions of the bad pixel map. */
  int xwidth, ywidth;

  /** 2-dimensional pixel array containing the bad pixel map. */
  float** pixels;

  /** This flag specifies if the column contains any bad pixels (value
      1). If not (value 0), the bad pixel application routine can jump
      to the next column. */
  int* anybadpix;

} BadPixMap;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new empty BadPixMap data
    structure. */
BadPixMap* newBadPixMap(int* const status);

/** Allocate memory for a new BadPixMap data structure via the
    constructor and load the bad pixel data from the specified
    file. */
BadPixMap* loadBadPixMap(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the BadPixMap data structure to NULL. */
void destroyBadPixMap(BadPixMap** const map);


#endif /* BADPIXMAP_H */
