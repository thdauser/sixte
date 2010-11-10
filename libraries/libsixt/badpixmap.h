#ifndef BADPIXMAP_H 
#define BADPIXMAP_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/*typedef enum {
  BADPIX_NONE        =  0,
  BADPIX_INSENSITIVE = -1,
  BADPIX_COLD        = -2,
  BADPIX_HOT         =  1
  // TODO Add type for flickering pixels.
  } BadPixType;*/


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

/** Apply the bad pixel map on the detector pixel array. The second
    argument has to be a pointer to the function that is called for
    bad pixel and the 3rd argument is handled to this function as a
    parameter. The function 'func' has the following syntax: void
    func(void* const data, const int x, const int y, const float
    value) with 'x' and 'y' the coordinates of the bad pixel and
    'value' the value in the bad pixel map image. */
void applyBadPixMap(const BadPixMap* const map, const double timespan,
		    void (*encounter) (void* const data, 
				       const int x, const int y, 
				       const float value),
		    void* const data);


#endif /* BADPIXMAP_H */
