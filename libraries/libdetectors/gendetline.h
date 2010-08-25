#ifndef GENDETLINE_H 
#define GENDETLINE_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** One line in the GenDet model for a generic pixelized X-ray
    detector. */
typedef struct {
  
  /** Number of pixels in this line. */
  int xwidth;

  /** Charges contained in the individual pixels of this line. */
  float* charge;

  /** This flag specifies if the line contains any charges (value
      1). If not (value 0), the read-out does not have to be
      performed. */
  int anycharge;
  
} GenDetLine;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to the newly allocated and cleared
    GenDetLine data structure. */
GenDetLine* newGenDetLine(int xwidth, int* status);

/** Destructor. Releases the allocated memory and sets the pointer to
    the GenDetLine data structure to NULL. */
void destroyGenDetLine(GenDetLine** line);

/** Clear all pixels in the GenDetLine. If the anycharge flag of the
    GenDetLine is set to 0, the clearing will be skipped, because in
    this cause all of the pixels should already be set to 0 charge. */
void clearGenDetLine(GenDetLine* line);


#endif /* GENDETLINE_H */
