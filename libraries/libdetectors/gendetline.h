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
GenDetLine* newGenDetLine(const int xwidth, int* const status);

/** Destructor. Releases the allocated memory and sets the pointer to
    the GenDetLine data structure to NULL. */
void destroyGenDetLine(GenDetLine** line);

/** Clear all pixels in the GenDetLine. If the anycharge flag of the
    GenDetLine is set to 0, the clearing will be skipped, because in
    this cause all of the pixels should already be set to 0 charge. */
void clearGenDetLine(GenDetLine* const line);

/** Add the charges in line 1 to line 0. The line must have the same
    width. */
void addGenDetLine(GenDetLine* const line0, GenDetLine* const line1);

/** Switch the 2 specified lines. */
void switchGenDetLines(GenDetLine** const line0, GenDetLine** const line1);

/** Read-out a charge from the specified line. As long as there are
    any charges, the return value is 1. If the line contains no more
    charges, the function return value is 0. */
int readoutGenDetLine(GenDetLine* const line, float* charge, int* x);

/** Add a charge (photon energy [keV]) to a particular pixel in the
    specified GenDetLine. The routine sets the anycharge flag of the
    affected line. */
inline void addGenDetCharge2Pixel(GenDetLine* const line, const int column, 
				  float energy);


#endif /* GENDETLINE_H */
