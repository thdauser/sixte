#ifndef GENDETLINE_H 
#define GENDETLINE_H 1

#include "sixt.h"
#include "event.h"


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

  /** Photon IDs corresponding to the charges in the individual
      pixels. */
  long** ph_id;

  /** Source IDs corresponding to the charges in the individual
      pixels. */
  long** src_id;

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
void destroyGenDetLine(GenDetLine** const line);

/** Clear all pixels in the GenDetLine. If the anycharge flag of the
    GenDetLine is set to 0, the clearing will be skipped, because in
    this cause all of the pixels should already be set to 0 charge. */
void clearGenDetLine(GenDetLine* const line);

/** Add the charges in line 1 to line 0. The line must have the same
    width. The line 1 is not modified, i.e. the contained charges
    remain in there and have to be cleared separately. */
void addGenDetLine(GenDetLine* const line0, const GenDetLine* const line1);

/** Add a charge (photon energy [keV]) to a particular pixel in the
    specified GenDetLine. The routine sets the anycharge flag of the
    affected line. */
void addGenDetCharge2Pixel(GenDetLine* const line, const int column, 
			   const float energy,
			   const long ph_id, const long src_id);

#endif /* GENDETLINE_H */
