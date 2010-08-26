#ifndef GENDET_H 
#define GENDET_H 1

#include "sixt.h"
#include "sixt_string.h"
#include "gendetline.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic pixelized X-ray detector. The data structure is designed
    in a generic way. The characteristic properties for a particular
    detector are defined in a detector specific XML file. */
typedef struct {

  /** Detector dimensions. Width and height [pixels]. */
  int xwidth, ywidth;

  /** Array of pointers to pixel lines. */
  GenDetLine** line;

  /** Detector response matrix. The RSP file that is originally loaded
      may also contain ARF contributions. But as they already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix. */
  struct RMF* rmf;

} GenDet;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenDet data structure. */
GenDet* newGenDet(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenDet data structure to NULL. */
void destroyGenDet(GenDet** det);


#endif /* GENDET_H */
