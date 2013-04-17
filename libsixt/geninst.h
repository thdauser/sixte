#ifndef GENINST_H 
#define GENINST_H 1

#include "sixt.h"
#include "gendet.h"
#include "gentel.h"
#include "eventlistfile.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic X-ray instrument with a pixelized detector. The data
    structure is designed in a generic way in order to provide the
    capability to model different instruments. The characteristic
    properties for a particular instruments are defined in a specific
    XML file. */
typedef struct {

  /** Telescope data structure. */
  GenTel* tel;

  /** Detector data structure. */
  GenDet* det;

  /** File name (without path contributions) of the FITS file
      containing the XML detector definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;

  /** FITS header keyword TELESCOP. */
  char* telescop;

  /** FITS header keyword INSTRUME. */
  char* instrume;

} GenInst;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenInst data structure and
    initializes it with the values from the specified XML definition
    file. The second parameter determines, whether the cosmic ray
    detector background model should be activated. */
GenInst* newGenInst(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenInst data structure to NULL. */
void destroyGenInst(GenInst** const det, int* const status);

/** Parse the GenInst definition from an XML file. */
GenInst* loadGenInst(const char* const filename, int* const status);


#endif /* GENINST_H */
