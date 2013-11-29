#ifndef GENTEL_H 
#define GENTEL_H 1

#include "sixt.h"
#include "arf.h"
#include "phabkg.h"
#include "psf.h"
#include "vignetting.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic X-ray telescope. The characteristic properties for a
    particular telescope are defined in a specific XML file. */
typedef struct {

  /** Detector and telescope ARF containing the effective area. */
  char* arf_filename;
  struct ARF* arf;

  /** Telescope PSF. */
  PSF* psf;

  /** Telescope vignetting function. */
  Vignetting* vignetting;

  /** Focal length of the X-ray telescope [m]. */
  float focal_length;

  /** Diameter of the FoV [rad]. In the XML file the diameter is given
      in [deg], but it is converted to [rad] when parsing the XML
      file. */
  float fov_diameter;

} GenTel;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenTel data structure and
    initializes it with the values from the specified XML definition
    file. */
GenTel* newGenTel(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenTel data structure to NULL. */
void destroyGenTel(GenTel** const tel);


#endif /* GENTEL_H */
