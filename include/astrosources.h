/** 
 * Contains all FUNCTION declarations needed to handle source catalogs.
 * The variable and TYPE definitions can be found in 'astrosources.types.h'.
 */

#ifndef ASTROSOURCES_H
#define ASTROSOURCES_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "fits_pha.h"
#include "fits_ctlg.h"
#include "spectrum.h"
#include "detectors.h"
#include "psf.h"
#include "vector.h"


#include "astrosources.types.h"


/** Category of input sources: 1=Point sources, 2=Extended Source Images */
typedef enum {
  POINT_SOURCES    =1,
  SOURCE_IMAGES    =3
} SourceCategory;

				   

#endif /* ASTROSOURCES_H */

