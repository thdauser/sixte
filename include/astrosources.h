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


#define POINT_SOURCES    1
#define EXTENDED_SOURCES 2


// Function opens the specified point-source catalog files.
PointSourceFiles* get_PointSourceFiles(int nfiles, char** filename, int* status);
// Close the point.source catalog files.
void free_PointSourceFiles(PointSourceFiles* psf, int* status);

// Functions scans the point-source catalog and returns the sources 
// close to the FOV.
int get_PointSourceCatalog(PointSourceFiles*, PointSourceCatalog**,
			   struct vector normal_vector, const double max_align,
			   struct Spectrum_Store);
				   

#endif /* ASTROSOURCES_H */

