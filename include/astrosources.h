/** 
 * Contains basic definitions for astronomical x-ray sources.
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


/** Category of input sources: 1=Point sources, 2=Extended Source Images */
typedef enum {
  POINT_SOURCES    =1,
  SOURCE_IMAGES    =3
} SourceCategory;

				   

#endif /* ASTROSOURCES_H */

