#ifndef PHPROJ_H
#define PHPROJ_H 1

#include "sixt.h"
#include "attitude.h"
#include "pattern.h"
#include "patternfile.h"
#include "geninst.h"
#include "point.h"
#include "vector.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phproj(GenInst* const inst,
	    Attitude* const ac,
	    PatternFile* const plf,
	    const double t0,
	    const double exposure,
	    int* const status);


#endif /* PHPROJ_H */
