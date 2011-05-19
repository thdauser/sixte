#ifndef PHPROJ_H
#define PHPROJ_H 1

#include "sixt.h"
#include "attitudecatalog.h"
#include "event.h"
#include "eventlistfile.h"
#include "gendet.h"
#include "point.h"
#include "vector.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phproj(GenDet* const det,
	    AttitudeCatalog* const ac,
	    EventListFile* const elf,
	    const double t0,
	    const double exposure,
	    int* const status);


#endif /* PHPROJ_H */
