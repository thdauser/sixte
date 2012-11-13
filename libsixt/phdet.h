#ifndef PHDET_H
#define PHDET_H 1

#include "sixt.h"
#include "geninst.h"
#include "impact.h"
#include "impactlistfile.h"
#include "eventlistfile.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phdetGenInst(GenInst* const inst,
		  Impact* const impact,
		  const double tend,
		  int* const status);


#endif /* PHDET_H */
