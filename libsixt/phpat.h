#ifndef PHPAT_H
#define PHPAT_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "gendet.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phpat(GenDet* const det,
	   const EventFile* const src,
	   EventFile* const dest,
	   const char skip_invalids,
	   int* const status);


#endif /* PHPAT_H */
