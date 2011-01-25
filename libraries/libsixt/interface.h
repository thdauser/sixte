#ifndef INTERFACE_H
#define INTERFACE_H 1

#include "sixt.h"
#include "attitudecatalog.h"
#include "gendet.h"
#include "photon.h"
#include "photonlistfile.h"
#include "xraysourcecatalog.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void photon_generation(const char* const xml_filename,
		       const char* const attitude_filename,
		       const char* const simput_filename,
		       const char* const photon_filename,
		       const double t0, const double t1,
		       int* const status);


#endif /* INTERFACE_H */
