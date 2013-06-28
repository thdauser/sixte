#ifndef PHIMG_H
#define PHIMG_H 1

#include "sixt.h"
#include "attitude.h"
#include "check_fov.h"
#include "gentel.h"
#include "photon.h"
#include "impact.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


int phimg(const GenTel* const tel,
	  Attitude* const ac,
	  Photon* const ph,
	  Impact* const imp,
	  int* const status);


#endif /* PHIMG_H */