#ifndef PHOGEN_H
#define PHOGEN_H 1

#include "sixt.h"

#include "attitudecatalog.h"
#include "gendet.h"
#include "phgen.h"
#include "photonlistfile.h"
#include "sourceimage.h"
#include "vector.h"
#include "xraysourcecatalog.h"

#define TOOLSUB phogen_main
#include "headas_main.c"


struct Parameters {
  char xml_filename[MAXFILENAME];
  char attitude_filename[MAXFILENAME];
  char simput_filename[MAXFILENAME];
  char photonlist_filename[MAXFILENAME];
  char photonlist_template[MAXFILENAME];

  double t0, timespan;
};


int phogen_getpar(struct Parameters* parameters);


#endif /* PHOGEN_H */

