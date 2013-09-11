#ifndef XMS_PIXTEMP_H
#define XMS_PIXTEMP_H 1

#include "sixt.h"
#include "eventfile.h"
#include "event.h"


#define TOOLSUB xms_pixtemp_main
#include "headas_main.c"


struct Parameters{
  /** Filename of the XMS event file. */
  char EventList[MAXMSG];
  /** Filename of the output file. */
  char OutputFile[MAXMSG];
  /** Filename of the detector response file containing the EBOUDNS
      table. */
  char RSP[MAXMSG];

  /** X- and y- coordinate of the pixel to be investigated. */
  int pixx, pixy;

};


//////////////////////////////////////////////////////////////////


int xms_pixtemp_getpar(struct Parameters*);


#endif /* XMS_PIXTEMP_H */
