#ifndef XMS_PIXTEMP_H
#define XMS_PIXTEMP_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "xmseventfile.h"
#include "xmsevent.h"


#define TOOLSUB xms_pixtemp_main
#include "headas_main.c"


struct Parameters{
  /** Filename of the XMS event file. */
  char eventlist_filename[MAXMSG];
  /** Filename of the output file. */
  char output_filename[MAXMSG];
  /** Filename of the detector response file containing the EBOUDNS table. */
  char rsp_filename[MAXMSG];

  /** X- and y- coordinate of the pixel to be investigated. */
  int pixx, pixy;

};


//////////////////////////////////////////////////////////////////


int xms_pixtemp_getpar(struct Parameters*);


#endif /* XMS_PIXTEMP_H */
