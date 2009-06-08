#ifndef EVENTLIST_TYPES_H
#define EVENTLIST_TYPES_H (1)


#include "detectors.enum.h"
#include "fitsio.h"
#include "sixt.h"


/** Event. */
struct Event {
  double time;
  long pha;
  int xi, yi;
  long frame;
  long patnum, patid;
  long pileup;
  double ra, dec;      /**< Right ascension and declination [degree]. */
  long sky_xi, sky_yi; /**< Integer sky pixel coordinates. */
};



/** Structure that contains all information, which is necessary to access
 * an event list FITS file. */
struct Eventlist_File {
  fitsfile *fptr;

  DetectorTypes detectortype;

  int ncolumns; /**< Number of columns in the FITS event list table. */
  long row;     /**< Current row in the table (starting at 0). */
  long nrows;   /**< Total number of rows in the table. */

  /* Column numbers of the individual event list entries. 
   * The numbers start at 1. The number 0 means, that there 
   * is no corresponding column in the table. */
  int ctime, cpha, crawx, crawy, cframe;
  int cra, cdec, cskyx, cskyy;
  int cpatnum, cpatid, cpileup;
};


#endif /* EVENTLIST_TYPES_H */

