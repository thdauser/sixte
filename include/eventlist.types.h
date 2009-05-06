#ifndef EVENTLIST_TYPES_H
#define EVENTLIST_TYPES_H (1)


#include "fitsio.h"
#include "sixt.h"


// Number of fields in an event list table:
// TIME, PHA, GRADE, RAWX/COLUMN, RAWY/ROW,
// FRAME, PATNUM, PATID
#define N_EVENT_FIELDS 13


// Structure that contains all information, which is necessary to access
// an event list FITS file.
struct Eventlist_File {
  char filename[FILENAME_LENGTH];
  fitsfile *fptr;

  long row;               // current row in the table (starting at 0)
  long nrows;             // total number of rows in the table
};




//
struct Event {
  double time;
  long pha;
  int grade;
  int xi, yi;
  long frame;
  long patnum, patid;
  long pileup;
  double ra, dec;      // right ascension and declination [degree]
  long sky_xi, sky_yi; // integer sky pixel coordinates
};



#endif /* EVENTLIST_TYPES_H */

