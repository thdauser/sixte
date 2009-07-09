#ifndef IMPACTLIST_H
#define IMPACTLIST_H 1

#include "sixt.h"
#include "point.h"


/** Structure that contains all information, which is necessary to access
 * an impact list FITS file. */
struct ImpactlistFile {
  fitsfile *fptr;

  long row; /**< Current row in the table (starting at 0). */
  long nrows; /**< Total number of rows in the table. */

  /* Column numbers of the individual impact list entries. 
   * The numbers start at 1. The number 0 means, that there 
   * is no corresponding column in the table. */
  int ctime, cenergy, cx, cy;

  /** Filename of the attitude file that was used for the generation of the impact list. */
  char attitude_filename[MAXMSG];
};


/** Impact of a photon on the detector plane. */
struct Impact {
  double time;
  float energy;
  struct Point2d position;
};


////////////////////////////////////////////////////////////////////////////////////


/** Open an existing impact list FITS file. The access mode can be specified by using
 * the standard CFITSIO access constants READONLY or READWRITE. */
int impactlist_openFile(struct ImpactlistFile* imf, char* filename, int access_mode);

/** Reads the next row from the impact list FITS file. The function increases the internal
 * row counter of the ImpactlistFile data structure. E.g. if 'row==0' at the beginning of
 * the function call, the first row from the FITS table is read and the counter is 
 * increased to 'row==1'. */
int impactlist_getNextRow(struct ImpactlistFile* imf, struct Impact* impact);


#endif /* IMPACTLIST_H */
