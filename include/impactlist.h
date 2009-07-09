#ifndef IMPACTLIST_H
#define IMPACTLIST_H 1

#include "sixt.h"


/** Structure that contains all information, which is necessary to access
 * an impact list FITS file. */
struct Impactlist_File {
  fitsfile *fptr;

  long row; /**< Current row in the table (starting at 0). */
  long nrows; /**< Total number of rows in the table. */

  /* Column numbers of the individual impact list entries. 
   * The numbers start at 1. The number 0 means, that there 
   * is no corresponding column in the table. */
  int ctime, cenergy, cx, cy;
};


#endif /* IMPACTLIST_H */
