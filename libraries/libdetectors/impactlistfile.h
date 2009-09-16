#ifndef IMPACTLISTFILE_H
#define IMPACTLISTFILE_H 1

#include "sixt.h"
#include "impact.h"
#include "point.h"


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;
  
  /** Total number of rows in the FITS file. */
  long nrows;

  /** Number of the current row in the FITS file. 
   * Numbering starts at 1 for the first line. 
   * If row is equal to 0, no row has been read or written so far. */
  long row;

  /** Column numbers in the FITS binary table. */
  int ctime, cenergy, cx, cy;

} ImpactListFile;


//////////////////////////////////////////////////////////////////////////


/** Create new impact list FITS file according to a given template and
 * return the ImpactListFile object for the open file. */
int openImpactListFile(ImpactListFile* plf, char* filename, int access_mode);

/** Open an existing photon list FITS file and return the corresponding 
 * ImpactListFile object. */
int openNewImpactListFile(ImpactListFile* plf, char* filename, char* template);

/** Close open ImpactListFile (FITS file). */
int closeImpactListFile(ImpactListFile*) ;

/** Reads the next row from the impact list FITS file. The function increases the internal
 * row counter of the ImpactListFile data structure. E.g. if 'row==0' at the beginning of
 * the function call, the first row from the FITS table is read and the counter is 
 * increased to 'row==1'. */
int getNextImpactListFileRow(ImpactListFile*, Impact*);

/** Checks whether the end of the impact list is reached. 
 * If the internal pointer of the ImpactListFile data structure points to the last line
 * in the file, i.e. this is the formerly read line, or has an even higher value, the 
 * function return value is 1, otherwise it is 0. */
int ImpactListFile_EOF(ImpactListFile*);


#endif /* IMPACTLISTFILE_H */
