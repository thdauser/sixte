
#ifndef PHOTONLISTFILE_H
#define PHOTONLISTFILE_H 1

#include "sixt.h"


typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;
  
  /** Total number of rows in the FITS file. */
  long nrows;

  /** Number of the current row in the FITS file. Numbering starts at
      1 for the first line. If row is equal to 0, no row has been
      read or written so far. */
  long row;

  /** Column numbers in the FITS binary table. */
  int ctime, cenergy, cra, cdec;

} PhotonListFile;


//////////////////////////////////////////////////////////////////////////


/** Open an existing photon list FITS file and return the
    corresponding PhotonListFile object. The access_mode parameter can
    be either READONLY or READWRITE. */
int openPhotonListFile(PhotonListFile* plf, char* filename, int access_mode);


/** Create new photon list FITS file according to a given template and
    return the PhotonListFile object for the open file.  */
int openNewPhotonListFile(PhotonListFile* plf, char* filename, char* template);


/** Close open PhotonListFile (FITS file). */
int closePhotonListFile(PhotonListFile*) ;


#endif /* PHOTONLISTFILE_H */
