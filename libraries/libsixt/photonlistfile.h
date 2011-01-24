#ifndef PHOTONLISTFILE_H
#define PHOTONLISTFILE_H 1

#include "sixt.h"
#include "photon.h"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


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


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Open an existing photon list FITS file and return the
    corresponding PhotonListFile object. The access_mode parameter can
    be either READONLY or READWRITE. */
int openPhotonListFile(PhotonListFile* const plf, 
		       const char* const filename, 
		       const int access_mode);


/** Create new photon list FITS file according to a given template and
    return the PhotonListFile object for the open file.  */
int openNewPhotonListFile(PhotonListFile* const plf, 
			  const char* const filename, 
			  const char* const template);


/** Close open PhotonListFile (FITS file). */
int closePhotonListFile(PhotonListFile* plf);

/** Read the next Photon from the PhotonListFile. This routine
    increases the internal counter of the PhotonListFile data
    structure. The return value is the error status. */
int PhotonListFile_getNextRow(PhotonListFile* plf, Photon* ph);

/** Read a specific row from the PhotonListFile. This routine does NOT
    increase the internal counter of the PhotonListFile data
    structure. The return value of the function is the error
    status. */
int PhotonListFile_getRow(PhotonListFile* plf, Photon* ph, long row);

/** Append a new photon to the to PhotonListFile. */
int addPhoton2File(PhotonListFile* plf, Photon* ph);


#endif /* PHOTONLISTFILE_H */
