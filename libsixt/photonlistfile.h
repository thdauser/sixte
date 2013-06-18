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
  int ctime, cenergy, cra, cdec, cph_id, csrc_id;

} PhotonListFile;


////////////////////////////////////////////////////////////////////////
// Function declarations.
////////////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory and initializes pointers and
    variables. */
PhotonListFile* newPhotonListFile(int* const status);

/** Close open PhotonListFile (FITS file). */
void freePhotonListFile(PhotonListFile** const plf, int* const status);

/** Open an existing photon list FITS file and return the
    corresponding PhotonListFile object. The access_mode parameter can
    be either READONLY or READWRITE. */
PhotonListFile* openPhotonListFile(const char* const filename, 
				   const int access_mode,
				   int* const status);


/** Create new photon list FITS file according to a given template and
    return the PhotonListFile object for the open file. If the clobber
    parameter has a value different from 0, any existing files will be
    overwritten. */
PhotonListFile* openNewPhotonListFile(const char* const filename, 
				      char* const telescop,
				      char* const instrume,
				      char* const filter,
				      const double mjdref,
				      const double timezero,
				      const double tstart,
				      const double tstop,
				      const char clobber,
				      int* const status);

/** Read the next Photon from the PhotonListFile. This routine
    increases the internal counter of the PhotonListFile data
    structure. The return value is the error status. */
int PhotonListFile_getNextRow(PhotonListFile* const plf, Photon* const ph);

/** Read a specific row from the PhotonListFile. This routine does NOT
    increase the internal counter of the PhotonListFile data
    structure. The return value of the function is the error
    status. */
int PhotonListFile_getRow(PhotonListFile* const plf, 
			  Photon* const ph, const long row);

/** Append a new photon to the to PhotonListFile. */
int addPhoton2File(PhotonListFile* const plf, Photon* const ph);


#endif /* PHOTONLISTFILE_H */
