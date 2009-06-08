#ifndef ATTITUDEFILE_H
#define ATTITUDEFILE_H 1

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"


#define N_ATTITUDE_FIELDS 6  // number of fields in the attitude table



typedef struct {
  fitsfile* fptr; /**< File pointer to the FITS file. */

  long row; /**< Current row in the attitude FITS table (starting at 0). */
  long nrows; /**< Number of rows in the attitude table. */

  /* Column numbers of the individual attitude file entries. 
   * The numbers start at 1. The number 0 means, that there 
   * is no corresponding column in the table. */
  int ctime, cviewra, cviewdec, crollang;

} AttitudeFile;


typedef struct {
  double time; /**< Time for which the AttitudeFileEntry is valid. */
  double viewra; /**< Right ascension of telescope pointing direction [deg]. */
  double viewdec; /**< Declination of telescope pointing direction [deg].*/
  double rollang; /**< Rollangle [deg]. */
} AttitudeFileEntry;



///////////////////////////////////////////////////////////////////////////


/** Reads a line of data from the attitude table in a FITS file. 
 * The routine does NOT increment the row counter of the AttitudeFile object. */
AttitudeFileEntry read_AttitudeFileEntry(AttitudeFile* af, int* status);

/** Opens an existing attitude file. 
 * The access_mode parameter can be either READONLY or READWRITE. */
AttitudeFile* open_AttitudeFile(const char filename[], int access_mode, int* status);



/////////////////////
// Old routines: 
// writes a row of attitude data in to the FITS file
int add_attitudetbl_row(fitsfile *fptr, long row, char valtime[], 
			double time, double view_ra, double view_dec, 
			double rollangle, double aspangle, int fitsstatus);

// creates the necessary parameters to generate the table in the attitude FITS file
void create_attitudetbl_parameter(char *ftype[N_ATTITUDE_FIELDS], 
				  char *fform[N_ATTITUDE_FIELDS], 
				  char *funit[N_ATTITUDE_FIELDS]);


#endif /* ATTITUDEFILE_H */
