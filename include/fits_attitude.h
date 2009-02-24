#ifndef FITS_ATTITUDE_H
#define FITS_ATTITUDE_H

#include <stdlib.h>
#include <malloc.h>

#include "fitsio.h"

#define N_ATTITUDE_FIELDS 6  // number of fields in the attitude table


// writes a row of attitude data in to the FITS file
int add_attitudetbl_row(fitsfile *fptr, long row, char valtime[], 
			double time, double view_ra, double view_dec, 
			double rollangle, double aspangle, int fitsstatus);

// creates the necessary parameters to generate the table in the attitude FITS file
void create_attitudetbl_parameter(char *ftype[N_ATTITUDE_FIELDS], 
				  char *fform[N_ATTITUDE_FIELDS], 
				  char *funit[N_ATTITUDE_FIELDS]);

// reads a row of data from the attitude FITS file
int get_atttbl_row(fitsfile *fptr, long row, char valtime[], double *time, 
		   double *view_ra, double *view_dec, double *rollangle, 
		   int *status);


#endif
