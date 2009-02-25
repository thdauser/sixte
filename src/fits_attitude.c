#include "fits_attitude.h"


// creates the necessary parameters to generate the attitude table in the FITS file
void create_attitudetbl_parameter(char *ftype[N_ATTITUDE_FIELDS], 
				  char *fform[N_ATTITUDE_FIELDS], 
				  char *funit[N_ATTITUDE_FIELDS]) 
{
  int counter;

  /* determine field types of the table in the FITS file */
  for(counter = 0; counter < N_ATTITUDE_FIELDS; counter++) {
    ftype[counter] = (char *) malloc(10 * sizeof(char));
    fform[counter] = (char *) malloc(5 * sizeof(char));
    funit[counter] = (char *) malloc(21 * sizeof(char));
  }

  // 1. time in human readable format
  ftype[0] = "VALTIME";
  fform[0] = "19A";
  funit[0] = "yyyy-mm-ddThh:mm:ss";

  // 2. time in MJD format
  ftype[1] = "Time";
  fform[1] = "D";
  funit[1] = "s";

  // 3. right ascension of the viewing direction 
  ftype[2] = "VIEWRA";
  fform[2] = "D";
  funit[2] = "degree";

  // 4. declination of the viewing direction
  ftype[3] = "VIEWDECL";
  fform[3] = "D";
  funit[3] = "degree";

  // 5. roll angle
  ftype[4] = "ROLLANG";
  fform[4] = "D";
  funit[4] = "degree";

  // 6. solar aspect angle
  ftype[5] = "ASPANGLE";
  fform[5] = "D";
  funit[5] = "degree";
}



// writes a row of attitude data in to the FITS file
int add_attitudetbl_row(fitsfile *fptr, long row, char valtime[], double time, double view_ra, double view_dec, double rollangle, double aspangle, int fitsstatus)
{
  int status = fitsstatus;

  fits_insert_rows(fptr, row++, 1, &status);
  fits_write_col(fptr, TSTRING, 1, row, 1, 1, &valtime, &status);
  fits_write_col(fptr, TDOUBLE, 2, row, 1, 1, &time, &status);
  fits_write_col(fptr, TDOUBLE, 3, row, 1, 1, &view_ra, &status);
  fits_write_col(fptr, TDOUBLE, 4, row, 1, 1, &view_dec, &status);
  fits_write_col(fptr, TDOUBLE, 5, row, 1, 1, &rollangle, &status);
  fits_write_col(fptr, TDOUBLE, 6, row, 1, 1, &aspangle, &status);

  return(status);
}




// reads a row of attitude data from the FITS file
int get_atttbl_row(fitsfile *fptr, long row, char valtime[], double *time, 
		   double *view_ra, double *view_dec, double *rollangle, int *status)
{
  int anynul = 0;
  double dbuffer[1];
  
  // TODO implement valtime
  valtime = "";

  // time
  dbuffer[0] = 0.;
  fits_read_col(fptr, TDOUBLE, 2, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *time = dbuffer[0];

  // right ascension of telescope direction
  dbuffer[0] = 0.;
  fits_read_col(fptr, TDOUBLE, 3, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *view_ra = dbuffer[0];

  // declination of telescope direction
  dbuffer[0] = 0.;
  fits_read_col(fptr, TDOUBLE, 4, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *view_dec = dbuffer[0];
  
  // roll angle
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 5, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *rollangle = dbuffer[0];

  return(anynul);
}

