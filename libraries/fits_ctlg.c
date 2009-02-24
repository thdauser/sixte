#include <malloc.h>
#include "fits_ctlg.h"



/////////////////////////
// SOURCEFILES:

void create_srctbl_parameters(char *ftype[N_SOURCE_FIELDS], char *fform[N_SOURCE_FIELDS], 
			      char *funit[N_SOURCE_FIELDS]) {
  int counter;

  /* determine field types of the table in the FITS file */
  for(counter = 0; counter < N_SOURCE_FIELDS; counter++) {
    ftype[counter] = (char *) malloc(10 * sizeof(char));
    fform[counter] = (char *) malloc(5 * sizeof(char));
    funit[counter] = (char *) malloc(32* sizeof(char));
  }

  /* 1. : right ascension */
  ftype[0] = "R_A_";
  fform[0] = "D";
  funit[0] = "decimal degrees";

  /* 2. : declination */
  ftype[1] = "Dec_";
  fform[1] = "D";
  funit[1] = "decimal degrees";

  /* 3. : source countrate */
  ftype[2] = "src_cps";
  fform[2] = "E";
  funit[2] = "";

  /* 4. : E_min */
  ftype[3] = "E_min";
  fform[3] = "E";
  funit[3] = "";

  /* 5. : E_max */
  ftype[4] = "E_max";
  fform[4] = "E";
  funit[4] = "";

  /* 6. : file containing a source model */
  ftype[5] = "SrcModel";
  fform[5] = "10A";
  funit[5] = "";
}





// writes a row of data into the FITS file
void add_srctbl_row(fitsfile *fptr, long row, const double rasc, const double dec, const float countrate, int *status) 
{
  double dbuffer[1];
  float fbuffer[1];

  // insert new row to binary FITS table
  fits_insert_rows(fptr, row++, 1, status);

  // right ascension
  dbuffer[0] = rasc;
  fits_write_col(fptr, TDOUBLE, 1, row, 1, 1, dbuffer, status);

  // declination
  dbuffer[0] = dec;
  fits_write_col(fptr, TDOUBLE, 2, row, 1, 1, dbuffer, status);

  // source countrate
  fbuffer[0] = countrate;
  fits_write_col(fptr, TFLOAT, 3, row, 1, 1, fbuffer, status);

}



// reads a row of data from the FITS file
int get_srctbl_row(fitsfile *fptr, long row, const int columns[], double *rasc, double *dec, float *countrate,int *status) 
{
  double dbuffer[1];
  float fbuffer[1];
  int anynul;

  // right ascension
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, columns[0], row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  // int fits_read_col / ffgcv
  //    (fitsfile *fptr, int datatype, int colnum, LONGLONG firstrow, LONGLONG firstelem,
  //     LONGLONG nelements, DTYPE *nulval, DTYPE *array, int *anynul, int *status) 
  *rasc = dbuffer[0];

  // declination
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, columns[1], row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *dec = dbuffer[0];

  // countrate
  fbuffer[0] = 0.0;
  fits_read_col(fptr, TFLOAT, columns[2], row+1, 1, 1, fbuffer, fbuffer, &anynul, status);
  *countrate = fbuffer[0];

  return(anynul);
}






////////////////////////////////////////////
// ORBITFILES:


// This function sets the necessary parameters to create a table in the orbit FITS file.
void create_orbtbl_parameter(char *ftype[N_ORBIT_FIELDS], char *fform[N_ORBIT_FIELDS], char *funit[N_ORBIT_FIELDS]) {
  int counter;

  /* determine field types of the table in the FITS file */
  for(counter = 0; counter < N_ORBIT_FIELDS; counter++) {
    ftype[counter] = (char *) malloc(10 * sizeof(char));
    fform[counter] = (char *) malloc(5 * sizeof(char));
    funit[counter] = (char *) malloc(10 * sizeof(char));
  }

  // 1. time
  ftype[0] = "Time";
  fform[0] = "D";
  funit[0] = "s";

  // 2. x coordinate
  ftype[1] = "X";
  fform[1] = "D";
  funit[1] = "km";

  // 3. y coordinate
  ftype[2] = "Y";
  fform[2] = "D";
  funit[2] = "km";

  // 4. z coordinate
  ftype[3] = "Z";
  fform[3] = "D";
  funit[3] = "km";

  // 5. velocity component in x direction
  ftype[4] = "Vx";
  fform[4] = "D";
  funit[4] = "km/s";

  // 6. velocity component in y direction
  ftype[5] = "Vy";
  fform[5] = "D";
  funit[5] = "km/s";

  // 7. velocity component in z direction
  ftype[6] = "Vz";
  fform[6] = "D";
  funit[6] = "km/s";
}




// writes a row of orbit data in to the orbit FITS file
void add_orbtbl_row(fitsfile *fptr, long row, double time, struct vector r, struct vector v, int *status) {
  fits_insert_rows(fptr, row++, 1, status);
  fits_write_col(fptr, TDOUBLE, 1, row, 1, 1, &time, status);
  fits_write_col(fptr, TDOUBLE, 2, row, 1, 1, &r.x, status);
  fits_write_col(fptr, TDOUBLE, 3, row, 1, 1, &r.y, status);
  fits_write_col(fptr, TDOUBLE, 4, row, 1, 1, &r.z, status);
  fits_write_col(fptr, TDOUBLE, 5, row, 1, 1, &v.x, status);
  fits_write_col(fptr, TDOUBLE, 6, row, 1, 1, &v.y, status);
  fits_write_col(fptr, TDOUBLE, 7, row, 1, 1, &v.z, status);
}




// reads a row of data from the orbit FITS file
int get_orbtbl_row(fitsfile *fptr, long row, double *time, struct vector *r, struct vector *v, int *status)
{
  int anynul = 0;
  double dbuffer[1];
  
  // time (1st column)
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 1, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  *time = dbuffer[0];

  // satellite position (columns 2-4)
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 2, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  r->x = dbuffer[0];
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 3, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  r->y = dbuffer[0];
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 4, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  r->z = dbuffer[0];

  // satellite velocity (columns 5-7)
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 5, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  v->x = dbuffer[0];
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 6, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  v->y = dbuffer[0];
  dbuffer[0] = 0.0;
  fits_read_col(fptr, TDOUBLE, 7, row+1, 1, 1, dbuffer, dbuffer, &anynul, status);
  v->z = dbuffer[0];

  return(anynul);
}



