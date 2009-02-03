#include "fits_psf.h"



///////////////////////////////////////////////////////////////////////////////////
// creates the necessary data for the FITS-table layout
void psf_create_tbl_parameter(
			      char *ftype[PSF_NFIELDS], 
			      char *fform[PSF_NFIELDS], 
			      char *funit[PSF_NFIELDS],
			      int width      // width of the PSF array in [pixel]
			      ) 
{
    int counter;

    // determine field types of the table in the FITS file
    for(counter = 0; counter < PSF_NFIELDS; counter++) {
      ftype[counter] = (char *) malloc(20 * sizeof(char));
      fform[counter] = (char *) malloc(20 * sizeof(char));
      funit[counter] = (char *) malloc(30 * sizeof(char));
    }

    // off-axis angle
    ftype[0] = "OFFAXANG";
    fform[0] = "D";
    funit[0] = "degrees";

    ftype[1] = "ENERGY";
    fform[1] = "D";
    funit[1] = "keV";

    // x coordinate of PSF rectangle
    ftype[2] = "X";
    fform[2] = "I";
    funit[2] = "";

    ftype[3] = "Y";
    fform[3] = "I";
    funit[3] = "";

    // PSF data
    ftype[4] = "PSF_DATA";
    sprintf(fform[4], "%dD", width*width);
    funit[4] = "";

}



//////////////////////////////////////////////////////////////////////////////////
// writes PSF data to a binary FITS table
int insert_psf_fitsrow(
		       double angle,          // off-axis angle
		       double energy,         // photon energy
		       int x, int y,          // coordinates of PSF sub-rectangle
		       double *data, // PSF data (probabilities within sub-rectangle)
		       long size,    // size of the sub-rectangle (width*height)
		       fitsfile *output_fptr, // FITS file pointer to output file
		       long row               // actual row in the FITS file
		       )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)

    // insert new row to binary FITS table
    if (fits_insert_rows(output_fptr, row-1, 1, &status)) break;

    // insert table data
      
    // write column-entries
    // fits_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
    //         long firstelem, long nelements, void *array, int *status)
    if (fits_write_col(output_fptr, TDOUBLE, 1, row, 1, 1, &angle, &status)) break;
    if (fits_write_col(output_fptr, TDOUBLE, 2, row, 1, 1, &energy, &status)) break;
    if (fits_write_col(output_fptr, TINT, 3, row, 1, 1, &x, &status)) break;
    if (fits_write_col(output_fptr, TINT, 4, row, 1, 1, &y, &status)) break;
    if (fits_write_col(output_fptr, TDOUBLE, 5, row, 1, size, data, &status)) break;

  } while (0);  // end of error loop


  return(status);
}




// read PSF data from a FITS table
int read_psf_fitsrow(
		     double *angle, 
		     double *energy, 
		     int *x, int *y, 
		     double *data, 
		     long size, 
		     fitsfile *fptr, 
		     long row
		     )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)

    int anynul;
    if (fits_read_col(fptr, TDOUBLE, 1, row, 1, 1, angle, angle, &anynul, &status)) 
      break;
    if (fits_read_col(fptr, TDOUBLE, 2, row, 1, 1, energy, energy, &anynul, &status)) 
      break;
    if (fits_read_col(fptr, TINT, 3, row, 1, 1, x, x, &anynul, &status)) break;
    if (fits_read_col(fptr, TINT, 4, row, 1, 1, y, y, &anynul, &status)) break;
    if (fits_read_col(fptr, TDOUBLE, 5, row, 1, size, data, data, &anynul, &status)) 
      break;
    
  } while(0);  // END of error handling loop

  return(status);
}





