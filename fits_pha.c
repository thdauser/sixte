#include "fits_pha.h"


////////////////////////////////////////////////////////////////////////////////
// Spectra

// creates the necessary data for the FITS-table layout
void spectrum_create_tbl_parameter(
				   char *ftype[SPECTRUM_NFIELDS], 
				   char *fform[SPECTRUM_NFIELDS], 
				   char *funit[SPECTRUM_NFIELDS]
				   ) 
{
  int counter;

    // determine field types of the table in the FITS file
    for(counter = 0; counter < SPECTRUM_NFIELDS; counter++) {
      ftype[counter] = (char *) malloc(20 * sizeof(char));
      fform[counter] = (char *) malloc(20 * sizeof(char));
      funit[counter] = (char *) malloc(20 * sizeof(char));
    }

    // PHA channel
    ftype[0] = "CHANNEL";
    fform[0] = "J";
    funit[0] = "";

    // lower bound of energy bin
    ftype[1] = "COUNTS";
    fform[1] = "E";
    funit[1] = "";
}


// writes spectrum to a binary FITS table
int insert_spectrum_fitsrow(
			    long channel,         // PHA channel
			    float p,
			    fitsfile *fptr,       // FITS file pointer to output file
			    long row              // actual row in the FITS file
			    )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)

    // insert new row to binary FITS table
    if (fits_insert_rows(fptr, row-1, 1, &status)) break;

    // insert data to binary table
      
    // write column-entries
    // fits_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
    //         long firstelem, long nelements, void *array, int *status)
    if (fits_write_col(fptr, TLONG, 1, row, 1, 1, &channel, &status)) break;
    if (fits_write_col(fptr, TFLOAT, 2, row, 1, 1, &p, &status)) break;

  } while (0);  // end of error loop

  return(status);
}




// read spectrum line from a FITS table (PHA channel)
int read_spec_fitsrow(
		      long *channel,
		      float *probability,
		      fitsfile *fptr,
		      long row
		      )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)
    int anynul;
    if (fits_read_col(fptr, TLONG, 1, row, 1, 1, channel, channel, &anynul, &status)) break;
    if (fits_read_col(fptr, TFLOAT, 2, row, 1, 1, probability, probability, &anynul, &status)) break;
  } while(0);  // end of error handling loop

  return(status);
}



////////////////////////////////////////////////////////////////////
// Detector Redistribution Matrix

// creates the necessary data for the FITS-table layout
void rmf_create_tbl_parameter(
				  char *ftype[PHA_NFIELDS], 
				  char *fform[PHA_NFIELDS], 
				  char *funit[PHA_NFIELDS]
				  ) 
{
  int counter;

    // determine field types of the table in the FITS file
    for(counter = 0; counter < PHA_NFIELDS; counter++) {
      ftype[counter] = (char *) malloc(20 * sizeof(char));
      fform[counter] = (char *) malloc(20 * sizeof(char));
      funit[counter] = (char *) malloc(20 * sizeof(char));
    }

    // PHA channel
    ftype[0] = "CHANNEL";
    fform[0] = "J";
    funit[0] = "";

    // lower bound of energy bin
    ftype[1] = "E_MIN";
    fform[1] = "E";
    funit[1] = "keV";

    // upper bound of energy bin
    ftype[2] = "E_MAX";
    fform[2] = "E";
    funit[2] = "keV";
}




// writes detector response data to a binary FITS table
int insert_rmf_fitsrow(
			   long channel,         // PHA channel
			   float Emin,
			   float Emax,
			   fitsfile *fptr,       // FITS file pointer to output file
		           long row              // actual row in the FITS file
			   )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)

    // insert new row to binary FITS table
    if (fits_insert_rows(fptr, row-1, 1, &status)) break;

    // insert data to binary table
      
    // write column-entries
    // fits_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
    //         long firstelem, long nelements, void *array, int *status)
    if (fits_write_col(fptr, TLONG, 1, row, 1, 1, &channel, &status)) break;
    if (fits_write_col(fptr, TFLOAT, 2, row, 1, 1, &Emin, &status)) break;
    if (fits_write_col(fptr, TFLOAT, 3, row, 1, 1, &Emax, &status)) break;

  } while (0);  // end of error loop

  return(status);
}



// read detector response data line from a FITS table (one PHA channel)
int read_rmf_fitsrow(
		     float *Emin, 
		     float *Emax, 
		     int *Ngrp, 
		     int *Fchan, 
		     int *Nchan, 
		     float *matrix, 
		     fitsfile *fptr, 
		     long row
		     )
{
  int count;
  int status=EXIT_SUCCESS;  // error status


  do {  // beginning of error handling loop (is only run once)
    int anynul;
    if (fits_read_col(fptr, TFLOAT, 1, row, 1, 1, Emin, Emin, &anynul, &status)) break;
    if (fits_read_col(fptr, TFLOAT, 2, row, 1, 1, Emax, Emax, &anynul, &status)) break;
    if (fits_read_col(fptr, TINT, 3, row, 1, 1, Ngrp, Ngrp, &anynul, &status)) break;
    if (fits_read_col(fptr, TINT, 4, row, 1, *Ngrp, Fchan, Fchan, &anynul, &status)) break;
    if (fits_read_col(fptr, TINT, 5, row, 1, *Ngrp, Nchan, Nchan, &anynul, &status)) break;

    // determine number of columns in matrix for current row
    int ncols=0;
    for (count=0; count<*Ngrp; count++) {
      ncols += Nchan[count];
    }

    if (fits_read_col(fptr, TFLOAT, 6, row, 1, ncols, matrix, matrix, &anynul, &status)) break;
  } while(0);  // end of error handling loop

  return(status);
}




// read EBOUNDS data line from a FITS table (one PHA channel)
int read_ebounds_fitsrow(
			 long *channel, 
			 float *Emin,
			 float *Emax,
			 fitsfile *fptr, 
			 long row
			 )
{
  int status=EXIT_SUCCESS;  // error status

  do {  // beginning of error handling loop (is only run once)
    int anynul;
    if (fits_read_col(fptr, TLONG, 1, row, 1, 1, channel, channel, &anynul, &status)) break;
    if (fits_read_col(fptr, TFLOAT, 2, row, 1, 1, Emin, Emin, &anynul, &status)) break;
    if (fits_read_col(fptr, TFLOAT, 3, row, 1, 1, Emax, Emax, &anynul, &status)) break;
  } while(0);  // end of error handling loop

  return(status);
}


