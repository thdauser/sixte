#include "attitudefile.h"


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



// Writes a row of attitude data in to the FITS file.
int add_attitudetbl_row(fitsfile *fptr, 
			long row, 
			char valtime[], 
			double time, 
			double view_ra, 
			double view_dec, 
			double rollangle, 
			double aspangle, 
			int fitsstatus)
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



AttitudeFileEntry read_AttitudeFileEntry(AttitudeFile* af, int* status)
{
  AttitudeFileEntry afe = { .time = 0. };
  int anynul = 0;
 
  // TIME
  fits_read_col(af->fptr, TDOUBLE, af->ctime, af->row+1, 1, 1, 
		&afe.time, &afe.time, &anynul, status);
  // VIEWRA
  fits_read_col(af->fptr, TDOUBLE, af->cviewra, af->row+1, 1, 1, 
		&afe.viewra, &afe.viewra, &anynul, status);
  // VIEWDEC
  fits_read_col(af->fptr, TDOUBLE, af->cviewdec, af->row+1, 1, 1, 
		&afe.viewdec, &afe.viewdec, &anynul, status);
  // ROLLANG
  fits_read_col(af->fptr, TDOUBLE, af->crollang, af->row+1, 1, 1, 
		&afe.rollang, &afe.rollang, &anynul, status);

  return(afe);
}



AttitudeFile* open_AttitudeFile(const char filename[], int access_mode, int* status)
{
  AttitudeFile* af=NULL;
  char msg[MAXMSG];  // buffer for error messages

  do { // Beginning of ERROR handling loop

    af = (AttitudeFile*)malloc(sizeof(AttitudeFile));
    if (NULL==af) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation in attitude file "
	      "open routine failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    
    // Open the FITS file table for reading:
    af->fptr=NULL;
    if (fits_open_table(&af->fptr, filename, 
			access_mode, status)) break;

    // Get the HDU type
    int hdutype;
    if (fits_get_hdu_type(af->fptr, &hdutype, status)) break;

    // Image HDU results in an error message.
    if (hdutype==IMAGE_HDU) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: no table extension available in attitude "
	      "FITS file '%s'!\n", filename);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Determine the number of rows in the attitude file.
    if (fits_get_num_rows(af->fptr, &af->nrows, status)) 
      break;

    // Set internal row counter to first row (starting at 0).
    af->row = 0;


    // Determine the individual column numbers:
    // REQUIRED columns:
    if(fits_get_colnum(af->fptr, CASEINSEN, "TIME", &af->ctime, status)) break;
    if(fits_get_colnum(af->fptr, CASEINSEN, "VIEWRA", &af->cviewra, status)) break;
    if(fits_get_colnum(af->fptr, CASEINSEN, "VIEWDECL", &af->cviewdec, status)) break;
    if(fits_get_colnum(af->fptr, CASEINSEN, "ROLLANG", &af->crollang, status)) break;
  
  } while(0); // END of ERROR handling loop

  if(EXIT_SUCCESS!=*status) af=NULL;
  return(af);
}





