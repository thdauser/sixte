/** 
 * This file belongs to 'measurement.c' and
 * contains all source code for orbit and attitude handling.
 */

#include "sixt.h"
#include "orbatt.h"



/** Constructor for the AttitudeCatalog. */
AttitudeCatalog* get_AttitudeCatalog(
				     const char attitude_filename[],
				     double t0, double timespan,
				     int* status
				     )
{
  AttitudeCatalog* ac=NULL;

  fitsfile *att_fptr=NULL; // FITS file pointer to attitude file
  int hdunum, hdutype;     // HDU number and type for FITS access
  long nrows;              // number of lines in the FITS file
  long counter, entry=0;   // Counters for the AttitudeEntry elements in the FITS
                           // file and int the AttitudeCatalog respectively. 
  char msg[MAXMSG];

  do {  // beginning of ERROR handling loop

    // Read-in the attitude data from the FITS file and store 
    // them in the AttitudeCatalog.
    // Open the attitude file:
    if (fits_open_table(&att_fptr, attitude_filename, READONLY, status)) break;

    // After opening the attitude file, get the number of the current HDU
    if (fits_get_hdu_num(att_fptr, &hdunum) == 1) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table
      if (fits_movabs_hdu(att_fptr, 2, &hdutype, status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(att_fptr, &hdutype, status)) break;
    }

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in attitude file '%s' is not a table "
	      "but an image (HDU number: %d)\n", attitude_filename, hdunum);
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Get the number of rows in the attitude file 
    // (number of given points of time)
    fits_get_num_rows(att_fptr, &nrows, status);


    // Get memory for the AttitudeCatalog:
    ac = (AttitudeCatalog*)malloc(sizeof(AttitudeCatalog));
    if (NULL==ac) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the AttitudeCatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    ac->entry = (AttitudeEntry*)malloc(nrows*sizeof(AttitudeEntry));
    if (NULL==ac->entry) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the AttitudeCatalog!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    
    /* Point of time at which the attitude data are valid. */
    double time;
    char valtime[19];
    /* Angles specifying the direction of the telescope (attitude data). */
    double view_ra, view_dec, rollangle;  

    // Read all lines from attitude file subsequently.
    for (counter=0, entry=0; (counter<nrows)&&(!*status); counter++) {
      if (get_atttbl_row(att_fptr, counter, valtime, &time, &view_ra, 
			 &view_dec, &rollangle, status)) break;

      if (time >= t0) {
	// Check if we are already at the end of the attitude file:
	if (counter>=nrows) {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: Not enough attitude data available for the "
		  "specified period from %lf to %lf!", t0, timespan);
	  HD_ERROR_THROW(msg, *status);
	  break;
	}
	  
	// Calculate and store attitude data:
	ac->entry[entry].time = time;
	// Rescale from degrees to radians:
	rollangle = rollangle*M_PI/180.;

	// Telescope pointing direction nz:
	ac->entry[entry].nz = unit_vector(view_ra*M_PI/180., view_dec*M_PI/180.);

	entry++;
      }
	  
      // Check whether the end of the reading process 
      // (end of specified period) is reached.
      if (time > t0+timespan) {
	break;
      }
    }  // End of the attitude readout loop
    
    // Save the number of AttitudeEntry elements.
    ac->nentries = entry;


    // Loop over all AttitudeEntry elements in the AttitudeCatalog in
    // order to determine the telescope nx-direction (so far only the
    // nz direction is set).
    // Check the change of the telescope pointing direction between two subsequent
    // AttitudeEntry elements.
    Vector dnz = 
      vector_difference(ac->entry[1].nz, ac->entry[0].nz);
    if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { 
      // Change of the telescope axis is too small to be significant.
      Vector ny = {0., 1., 0.};
      ac->entry[0].nx=vector_product(ac->entry[0].nz, ny); // TODO
    } else {
      // nx = (nz_0 x nz_1) x nz_0
      ac->entry[0].nx=
	normalize_vector(vector_product(vector_product(ac->entry[0].nz,
						       ac->entry[1].nz),
					ac->entry[0].nz));
    }

    for (entry=1; entry<ac->nentries; entry++) {

      Vector dnz = 
	vector_difference(ac->entry[entry].nz, ac->entry[entry-1].nz);
      if (sqrt(scalar_product(&dnz, &dnz))<1.e-7) { 
	// Change of the telescope axis is too small to be significant.
	Vector ny = {0., 1., 0.};
	ac->entry[entry].nx=vector_product(ac->entry[entry].nz, ny); // TODO
      } else {
	ac->entry[entry].nx=
	  normalize_vector(vector_difference(ac->entry[entry].nz,
					     ac->entry[entry-1].nz));
      }

    } // END of loop over all AttitudeEntry elements for the calculation of nx.

  } while (0); // End of error handling loop


  // --- clean up ---

  // close FITS files
  if (att_fptr) fits_close_file(att_fptr, status);

  if (EXIT_SUCCESS != *status) ac = NULL;
  return(ac);
}




/** Destructor for the AttitudeCatalog */
void free_AttitudeCatalog(AttitudeCatalog* ac)
{
  if (NULL != ac) {
    if (NULL != ac->entry) {
      free(ac->entry);
    }
    free(ac);
  }
}


