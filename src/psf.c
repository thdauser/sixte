/** 
 * This file contains all source code for PSF calculations.
 */
#include "psf.h"




// This function determines from the sky position of the source and the photon
// energy, which PSF data should be used to calculate the photon-detector hitting 
// point. It returns the corresponding PSF data structure.
// IMPORTANT: The function assumes that the individual PSF data sets lie on a 
// regular pattern: energy_{i,j} = energy_{i,k} and angle_{i,j} = angle_{k,j} !
static inline PSF_Item *get_best_psf_item(
					  double offaxis_angle, 
					  double energy,        
					  PSF *psf  // PSF (all angles & energies)
					  )
{
  // In order to find the PSF that matches the required position best, perform a loop
  // over all available PSFs and remember the best one.
  int count, index=0;
  double best_angle = psf->item[0].angle;
  double best_energy = psf->item[0].energy;
  for (count=1; count<psf->N_elements; count++) {
    assert(offaxis_angle >= 0.); // TODO

    // Check whether the current PSF parameters (off-axis angle and energy)
    // are better than the best values found so far.
    if ((fabs(psf->item[count].angle-offaxis_angle) - fabs(best_angle-offaxis_angle) 
	 < -0.000001) ||
	(fabs(psf->item[count].energy-energy) - fabs(best_energy-energy) 
	 < -0.000001)) {
      // Set is better => set new optimum values.
      index = count;
      best_angle  = psf->item[index].angle;
      best_energy = psf->item[index].energy;
    }
  }

  return(&psf->item[index]);
}






// Calculates the position on the detector, where a photon at given sky 
// position with specified energy hits the detector according to the PSF 
// data and a random number generator (randomization over one PSF pixel).
// Return value is '1', if the photon hits the detector. If it does not 
// fall onto the detector, the function returns '0'.
// The output detector position is stored in [mu m] in the first 2 parameters 
//of the function.
int get_psf_pos(
		// output: coordinates of the photon on the detector [mu m]
		struct Point2d* position,
		struct Photon photon,       // incident photon
		// telescope information (focal length, pointing direction)
		struct Telescope telescope, 
		PSF* psf
		)
{
  // Calculate the off-axis angle
  double offaxis_angle = acos(scalar_product(telescope.nz, photon.direction));
  // and the azimuth of the source position.
  double azimuth = atan2(scalar_product(telescope.ny,photon.direction), 
			 scalar_product(telescope.nx,photon.direction));

  // Determine, which PSF should be used for that particular source 
  // direction and photon energy.
  PSF_Item* psf_item = get_best_psf_item(offaxis_angle, photon.energy, psf);



  // Get a position from this closest PSF image using randomization.

  // Detector coordinates [pixel] of position obtained from closest PSF image.
  // All coordinates given in [pixel] are related to a origin in the corner 
  // of the detector.
  int x1, y1;   

  // get a random number to determine a random hitting position
  double rnd = get_random_number();
  if (rnd > psf_item->data[psf->width-1][psf->width-1]) {
    // The photon does not hit the detector at all (e.g. it is absorbed).
    n_outside++;
    return(0);
  }
  // Otherwise the photon hits the detector.
  // Perform a binary search to determine the position:
  // -> one binary search for each of the 2 coordinates x and y
  int high = psf->width-1;
  int low = 0;
  while (high-low > 1) {
    if (psf_item->data[(low+high)/2][0] < rnd) {
      low = (low+high)/2;
    } else {
      high = (low+high)/2;
    }
  }
  if (psf_item->data[low][psf->width-1] > rnd) {
    x1 = low;
  } else {
    x1 = high;
  }
    
  // Search for the y coordinate:
  high = psf->width-1;
  low = 0;
  while (high-low > 1) {
    if (psf_item->data[x1][(low+high)/2] < rnd) {
      low = (low+high)/2;
    } else {
      high = (low+high)/2;
    }
  }
  if (psf_item->data[x1][low] < rnd) {
    y1 = high;
  } else {
    y1 = low;
  }
  // Now x1 and y1 have pixel positions [integer pixel].
 

  // Randomize the [pixel] position (x1,y1), add the shift resulting 
  // from the off-axis angle difference (between actuall angle and the 
  // available angle in the PSF_Store), and transform all coordinates to [mu m]:
  double x2 = ((double)(x1-psf->width/2) + get_random_number()) *psf->pixelwidth -
              tan(offaxis_angle-psf_item->angle)*telescope.focal_length;
  double y2 = ((double)(y1-psf->width/2) + get_random_number()) *psf->pixelwidth;


  // Rotate the PSF postition [mu m] according to the azimuth angle.
  position->x = cos(azimuth)*x2 - sin(azimuth)*y2;
  position->y = sin(azimuth)*x2 + cos(azimuth)*y2;

  return(1);  
}





///////////////////////////////////////////////////////////////////////////
// Releases the memory which has been allocated to store the PSF data.
void free_psf(
	      PSF *psf  // pointer to the PSF data structure
	      )
{
  int count1, count2;

  if (psf->item) {
    for (count1=0; count1<psf->N_elements; count1++) {
      if (psf->item[count1].data) {
	for (count2=0; count2<psf->width; count2++) {
	  if (psf->item[count1].data[count2]) {
	    free(psf->item[count1].data[count2]);
	  }
	}
	free(psf->item[count1].data);
      }
    }
    free(psf->item);
    psf->item=NULL;
  }
}



/*
////////////////////////////////////////////////////////////
// This routine stores a PSF to a FITS file.
int save_psf_to_fits(
		     PSF *psf,
		     const char filename[],
		     int *status
		     )
{
  int count1, count2, count3;
  // 1D array used as buffer to store the 2D psf arrays (or a part of them) 
  // to the FITS file
  double *sub_psf=NULL;       
  fitsfile *output_fptr=NULL;

  char msg[MAXMSG];           // buffer for error output messages


  do {   // beginning of the error handling loop

    // Determine size of PSF sub-rectangles (don't save entire PSF but only 
    // the relevant region around the central peak, which has a probability 
    // greater than 0).
    int n = 0;     // width and
    int m = 0;     // height of sub-rectangle
    n = psf->width; // TODO
    m = psf->width;

    // delete old FITS output file
    remove(filename);

    // create new FITS output file
    headas_chat(5, "create FITS file '%s' ...\n", filename);
    if (fits_create_file(&output_fptr, filename, status)) break;

    // variables for the creation of the fits-table 
    // (field-type and -format of columns)
    char *ftype[PSF_NFIELDS];
    char *fform[PSF_NFIELDS];
    char *funit[PSF_NFIELDS];
    // create a binary table in the FITS file
    psf_create_tbl_parameter(ftype, fform, funit, psf->width);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, PSF_NFIELDS, ftype, fform, 
			funit, "PSF" , status)) break;
  
    // write headers
    fits_write_key(output_fptr, TSTRING, "TELESCOP", "eROSITA", 
		   "name of the telescope", status);
    fits_write_key(output_fptr, TSTRING, "COMMENT", "DESCRIPT", 
		   "simulated PSF for the eROSITA Wolter telescope",status);
    fits_write_key(output_fptr, TINT, "WIDTH", &psf->width, 
		   "width of the entire PSF [pixel]",status);
    fits_write_key(output_fptr, TDOUBLE, "PIXWIDTH", &psf->pixelwidth, 
		   "width of the PSF pixels in [mu m]",status);
    fits_write_key(output_fptr, TINT, "n", &n, "width of PSF sub-rectangle",status);
    fits_write_key(output_fptr, TINT, "m", &m, "height of PSF sub-rectangle",status);

    // If desired by the user, print all program parameters to HISTORY 
    // of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, status);

    // check, if errors have occurred on writing the headers
    if (*status) break;
    headas_chat(5, "headers written into FITS file '%s' ...\n", filename);


    // create relevant PSF sub-rectangle
    sub_psf = (double *) malloc((long)n*(long)m*sizeof(double));
    if (!sub_psf) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error allocating memory!\n");
      HD_ERROR_THROW(msg,*status);
      break;
    }

    // Copy sub-rectangle to 1d-array for each energy and off-axis angle.
    long row=0;      // row in the FITS table
    for (count1=0; count1<psf->N_elements; count1++) {
      int x=0;       // coordinates of upper left corner of sub-rectangle
      int y=0;

      for (count2=x; count2<x+n; count2++) {
	for (count3=y; count3<y+m; count3++) {
	  sub_psf[((count2-x)*n+count3-y)] = psf->item[count1].data[count2][count3];
	}
      }
    
      // write data row to FITS file                                
      if ((*status=insert_psf_fitsrow(psf->item[count1].angle, 
				      psf->item[count1].energy, x, y, sub_psf, n*m, 
				      output_fptr, ++row))) break;      
    }
    
  } while (0);   // end of the error handling loop


  // ------- cleaning up -----------

  // release the memory of the buffer array
  if (!sub_psf) free(sub_psf);

  // close the FITS file
  if(output_fptr) fits_close_file(output_fptr, status);

  return(*status);
}
*/






////////////////////////////////////////////////////////////////////////
// This routine creates the necessary data for the FITS-table layout.
static void psf_create_tbl_parameter(
				     char *ftype[PSF_NFIELDS], 
				     char *fform[PSF_NFIELDS], 
				     char *funit[PSF_NFIELDS],
				     int width // width of the PSF array in [pixel]
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








////////////////////////////////////////////////
// Routine reads PSF data from a file with FITS images.
PSF* get_psf(
	     const char* filename,
	     int* status
	     )
{
  PSF* psf;
  fitsfile* fptr=NULL;   // FITSfile-pointer to PSF file
  double* data;          // input buffer (1D array)
  long count, count2, count3;
  
  char msg[MAXMSG];      // error message output buffer


  do {  // beginning of error handling loop

    // Allocate memory for PSF data structure:
    psf = (PSF*)malloc(sizeof(PSF));
    if (psf == NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: memory allocation for PSF structure failed!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Open PSF FITS file
    headas_chat(5, "open PSF FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;
    
    // Get the number of HDUs in the FITS file.
    if (fits_get_num_hdus(fptr, &psf->N_elements, status)) break;
    
    // Determine the width of the PSF image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    if (naxes[0] != naxes[1]) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: PSF image must be a square!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      psf->width = (int)naxes[0];
    }
    
    // Determine the pixelwidth of the PSF array from the header keywords.
    char comment[MAXMSG]; // buffer 
    if (fits_read_key(fptr, TDOUBLE, "PIXWIDTH", &psf->pixelwidth, comment, 
		      status)) break;

    
    // Get memory for the PSF.
    psf->item = (PSF_Item *) malloc(psf->N_elements * sizeof(PSF_Item));
    if (psf->item) {   // memory was allocated successfully
      for (count=0; count<psf->N_elements; count++) {
	psf->item[count].data = (double **) malloc(psf->width * sizeof(double *));
	if (psf->item[count].data) {
	  for (count2=0; count2<psf->width; count2++) {
	    psf->item[count].data[count2] = (double *) 
	      malloc(psf->width * sizeof(double));
	    if (!psf->item[count].data[count2]) {
	      *status = EXIT_FAILURE;
	    }
	  }
	} else { *status = EXIT_FAILURE; }
      }
    } else { *status = EXIT_FAILURE; }
    // Check if all memory was allocated successfully
    if (*status != EXIT_SUCCESS) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg, *status);  
      break;
    }

    // Allocate memory for input buffer (1D array)
    data=(double*)malloc(psf->width*psf->width*sizeof(double));
    if (!data) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory for PSF input buffer!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }



    // Loop over the individual PSFs for the storage.
    for (count=0; (count<psf->N_elements)&&(*status==EXIT_SUCCESS); count++) {
      
      int hdutype;
      if (fits_movabs_hdu(fptr, count+1, &hdutype, status)) break;

      // Read the PSF information (energy, off-axis angle) from the
      // header keywords in each FITS HDU.
      if (fits_read_key(fptr, TDOUBLE, "ENERGY", &psf->item[count].energy, comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "OFFAXANG", &psf->item[count].angle, comment, 
			status)) break;      
      // convert the off-axis angle from [degree] to [rad]
      psf->item[count].angle = psf->item[count].angle * M_PI/180.;


      int anynul;
      double null_value=0.;
      long fpixel[2] = {1, 1};   // lower left corner
      //                |--|--> FITS coordinates start at (1,1)
      long lpixel[2] = {psf->width, psf->width};  // upper right corner
      long inc[2] = {1, 1};

      if (fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, 
			   data, &anynul, status)) break;


      // Create a partition function from the 1D PSF data array,
      // i.e., sum up the individual probabilites.
      // The partition function is more adequate for determining a random 
      // photon impact position on the detector.
      double sum=0.;
      for (count2=0; count2<psf->width; count2++) {
	for (count3=0; count3<psf->width; count3++) {
	  sum += data[count2*psf->width+count3];
	  psf->item[count].data[count2][count3] = sum;
	}
      }

      // Store the integrated on-axis PSF for each energy band. (TODO)
      psf->item[count].scaling_factor = 1.; // sum;
      

      // Renormalize the PSF partition function to the integrated on-axis PSF.
      for (count2=0; count2<psf->width; count2++) {
	for (count3=0; count3<psf->width; count3++) {
	  psf->item[count].data[count2][count3] = 
	    psf->item[count].data[count2][count3] / psf->item[count].scaling_factor;
	}
      }

      // Plot normalization of PSF for current off-axis angle and energy
      headas_chat(5, "PSF: %lf of incident photons at (%lf rad, %lf keV), "
		  "normalized to %lf, factor 1/%lf\n",  sum, 
		  psf->item[count].angle, psf->item[count].energy, 
		  sum/psf->item[count].scaling_factor, 
		  psf->item[count].scaling_factor);

    } // END of loop over individual PSF items
  } while(0);  // END of error handling loop


  // Close PSF file.
  if (fptr) fits_close_file(fptr, status);

  // free memory of input buffer
  if (data) free(data);  

  return(psf);
}





/////////////////////////////////////////////
int save_psf_image(
		   PSF* psf,
		   const char *filename,
		   int *status
		   )
{
  int count;
  double *sub_psf=NULL;
  fitsfile *fptr;

  char msg[MAXMSG];           // buffer for error output messages


  do { // ERROR handling loop

    // Create a new FITS-file:
    if (fits_create_file(&fptr, filename, status)) break;


    // Loop over the different PSFs in the storage:
    for (count=0; count<psf->N_elements; count++) {

      // Determine size of PSF sub-rectangles (don't save entire PSF but only 
      // the relevant region around the central peak, which has a probability 
      // greater than 0).
      int n = 0;     // width and
      int m = 0;     // height of sub-rectangle
      n = psf->width; // TODO
      m = psf->width;

      // Create the relevant PSF sub-rectangle:
      sub_psf = (double *) malloc((long)n*(long)m*sizeof(double));
      if (!sub_psf) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error allocating memory!\n");
	HD_ERROR_THROW(msg, *status);
	break;
      }
      // Store the PSF in the 1D array to handle it to the FITS routine.
      int x0=0, y0=0; // coordinates of lower left corner of sub-rectangle
      int x, y;
      for (x=x0; x<x0+n; x++) {
	for (y=y0; y<y0+m; y++) {
	  sub_psf[((x-x0)*n+y-y0)] = psf->item[count].data[x][y];
	}
      }
    

      // Create an image in the FITS-file (primary HDU):
      long naxes[2] = {(long)(psf->width), (long)(psf->width)};
      if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
      //                                    |-> naxis
      int hdutype;
      if (fits_movabs_hdu(fptr, count+1, &hdutype, status)) break;



      // Write the header keywords for PSF FITS-files (CAL/GEN/92-027):
      double dbuffer;
      fits_write_key(fptr, TSTRING, "CTYPE1", "DETX",
		     "detector coordinate system", status);
      fits_write_key(fptr, TSTRING, "CTYPE2", "DETY",
		     "detector coordinate system", status);

      fits_write_key(fptr, TSTRING, "HDUCLASS", "OGIP",
		     "Extension is OGIP defined", status);
      fits_write_key(fptr, TSTRING, "HDUDOC", "CAL/GEN/92-020",
		     "Document containing extension definition", status);
      fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0",
		     "giving the version of the format", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS1", "Image",
		     "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS2", "PSF",
		     "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS3", "PREDICTED",
		     "", status);
      fits_write_key(fptr, TSTRING, "HDUCLAS4", "NET",
		     "", status);

      fits_write_key(fptr, TDOUBLE, "TSTRING", "pixel", 
		     "", status);
      fits_write_key(fptr, TDOUBLE, "TSTRING", "pixel", 
		     "", status);

      dbuffer = 0.;
      fits_write_key(fptr, TDOUBLE, "CRPIX1", &dbuffer, 
		     "", status);
      fits_write_key(fptr, TDOUBLE, "CRPIX2", &dbuffer, 
		     "", status);
      dbuffer = 0.;
      fits_write_key(fptr, TDOUBLE, "CRVAL1", &dbuffer, 
		     "", status);
      fits_write_key(fptr, TDOUBLE, "CRVAL2", &dbuffer, 
		     "", status);
      dbuffer = 0.;
      fits_write_key(fptr, TDOUBLE, "CDELT1", &dbuffer, 
		     "", status);
      fits_write_key(fptr, TDOUBLE, "CDELT2", &dbuffer, 
		     "", status);


      dbuffer = 0.0;
      fits_write_key(fptr, TDOUBLE, "BACKGRND", &dbuffer, 
		     "background count rate per pixel", status);


      // Mission
      fits_write_key(fptr, TSTRING, "TELESCOP", "", // eROSITA
		     "Satellite", status);
      fits_write_key(fptr, TSTRING, "INSTRUME", "", // pnCCD1
		     "Instrument", status);
      fits_write_key(fptr, TSTRING, "FILTER", "NONE",
		     "Filter", status);
    
      // TODO: Instead of using these keywords one should specifiy the 
      // physical units of the images.
      fits_write_key(fptr, TDOUBLE, "PIXWIDTH", &psf->pixelwidth, 
		     "width of the PSF pixels in [mu m]", status);
      fits_write_key(fptr, TDOUBLE, "ENERGY", &psf->item[count].energy, 
		     "photon energy for the PSF generation in [keV]", status);
      fits_write_key(fptr, TDOUBLE, "OFFAXANG", &psf->item[count].angle, 
		     "off-axis angle in [deg]", status);

      
      HDpar_stamp(fptr, count+1, status);
      
      // END of writing header keywords.


      // Write the image to the file:
      long fpixel[2] = {x0+1, y0+1};                // lower left corner
      //                   |-----|--> FITS coordinates start at (1,1)
      long lpixel[2] = {psf->width, psf->width};    // upper right corner
      fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, sub_psf, status);
      
    } // END of loop over individual PSF items in the storage.

  } while (0); // END of ERROR handling loop


  // close the FITS-file
  if(fptr) fits_close_file(fptr, status);

  return(*status);
}





/*
// read PSF data from a FITS table
int read_psf_fitsrow(
		     struct PSF *psf,
		     double *data,
		     long size,
		     fitsfile *fptr, 
		     long row
		     )
{
  int status=EXIT_SUCCESS;  // error status

  psf->angle=0.;
  psf->energy=0.;
  long count;
  for (count=0; count<size; count++) { data[count] = 0.; }

  do {  // beginning of error handling loop (is only run once)
    int anynul;
    if (fits_read_col(fptr, TDOUBLE, 1, row, 1, 1, &psf->angle, &psf->angle, 
		      &anynul, &status)) break;
    if (fits_read_col(fptr, TDOUBLE, 2, row, 1, 1, &psf->energy, &psf->energy, 
		      &anynul, &status)) break;
    //    if (fits_read_col(fptr, TINT, 3, row, 1, 1, x, x, &anynul, &status)) break;
    //    if (fits_read_col(fptr, TINT, 4, row, 1, 1, y, y, &anynul, &status)) break;
    if (fits_read_col(fptr, TDOUBLE, 5, row, 1, size, data, data, &anynul, &status)) 
      break;
  } while(0);  // end of error handling loop

  // convert the off-axis angle from [degree] to [rad]
  psf->angle = psf->angle * M_PI/180.;

  return(status);
}





////////////////////////////////////////////////
// Routine reads PSF data from a FITS image.
int read_psf_image(
		   struct PSF *psf,
		   double *data,
		   long width,
		   fitsfile *fptr
		   )
{
  int status=EXIT_SUCCESS;  // error status

  psf->angle=0.;
  psf->energy=0.;
  long count;
  for (count=0; count<width*width; count++) { data[count] = 0.; }

  do {  // beginning of error handling loop (is only run once)
    int anynul;
    double null_value=0.;
    long fpixel[2] = {1, 1};          // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    long lpixel[2] = {width, width};  // upper right corner
    long inc[2] = {1, 1};

    if (fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, data, 
			 &anynul, &status)) break;
  } while(0);  // END of error handling loop

  // convert the off-axis angle from [degree] to [rad]
  psf->angle = psf->angle * M_PI/180.;

  return(status);
}

*/



// Reads the PSF data file (ASCII eventlist), bins the events to detector pixels and
// stores this data in an array. 
// Additionally the function creates an array containing important information about 
// the PSF, i.e. the array indices which correspond to the desired  energy or 
// off-axis angle.

/*
int get_psf(
	    // Structure to store all PSF data 
	    // (for several discrete angles and energies)
	    struct PSF_Store *store,          
	    // FITS file containing images of the PSF
	    const char filename[FILENAME_LENGTH] 
	    )
{
  int count1, count2, count3;
  fitsfile *fptr=NULL;   // FITSfile-pointer to PSF file
  double *data=NULL;     // input buffer for PSF

  int status=0;          // error status
  char msg[MAXMSG];      // error message output buffer


  do {  // error handling loop (only run once)

    // first open PSF FITS file
    headas_chat(5, "open PSF FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, &status)) break;
  

    // Get memory for the PSF.
    store->psf = (struct PSF *) malloc(store->N_elements * sizeof(struct PSF));
    if (store->psf) {   // memory was allocated successfully
      for (count1=0; count1<store->N_elements; count1++) {
	store->psf[count1].data = (double **) malloc(store->width * sizeof(double *));
	if (store->psf[count1].data) {
	  for (count2=0; count2<store->width; count2++) {
	    store->psf[count1].data[count2] = (double *) 
	      malloc(store->width * sizeof(double));
	    if (!store->psf[count1].data[count2]) {
	      status = EXIT_FAILURE;
	    }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    
    // check if all necessary memory was allocated successfully
    if (status == EXIT_FAILURE) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg,status);  
      break;
    }

    // get memory for input buffer (1D array)
    data=(double*)malloc(store->width*store->width*sizeof(double));
    if (!data) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory for input buffer!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
 

    /////////////////////////////////////////////////////////
    // fill the psf.data array with data from the FITS file
    headas_chat(5, "load PSF data ...\n");
  
    // loop over all different off-axis angles (break, if an error has occured)
    for (count1=0; (count1<store->N_elements)&&(status==EXIT_SUCCESS); count1++) {
	
      // Create a partition function from the 1D PSF data array,
      // i.e., sum up the individual probabilites.
      // The partition function is more adequate for determining a random 
      // photon impact position on the detector.
      double sum=0.;
      for (count2=0; count2<store->width; count2++) {
	for (count3=0; count3<store->width; count3++) {
	  sum += data[count2*store->width+count3];
	  store->psf[count1].data[count2][count3] = sum;
	}
      }

      // Store the integrated on-axis PSF for each energy band. (TODO)
      store->psf[count1].scaling_factor = 1.; // sum;
      
      // Renormalize the PSF partition function to the integrated on-axis PSF.
      for (count2=0; count2<store->width; count2++) {
	for (count3=0; count3<store->width; count3++) {
	  store->psf[count1].data[count2][count3] = 
	    store->psf[count1].data[count2][count3]/store->psf[count1].scaling_factor;
	}
      }

      // Plot normalization of PSF for current off-axis angle and energy
      headas_chat(5, "PSF: %lf of incident photons at (%lf rad, %lf keV), "
		  "normalized to %lf, factor 1/%lf\n",  sum, 
		  store->psf[count1].angle, store->psf[count1].energy, 
		  sum/store->psf[count1].scaling_factor, 
		  store->psf[count1].scaling_factor);

    }  // END of loop over all rows in the FITS file
    
  } while (0);   // END of error handling loop


  // close psf file
  if (fptr) fits_close_file(fptr, &status);

  // free memory of input buffer
  if (data) free(data);  

  return(status);
}







// Reads the PSF data file (ASCII eventlist), bins the events to detector pixels and
// stores this data in an array. 
// Additionally the function creates an array containing important information about 
// the PSF, i.e. the array indices which correspond to the desired  energy or 
// off-axis angle.
int get_psf_old(
	    // structure containing all PSF data (for all angles and energies)
	    struct PSF_Store *store,          
	    const char filename[FILENAME_LENGTH] // FITS file containing the PSF
	    )
{
  int count1, count2, count3;
  fitsfile *fptr=NULL;   // FITSfile-pointer to PSF file
  double *data=NULL;     // input buffer for PSF

  int status=0;          // error status
  char msg[MAXMSG];      // error message output buffer


  do {  // error handling loop (only run once)

    // first open PSF FITS file
    headas_chat(5, "open PSF FITS file '%s' ...\n", filename);
    if (fits_open_file(&fptr, filename, READONLY, &status)) break;
  
    int hdunum, hdutype;
    // after opening the FITS file, get the number of the current HDU
    if (fits_get_hdu_num(fptr, &hdunum) == 1) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table
      if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(fptr, &hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in file '%s' is not a table but "
	      "an image (HDU number: %d)\n", filename, hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Get the number of rows in the FITS file (number of given points of time).
    long nrows;
    fits_get_num_rows(fptr, &nrows, &status);

    // The total number of different PSFs is equal to the number of rows in 
    // the FITS file.
    store->N_elements = (int)nrows;

    // Determine the width and the pixelwidth of the PSF arrays.
    char comment[MAXMSG];
    if (fits_read_key(fptr, TINT, "WIDTH", &store->width, comment, &status)) break;
    if (fits_read_key(fptr, TDOUBLE, "PIXWIDTH", &store->pixelwidth, comment, 
		      &status)) break;
    

    // Get memory for the PSF.
    store->psf = (struct PSF *) malloc(store->N_elements * sizeof(struct PSF));
    if (store->psf) {   // memory was allocated successfully
      for (count1=0; count1<store->N_elements; count1++) {
	store->psf[count1].data = (double **) malloc(store->width * sizeof(double *));
	if (store->psf[count1].data) {
	  for (count2=0; count2<store->width; count2++) {
	    store->psf[count1].data[count2] = (double *) 
	      malloc(store->width * sizeof(double));
	    if (!store->psf[count1].data[count2]) {
	      status = EXIT_FAILURE;
	    }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    
    // check if all necessary memory was allocated successfully
    if (status == EXIT_FAILURE) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg,status);  
      break;
    }

    // get memory for input buffer
    data=(double*)malloc(store->width*store->width*sizeof(double));
    if (!data) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory for input buffer!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
 

    /////////////////////////////////////////////////////////
    // fill the psf.data array with data from the FITS file
    headas_chat(5, "load PSF data ...\n");
  
    // loop over all different off-axis angles (break, if an error has occured)
    for (count1=0; (count1<store->N_elements)&&(status==EXIT_SUCCESS); count1++) {
      // input data row from FITS file
      //int x=0, y=0;
      if ((status=read_psf_fitsrow(&(store->psf[count1]), data, 
				   store->width * store->width, fptr, 
				   count1+1))!=EXIT_SUCCESS) break;
	
	
      // Create a partition function from the 1D PSF data array,
      // i.e., sum up the individual probabilites.
      // The partition function is more adequate for getting a random 
      // photon impact position on the detector.
      double sum=0.;
      for (count2=0; count2<store->width; count2++) {
	for (count3=0; count3<store->width; count3++) {
	  sum += data[count2*store->width+count3];
	  store->psf[count1].data[count2][count3] = sum;
	}
      }

      // Store the integrated on-axis PSF for each energy band. (TODO)
      store->psf[count1].scaling_factor = 1.; // sum;
      
      // Renormalize the PSF partition function to the integrated on-axis PSF.
      for (count2=0; count2<store->width; count2++) {
	for (count3=0; count3<store->width; count3++) {
	  store->psf[count1].data[count2][count3] = 
	    store->psf[count1].data[count2][count3]/store->psf[count1].scaling_factor;
	}
      }

      // Plot normalization of PSF for current off-axis angle and energy
      headas_chat(5, "PSF: %lf of incident photons at (%lf rad, %lf keV), "
		  "normalized to %lf, factor 1/%lf\n",  sum, 
		  store->psf[count1].angle, store->psf[count1].energy, 
		  sum/store->psf[count1].scaling_factor, 
		  store->psf[count1].scaling_factor);

    }  // END of loop over all rows in the FITS file
    
  } while (0);   // end of error handling loop


  // close psf file
  if (fptr) fits_close_file(fptr, &status);

  // free memory of input buffer
  if (data) free(data);  

  return(status);
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
*/


