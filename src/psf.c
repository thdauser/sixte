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
					  double offaxis_angle, // [rad]
					  double energy,        // [rad]
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






/** Calculates the position on the detector, where a photon at given sky 
 * position with specified energy hits the detector according to the PSF 
 * data and a random number generator (randomization over one PSF pixel).
 * Return value is '1', if the photon hits the detector. If it does not 
 * fall onto the detector, the function returns '0'.
 * The output detector position is stored in [m] in the first 2 parameters 
 * of the function. */
int get_psf_pos(
		// output: coordinates of the photon on the detector [m]
		struct Point2d* position,
		Photon photon,     // incident photon
		// telescope information (focal length, pointing direction)
		struct Telescope telescope, 
		Vignetting* vignetting,
		PSF* psf
		)
{
  // Calculate the off-axis angle
  double offaxis_angle = acos(scalar_product(&telescope.nz, &photon.direction));
  // and the azimuth of the source position.
  double azimuth = atan2(scalar_product(&telescope.ny, &photon.direction), 
			 scalar_product(&telescope.nx, &photon.direction));

  // Determine, which PSF should be used for that particular source 
  // direction and photon energy.
  PSF_Item* psf_item = get_best_psf_item(offaxis_angle, photon.energy, psf);


  // Get a position from this closest PSF image using randomization.

  // Detector coordinates [pixel] of position obtained from closest PSF image.
  // All coordinates given in [pixel] are related to a origin in the corner 
  // of the detector.
  int x1, y1;   

  // Get a random number to determine a random hitting position
  double rnd = get_random_number();
  //  if (rnd > psf_item->data[psf_item->naxis1-1][psf_item->naxis2-1]) {
  if (rnd > get_Vignetting_Factor(vignetting, photon.energy, offaxis_angle, azimuth)) {
    // The photon does not hit the detector at all (e.g. it is absorbed).
    return(0);
  }
  // Otherwise the photon hits the detector.
  // Perform a binary search to determine a random position:
  // -> one binary search for each of the 2 coordinates x and y
  rnd = get_random_number();
  int high = psf_item->naxis1-1;
  int low = 0;
  while (high-low > 1) {
    if (psf_item->data[(low+high)/2][0] < rnd) {
      low = (low+high)/2;
    } else {
      high = (low+high)/2;
    }
  }
  if (psf_item->data[low][psf_item->naxis1-1] > rnd) {
    x1 = low;
  } else {
    x1 = high;
  }
    
  // Search for the y coordinate:
  high = psf_item->naxis2-1;
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
  // available angle in the PSF_Store), and transform all coordinates to [m]:
  //double x2 = ((double)(x1-psf->width/2) + get_random_number()) *psf->pixelwidth -
  //              tan(offaxis_angle-psf_item->angle)*telescope.focal_length;
  //double y2 = ((double)(y1-psf->width/2) + get_random_number()) *psf->pixelwidth;
  double x2 = ((double)x1-psf_item->crpix1+1. + get_random_number())*psf_item->cdelt1 
    + psf_item->crval1 - tan(offaxis_angle-psf_item->angle)*telescope.focal_length;
  double y2 = ((double)y1-psf_item->crpix2+1.)*psf_item->cdelt2 + psf_item->crval2;

  // Rotate the PSF postition [mu m] according to the azimuth angle.
  position->x = cos(azimuth)*x2 - sin(azimuth)*y2;
  position->y = sin(azimuth)*x2 + cos(azimuth)*y2;

  return(1);  
}





////////////////////////////////////////////////////////////////////////
// Releases the memory which has been allocated to store the PSF data.
void free_psf(
	      PSF *psf  // pointer to the PSF data structure
	      )
{
  if (NULL!=psf) {
    int count1, count2;
    if (psf->item) {
      for (count1=0; count1<psf->N_elements; count1++) {
	if (psf->item[count1].data) {
	  for (count2=0; count2<psf->item[count1].naxis1; count2++) {
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
    free(psf);
  }
}





///////////////////////////////////////////////////////
/** Constructor for the PSF data structure.
 * Reads PSF data from a FITS file with image extensions. 
 * The file format is given by OGIP Calibration Memo 
 * CAL/GEN/92-027. */
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


  do {  // beginning of ERROR handling loop

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

    // Allocate memory for all PSF image extensions / PSF_Items in the
    // PSF data structure.
    psf->item = (PSF_Item *) malloc(psf->N_elements * sizeof(PSF_Item));
    if (NULL==psf->item) {   // memory was allocated successfully
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg, *status);  
      break;
    }


    // Loop over the individual PSF image extensions in the FITS file.
    for (count=0; (count<psf->N_elements)&&(*status==EXIT_SUCCESS); count++) {

      // Move to the right FITS extension.
      int hdutype;
      if (fits_movabs_hdu(fptr, count+1, &hdutype, status)) break;

      // Determine the width of the PSF image.
      long naxes[2];
      if (fits_get_img_size(fptr, 2, naxes, status)) break;
      psf->item[count].naxis1 = (int)naxes[0];
      psf->item[count].naxis2 = (int)naxes[1];

    
      // Determine the WCS keywords of the PSF array from the header.
      char comment[MAXMSG];
      if (fits_read_key(fptr, TDOUBLE, "CDELT1", &(psf->item[count].cdelt1), comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "CDELT2", &(psf->item[count].cdelt2), comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &(psf->item[count].crpix1), comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &(psf->item[count].crpix2), comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &(psf->item[count].crval1), comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &(psf->item[count].crval2), comment, 
			status)) break;

    
      // Get memory for the PSF_Item data.
      psf->item[count].data = (double **)
	malloc(psf->item[count].naxis1*sizeof(double *));
      if (NULL!=psf->item[count].data) {
	for (count2=0; count2<psf->item[count].naxis1; count2++) {
	  psf->item[count].data[count2] = (double *)
	    malloc(psf->item[count].naxis2 * sizeof(double));
	  if (NULL==psf->item[count].data[count2]) {
	    *status = EXIT_FAILURE;
	  }
	}
      } else { *status = EXIT_FAILURE; }
      // Check if all memory was allocated successfully
      if (*status != EXIT_SUCCESS) {
	sprintf(msg, "Error: not enough memory to store PSF data!\n");
	HD_ERROR_THROW(msg, *status);  
	break;
      }

      // Allocate memory for input buffer (1D array)
      data=(double*)malloc(psf->item[count].naxis1 * psf->item[count].naxis2
			   *sizeof(double));
      if (NULL==data) {
	*status = EXIT_FAILURE;
	sprintf(msg, "Error: not enough memory for PSF input buffer!\n");
	HD_ERROR_THROW(msg, *status);
	break;
      }
      
      // Read the PSF information (energy, off-axis angle) from the
      // header keywords in each FITS HDU.
      if (fits_read_key(fptr, TDOUBLE, "ENERGY", &psf->item[count].energy, comment, 
			status)) break;
      if (fits_read_key(fptr, TDOUBLE, "OFFAXANG", &psf->item[count].angle, comment, 
			status)) break;      
      // Convert the off-axis angle from [degree] to [rad]:
      psf->item[count].angle = psf->item[count].angle * M_PI/180.;


      int anynul;
      double null_value=0.;
      long fpixel[2] = {1, 1};   // lower left corner
      //                |--|--> FITS coordinates start at (1,1)
      // upper right corner
      long lpixel[2] = {psf->item[count].naxis1, psf->item[count].naxis2};  
      long inc[2] = {1, 1};

      if (fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, 
			   data, &anynul, status)) break;


      // Create a partition function from the 1D PSF data array,
      // i.e., sum up the individual probabilites.
      // The partition function is more adequate for determining 
      // a random photon impact position on the detector.
      double sum=0.;
      for (count2=0; count2<psf->item[count].naxis1; count2++) {
	for (count3=0; count3<psf->item[count].naxis2; count3++) {
	  // Take care of choosing x- and y-axis properly!
	  sum += data[count3*psf->item[count].naxis1+count2];  
	  psf->item[count].data[count2][count3] = sum;
	}
      }

      // Plot normalization of PSF for current off-axis angle and energy
      headas_chat(5, "PSF: images %.2lf%% of incident photons for "
		  "%.4lf deg, %.1lf keV\n", sum * 100., 
		  psf->item[count].angle*180./M_PI, 
		  psf->item[count].energy);

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

  char msg[MAXMSG];  // buffer for error output messages


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
      n = psf->item[count].naxis1; // TODO
      m = psf->item[count].naxis2;

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
      long naxes[2] = {(long)(psf->item[count].naxis1), (long)(psf->item[count].naxis2)};
      if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
      //                                    |-> naxis
      int hdutype;
      if (fits_movabs_hdu(fptr, count+1, &hdutype, status)) break;



      // Write the header keywords for PSF FITS-files (CAL/GEN/92-027):
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

      fits_write_key(fptr, TSTRING, "CUNIT1", "pixel", "", status);
      fits_write_key(fptr, TSTRING, "CUNIT2", "pixel", "", status);
      double dbuffer = psf->item[count].naxis1*0.5;
      fits_write_key(fptr, TDOUBLE, "CRPIX1", &dbuffer, 
		     "X axis reference pixel", status);
      fits_write_key(fptr, TDOUBLE, "CRPIX2", &dbuffer, 
		     "Y axis reference pixel", status);
      dbuffer = 0.;
      fits_write_key(fptr, TDOUBLE, "CRVAL1", &dbuffer, 
		     "coord of X ref pixel", status);
      fits_write_key(fptr, TDOUBLE, "CRVAL2", &dbuffer, 
		     "coord of Y ref pixel", status);
      fits_write_key(fptr, TDOUBLE, "CDELT1", &psf->item[count].cdelt1, 
		     "X axis increment", status);
      fits_write_key(fptr, TDOUBLE, "CDELT2", &psf->item[count].cdelt2, 
		     "Y axis increment", status);

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

      // Write the ENERGY and OFF-AXIS ANGLE for this particular PSF set.
      // This information is used to find the appropriate PSF for 
      // an incident photon with particular energy and off-axis angle.
      fits_write_key(fptr, TDOUBLE, "ENERGY", &psf->item[count].energy, 
		     "photon energy for the PSF generation in [keV]", status);
      fits_write_key(fptr, TDOUBLE, "OFFAXANG", &psf->item[count].angle, 
		     "off-axis angle in [deg]", status);

      
      HDpar_stamp(fptr, count+1, status);
      
      // END of writing header keywords.


      // Write the image to the file:
      long fpixel[2] = {x0+1, y0+1};                // lower left corner
      //                   |-----|--> FITS coordinates start at (1,1)
      // upper right corner
      long lpixel[2] = {psf->item[count].naxis1, psf->item[count].naxis2}; 
      fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, sub_psf, status);
      
    } // END of loop over individual PSF items in the storage.

  } while (0); // END of ERROR handling loop


  // close the FITS-file
  if(fptr) fits_close_file(fptr, status);

  return(*status);
}



