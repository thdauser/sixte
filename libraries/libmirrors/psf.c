#include "psf.h"


int get_psf_pos(struct Point2d* position, Photon photon, 
		struct Telescope telescope, 
		Vignetting* vignetting, PSF* psf)
{
  // Calculate the off-axis angle ([rad])
  double theta = acos(scalar_product(&telescope.nz, &photon.direction));
  // and the azimuthal angle ([rad]) of the source position.
  double phi = atan2(scalar_product(&telescope.ny, &photon.direction), 
		     scalar_product(&telescope.nx, &photon.direction));

  // Get a random number to determine a random hitting position
  double rnd = get_random_number();
  //  if (rnd > psf_item->data[psf_item->naxis1-1][psf_item->naxis2-1]) {
  if (rnd > get_Vignetting_Factor(vignetting, photon.energy, theta, phi)) {
    // The photon does not hit the detector at all (e.g. it is absorbed).
    return(0);
  }
  // Otherwise the photon hits the detector.

  // Determine, which PSF should be used for that particular source 
  // direction and photon energy.
  int index1, index2, index3;
  for (index1=0; index1<psf->nenergies-1; index1++) {
    if (psf->energies[index1+1] > photon.energy) break;
  }
  for (index2=0; index2<psf->nthetas-1; index2++) {
    if (psf->thetas[index2+1] > theta) break;
  }
  for (index3=0; index3<psf->nphis-1; index3++) {
    if (psf->phis[index3+1] > phi) break;
  }

  // Perform a linear interpolation between the next fitting PSF images.
  if (index1 < psf->nenergies-1) {
    if (get_random_number() < (photon.energy-psf->energies[index1])/
	(psf->energies[index1+1]-psf->energies[index1])) {
      index1++;
    }
  }
  if (index2 < psf->nthetas-1) {
    if (get_random_number() < (theta-psf->thetas[index2])/
	(psf->thetas[index2+1]-psf->thetas[index2])) {
      index2++;
    }
  }
  if (index3 < psf->nphis-1) {
    if (get_random_number() < (phi-psf->phis[index3])/
	(psf->phis[index3+1]-psf->phis[index3])) {
      index3++;
    }
  }

  // Set a pointer to the PSF image used for the determination of the photon
  // impact position in order to enable faster access.
  PSF_Item* psf_item = &psf->data[index1][index2][index3];

  // Get a random position from the best fitting PSF image.

  // Perform a binary search to determine a random position:
  // -> one binary search for each of the 2 coordinates x and y
  rnd = get_random_number();

  // This section is only necessary for PSFs that are not normalized to 1,
  // i.e., contain some vignetting effects. According to the OGIP recommmendation 
  // for PSF images, this should not be the case. So this part of code can be removed,
  // if only proper PSFs are used.
  if (rnd > psf_item->data[psf_item->naxis1-1][psf_item->naxis2-1]) {
    // The photon does not hit the detector at all (e.g. it is absorbed).
    headas_chat(1, "### Warning: PSF contains vignetting contributions!\n");
    return(0);
  }

  // PSF coordinates [pixel] of the position obtained from the best fitting PSF image.
  int x1, y1;   

  // Perform a binary search to obtain the x-coordinate.
  int high = psf_item->naxis1-1;
  int low = 0;
  while (high-low > 1) {
    if (psf_item->data[(low+high)/2][0] < rnd) {
      low = (low+high)/2;
    } else {
      high = (low+high)/2;
    }
  }
  if (psf_item->data[high][0] > rnd) {
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
 
  // Determine the distance ([m]) of the central reference position from the optical axis
  // according to the off-axis angle theta:
  double distance = telescope.focal_length * tan(theta); // TODO *(-1) ???

  // Add the relative position obtained from the PSF image (randomized pixel indices x1 and y1).
  double x2 = distance + 
    ((double)x1 -psf_item->crpix1 +0.5 +get_random_number())*psf_item->cdelt1 
    + psf_item->crval1; // [m]
  double y2 = 
    ((double)y1 -psf_item->crpix2 +0.5 +get_random_number())*psf_item->cdelt2 
    + psf_item->crval2; // [m]

  // Rotate the postition [m] according to the azimuth angle.
  position->x = cos(phi)*x2 - sin(phi)*y2;
  position->y = sin(phi)*x2 + cos(phi)*y2;

  return(1);  
}



void free_psf(PSF *psf /**< Pointer to the PSF data structure. */)
{
  if (NULL!=psf) {
    if (NULL!=psf->data) {
      int count1, count2, count3, xcount;
      for (count1=0; count1<psf->nenergies; count1++) {
	if (NULL!=psf->data[count1]) {
	  for (count2=0; count2<psf->nthetas; count2++) {
	    if (NULL!=psf->data[count1][count2]) {
	      for (count3=0; count3<psf->nphis; count3++) {
		if (NULL!=psf->data[count1][count2][count3].data) {
		  for (xcount=0; xcount<psf->data[count1][count2][count3].naxis1; xcount++) {
		    if (NULL!=psf->data[count1][count2][count3].data[xcount]) {
		      free(psf->data[count1][count2][count3].data[xcount]);
		    }
		  }
		  free(psf->data[count1][count2][count3].data);
		}
	      }
	      free(psf->data[count1][count2]);
	    }
	  }
	  free(psf->data[count1]);
	}
      }
      free(psf->data);
    }

    if (NULL!=psf->energies) free(psf->energies);
    if (NULL!=psf->thetas  ) free(psf->thetas  );
    if (NULL!=psf->phis    ) free(psf->phis    );

    free(psf);
  }
}



/** Add a double value to a list. Before adding the value, check
    whether it is already in the list. In that case it doesn't have to
    be added. The number of list entries is modified appropriately. */
static void addDValue2List(double value, double** list, int* nvalues, int* status)
{
  // Check whether the value is already in the list.
  int count=0;
  for (count=0; count<*nvalues; count++) {
    if (fabs((*list)[count]-value) <= fabs(value*1.e-6)) return;
  }
  // The value is not in the list. So continue with the following code.

  // Adapt the amount of allocated memory.
  *list = (double*)realloc(*list, ((*nvalues)+1)*sizeof(double));
  if (NULL==*list) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory for list in PSF data structure!\n",
		   *status);
    return;
  }
  (*nvalues)++; // Now we have one element more in the list.

  // Add the value at the end of the list.
  (*list)[*nvalues-1] = value;
}



/** Sort a list of double values with nvalues elements. Apply Bubble Sort algorithm. */
static void sortDList(double* list, int nvalues)
{
  int count, index;
  double buffer;
  int exchanged; // Flag if 2 elements have been exchanged during one loop cycle.
  for (count=0; count<nvalues-1; count++) {
    exchanged=0;
    for (index=0; index<nvalues-1-count; index++) {
      if (list[index] > list[index+1]) {
	// Exchange the 2 elements.
	buffer        = list[index+1];
	list[index+1] = list[index];
	list[index]   = buffer;
	// Set the exchanged flag.
	exchanged=1;
      }
    }
    // In this loop run no elements have been exchanged, so we can break the loop here.
    if (0==exchanged) break;
  }  
}



PSF* get_psf(const char* filename, int* status)
{
  PSF* psf;
  fitsfile* fptr=NULL;   // FITSfile-pointer to PSF file
  double* data=NULL;     // Input buffer (1D array)
  long count, count2, count3;
  
  do {  // Beginning of ERROR handling loop.

    // Allocate memory for PSF data structure:
    psf = (PSF*)malloc(sizeof(PSF));
    if (NULL==psf) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for PSF structure failed!\n", *status);
      break;
    }

    // Set default initial values.
    psf->data     = NULL;
    psf->energies = NULL;
    psf->thetas   = NULL;
    psf->phis     = NULL;
    psf->nenergies = 0;
    psf->nthetas   = 0;
    psf->nphis     = 0;

    // Open PSF FITS file
    headas_chat(5, "open PSF FITS file '%s' ...\n", filename);
    if (fits_open_file(&fptr, filename, READONLY, status)) break;
    // Get the total number of HDUs in the FITS file.
    int nhdus=0; // Number of HDUs (may also contain table extensions).
    if (fits_get_num_hdus(fptr, &nhdus, status)) break;

    // Loop over all extensions in the FITS file. Check whether
    // it's any image extensions and store the energy, off-axis angle,
    // and azimuthal angle.
    int hdu, hdu_type=0;
    for (hdu=0; hdu<nhdus; hdu++) {
      // Move to the right HDU.
      if (fits_movabs_hdu(fptr, hdu+1, &hdu_type, status)) break;
      // Check if the HDU is an image extension.
      if (IMAGE_HDU==hdu_type) {
	// Get the header keywords specifying the energy, off-axis angle,
	// on azimuthal angle. 
	double energy=0., theta=0., phi=0.;
	char comment[MAXMSG];
	if (fits_read_key(fptr, TDOUBLE, "ENERGY", &energy, comment, 
			  status)) break; // [eV]
	if (fits_read_key(fptr, TDOUBLE, "THETA", &theta, comment, 
			  status)) break; // [arc min]
	if (fits_read_key(fptr, TDOUBLE, "PHI", &phi, comment, 
			  status)) break; // [deg]
	// Convert to appropriate units.
	energy *= 1.e-3;         // [eV] -> [keV];
	theta  *= M_PI/180./60.; // [arc min] -> [rad]
	phi    *= M_PI/180.;     // [deg] -> [rad]
	
	// Add value to the list of available values.
	addDValue2List(energy, &(psf->energies), &psf->nenergies, status);
	if (EXIT_SUCCESS!=*status) break;
	addDValue2List(theta , &(psf->thetas)  , &psf->nthetas  , status);
	if (EXIT_SUCCESS!=*status) break;
	addDValue2List(phi   , &(psf->phis)    , &psf->nphis    , status);
	if (EXIT_SUCCESS!=*status) break;
      }
    }
    if (EXIT_SUCCESS!=*status) break;

    // Sort the lists. 
    sortDList(psf->energies, psf->nenergies);
    sortDList(psf->thetas  , psf->nthetas  );
    sortDList(psf->phis    , psf->nphis    );

    // Plot debug information about available energies, off-axis
    // angles, and azimuthal angles.
    headas_chat(5, "PSF - available energies:\n");
    for (count=0; count<psf->nenergies; count++) {
      headas_chat(5, " %.3lf keV\n", psf->energies[count]);
    }
    headas_chat(5, "PSF - available off-axis angles:\n");
    for (count=0; count<psf->nthetas; count++) {
      headas_chat(5, " %.3lf arc min\n", psf->thetas[count]/M_PI*180.*60.);
    }
    headas_chat(5, "PSF - available azimuthal angles:\n");
    for (count=0; count<psf->nphis; count++) {
      headas_chat(5, " %.3lf deg\n", psf->phis[count]/M_PI*180.);
    }

    // Allocate memory for the 3-dimensional PSF data array.
    psf->data=(PSF_Item***)malloc(psf->nenergies*sizeof(PSF_Item**));
    if (NULL==psf->data) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for PSF data failed!\n",
		     *status);
      break;
    }
    for(count=0; count<psf->nenergies; count++) {
      psf->data[count]=(PSF_Item**)malloc(psf->nthetas*sizeof(PSF_Item*));
      if (NULL==psf->data[count]) {
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error: memory allocation for PSF data failed!\n",
		       *status);
	break;
      }
      for (count2=0; count2<psf->nthetas; count2++) {
	psf->data[count][count2]=(PSF_Item*)malloc(psf->nphis*sizeof(PSF_Item));
	if (NULL==psf->data[count][count2]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: memory allocation for PSF data failed!\n",
			 *status);
	  break;
	}
	// Initialize the PSF_Item objects in the 3-dimensional array.
	for (count3=0; count3<psf->nphis; count3++) {
	  psf->data[count][count2][count3].data = NULL;
	  psf->data[count][count2][count3].naxis1 = 0;
	  psf->data[count][count2][count3].naxis2 = 0;
	}
      }
      if (EXIT_SUCCESS!=*status) break;
    }
    if (EXIT_SUCCESS!=*status) break;
    
    // Loop over all HDUs. Read the data and add PSF to the
    // 3-dimensional array.
    for (hdu=0; hdu<nhdus; hdu++) {
      // Move to the right HDU.
      if (fits_movabs_hdu(fptr, hdu+1, &hdu_type, status)) break;
      // Check if the HDU is an image extension.
      if (IMAGE_HDU==hdu_type) {
	// Get the header keywords specifying the energy, off-axis angle,
	// on azimuthal angle. 
	double energy=0., theta=0., phi=0.;
	char comment[MAXMSG];
	if (fits_read_key(fptr, TDOUBLE, "ENERGY", &energy, comment, 
			  status)) break; // [eV]
	if (fits_read_key(fptr, TDOUBLE, "THETA", &theta, comment, 
			  status)) break; // [arc min]
	if (fits_read_key(fptr, TDOUBLE, "PHI", &phi, comment, 
			  status)) break; // [deg]
	// Convert to appropriate units.
	energy *= 1.e-3;         // [eV] -> [keV];
	theta  *= M_PI/180./60.; // [arc min] -> [rad]
	phi    *= M_PI/180.;     // [deg] -> [rad]

	// Find the right indices for this parameter configuration.
	int index1, index2, index3;
	for (index1=0; index1<psf->nenergies; index1++) {
	  if (fabs(energy-psf->energies[index1]) <= fabs(energy*1.e-6)) break;
	}
	if (index1==psf->nenergies) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Could not find appropriate PSF entry!\n",
			 *status);
	  break;
	}
	for (index2=0; index2<psf->nthetas; index2++) {
	  if (fabs(theta-psf->thetas[index2]) <= fabs(theta*1.e-6)) break;	
	}
	if (index2==psf->nthetas) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Could not find appropriate PSF entry!\n",
			 *status);
	  break;
	}
	for (index3=0; index3<psf->nphis; index3++) {
	  if (fabs(phi-psf->phis[index3]) <= fabs(phi*1.e-6)) break;
	}
	if (index3==psf->nphis) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Could not find appropriate PSF entry!\n",
			 *status);
	  break;
	}

	// Determine the width of the PSF image.
	long naxes[2];
	if (fits_get_img_size(fptr, 2, naxes, status)) break;
	psf->data[index1][index2][index3].naxis1 = (int)naxes[0];
	psf->data[index1][index2][index3].naxis2 = (int)naxes[1];

	// Determine the WCS keywords of the PSF array from the header.
	if (fits_read_key(fptr, TDOUBLE, "CDELT1", &(psf->data[index1][index2][index3].cdelt1), 
			  comment, status)) break;
	if (fits_read_key(fptr, TDOUBLE, "CDELT2", &(psf->data[index1][index2][index3].cdelt2), 
			  comment, status)) break;
	if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &(psf->data[index1][index2][index3].crpix1), 
			  comment, status)) break;
	if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &(psf->data[index1][index2][index3].crpix2), 
			  comment, status)) break;
	if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &(psf->data[index1][index2][index3].crval1), 
			  comment, status)) break;
	if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &(psf->data[index1][index2][index3].crval2), 
			  comment, status)) break;
	// TODO Check whether units of PSF image are pixels.

	// Get memory for the PSF_Item data.
	psf->data[index1][index2][index3].data = (double **)
	  malloc(psf->data[index1][index2][index3].naxis1*sizeof(double *));
	if (NULL!=psf->data[index1][index2][index3].data) {
	  for (count=0; count<psf->data[index1][index2][index3].naxis1; count++) {
	    psf->data[index1][index2][index3].data[count] = (double *)
	      malloc(psf->data[index1][index2][index3].naxis2*sizeof(double));
	    if (NULL==psf->data[index1][index2][index3].data[count]) {
	      *status = EXIT_FAILURE;
	      HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", *status);  
	      break;
	    }
	  }
	} else { 
	  *status = EXIT_FAILURE; 
	  HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", *status);  
	  break;
	}
	// Check if all memory was allocated successfully
	if (*status != EXIT_SUCCESS) break;


	// Allocate memory for input buffer (1D array)
	data=(double*)realloc(data, psf->data[index1][index2][index3].naxis1 
			      *psf->data[index1][index2][index3].naxis2
			      *sizeof(double));
	if (NULL==data) {
	  *status = EXIT_FAILURE;
	  HD_ERROR_THROW("Error: not enough memory for PSF input buffer!\n", *status);
	  break;
	}


	// Read the PSF image to the input buffer.
	int anynul;
	double null_value=0.;
	long fpixel[2] = {1, 1};   // lower left corner
	//                |--|--> FITS coordinates start at (1,1)
	// upper right corner
	long lpixel[2] = {psf->data[index1][index2][index3].naxis1, 
			  psf->data[index1][index2][index3].naxis2};  
	long inc[2] = {1, 1};
	
	if (fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, 
			     data, &anynul, status)) break;

	// Create a partition function from the 1D PSF data array,
	// i.e., sum up the individual probabilites.
	// The partition function is more adequate for determining 
	// a random photon impact position on the detector.
	double sum=0.;
	for (count=0; count<psf->data[index1][index2][index3].naxis1; count++) {
	  for (count2=0; count2<psf->data[index1][index2][index3].naxis2; count2++) {
	    // Take care of choosing x- and y-axis properly!
	    sum += data[count2*psf->data[index1][index2][index3].naxis1+count];  
	    psf->data[index1][index2][index3].data[count][count2]=sum;
	  }
	}

	// Plot normalization of PSF for current off-axis angle and energy
	headas_chat(5, "PSF: images %.2lf%% of incident photons for "
		    "%.1lf keV, %.4lf arc min, %.4lf deg, \n", sum * 100., 
		    psf->energies[index1], 
		    psf->thetas[index2]/M_PI*180.*60.,
		    psf->phis[index3]/M_PI*180.);

      } // END of check if IMAGE_HDU.
    } // END of loop over all HDUs.
    if (EXIT_SUCCESS!=*status) break;
  } while(0);  // END of error handling loop

  // Close PSF file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  // Free memory of input buffer.
  if (NULL!=data) free(data);  

  return(psf);
}



int save_psf_image(PSF* psf, const char *filename, int *status)
{
  int nhdus=0; // Number of HDUs.
  double *sub_psf=NULL;
  fitsfile *fptr=NULL;

  do { // ERROR handling loop

    // Create a new FITS-file:
    if (fits_create_file(&fptr, filename, status)) break;

    // Loop over the different PSFs in the storage:
    int index1, index2, index3; // Counters for energies, off-axis angles, and azimuthal angles.
    for (index1=0; index1<psf->nenergies; index1++) {
      for (index2=0; index2<psf->nthetas; index2++) {
	for (index3=0; index3<psf->nphis; index3++) {
	  
	  // Determine size of PSF sub-rectangles (don't save entire PSF but only 
	  // the relevant region around the central peak, which has a probability 
	  // greater than 0).
	  int n = 0;     // width and
	  int m = 0;     // height of sub-rectangle
	  n = psf->data[index1][index2][index3].naxis1; // TODO
	  m = psf->data[index1][index2][index3].naxis2;

	  // Create the relevant PSF sub-rectangle:
	  sub_psf = (double *)realloc(sub_psf, n*m*sizeof(double));
	  if (!sub_psf) {
	    *status = EXIT_FAILURE;
	    HD_ERROR_THROW("Error allocating memory!\n", *status);
	    break;
	  }
	  // Store the PSF in the 1D array to handle it to the FITS routine.
	  int x0=0, y0=0; // coordinates of lower left corner of sub-rectangle
	  int x, y;
	  for (x=x0; x<x0+n; x++) {
	    for (y=y0; y<y0+m; y++) {
	      sub_psf[((x-x0)*n+y-y0)] = psf->data[index1][index2][index3].data[x][y];
	    }
	  }
    

	  // Create an image in the FITS-file (primary HDU):
	  long naxes[2] = {(long)(psf->data[index1][index2][index3].naxis1), 
			   (long)(psf->data[index1][index2][index3].naxis2)};
	  if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
	  //                                    |-> naxis
	  nhdus++;
	  int hdutype;
	  if (fits_movabs_hdu(fptr, nhdus, &hdutype, status)) break;


	  // Write the header keywords for PSF FITS-files (CAL/GEN/92-027):
	  if (fits_write_key(fptr, TSTRING, "CTYPE1", "DETX",
			     "detector coordinate system", status)) break;
	  if (fits_write_key(fptr, TSTRING, "CTYPE2", "DETY",
			     "detector coordinate system", status)) break;
	  
	  if (fits_write_key(fptr, TSTRING, "HDUCLASS", "OGIP",
			     "Extension is OGIP defined", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUDOC", "CAL/GEN/92-020",
			     "Document containing extension definition", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUVERS", "1.0.0",
			     "giving the version of the format", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUCLAS1", "Image",
			     "", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUCLAS2", "PSF",
			     "", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUCLAS3", "PREDICTED",
			     "", status)) break;
	  if (fits_write_key(fptr, TSTRING, "HDUCLAS4", "NET",
			     "", status)) break;
      
	  if (fits_write_key(fptr, TSTRING, "CUNIT1", "pixel", "", status)) break;
	  if (fits_write_key(fptr, TSTRING, "CUNIT2", "pixel", "", status)) break;
	  double dbuffer = psf->data[index1][index2][index3].naxis1*0.5;
	  if (fits_write_key(fptr, TDOUBLE, "CRPIX1", &dbuffer, 
			     "X axis reference pixel", status)) break;
	  if (fits_write_key(fptr, TDOUBLE, "CRPIX2", &dbuffer, 
			     "Y axis reference pixel", status)) break;
	  dbuffer = 0.;
	  if (fits_write_key(fptr, TDOUBLE, "CRVAL1", &dbuffer, 
			     "coord of X ref pixel", status)) break;
	  if (fits_write_key(fptr, TDOUBLE, "CRVAL2", &dbuffer, 
			     "coord of Y ref pixel", status)) break;
	  if (fits_write_key(fptr, TDOUBLE, "CDELT1", &psf->data[index1][index2][index3].cdelt1, 
			     "X axis increment", status)) break;
	  if (fits_write_key(fptr, TDOUBLE, "CDELT2", &psf->data[index1][index2][index3].cdelt2, 
			     "Y axis increment", status)) break;
      
	  dbuffer = 0.0;
	  if (fits_write_key(fptr, TDOUBLE, "BACKGRND", &dbuffer, 
			     "background count rate per pixel", status)) break;

	  // Mission
	  if (fits_write_key(fptr, TSTRING, "TELESCOP", "", // eROSITA
			     "Satellite", status)) break;
	  if (fits_write_key(fptr, TSTRING, "INSTRUME", "", // pnCCD1
			     "Instrument", status)) break;
	  if (fits_write_key(fptr, TSTRING, "FILTER", "NONE",
			     "Filter", status)) break;

	  // Write the ENERGY, THETA, and PHI for this particular PSF set.
	  // This information is used to find the appropriate PSF for 
	  // an incident photon with particular energy and off-axis angle.
	  dbuffer = psf->energies[index1]*1000.;
	  if (fits_write_key(fptr, TDOUBLE, "ENERGY", &dbuffer, 
			     "photon energy for the PSF generation in [eV]", status)) break;
	  dbuffer = psf->thetas[index2]*180.*60./M_PI;
	  if (fits_write_key(fptr, TDOUBLE, "THETA", &dbuffer, 
			     "off-axis angle in [arc min]", status)) break;
	  dbuffer = psf->phis[index3]*180./M_PI;
	  if (fits_write_key(fptr, TDOUBLE, "PHI", &dbuffer, 
			     "azimuthal angle in [degree]", status)) break;

	  HDpar_stamp(fptr, nhdus, status);
	  if (EXIT_SUCCESS!=*status) break;
	  // END of writing header keywords.


	  // Write the image to the file:
	  long fpixel[2] = {x0+1, y0+1};  // Lower left corner.
	  //                   |-----|--> FITS coordinates start at (1,1)
	  // Upper right corner.
	  long lpixel[2] = {psf->data[index1][index2][index3].naxis1, 
			    psf->data[index1][index2][index3].naxis2}; 
	  fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, sub_psf, status);
      
	}
      }
    } // END of loops over individual PSF items in the storage.

  } while (0); // END of ERROR handling loop


  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  if (NULL!=sub_psf) free(sub_psf);

  return(*status);
}



