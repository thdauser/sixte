#include "fudgexp.h"


int fudgexp_getpar(struct Parameters* par)
{
  char* sbuffer=NULL;

  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("PhotonList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of photon list!\n", status);
    return(status);
  } 
  strcpy(par->PhotonList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_file_name("ExposureMap", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the exposure map!\n", status);
    return(status);
  } 
  strcpy(par->ExposureMap, sbuffer);
  free(sbuffer);

  status=ape_trad_query_float("ExposureTime", &par->ExposureTime);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the exposure time!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


int fudgexp_main() {
  struct Parameters par;

  // Photon list.
  PhotonFile* plf=NULL;

  // Exposure map.
  float** map=NULL;
  long naxes[2];
  struct wcsprm* wcs;
  
  // FITS file-pointer to exposure map.
  fitsfile* fptr=NULL;   
  // String buffer for FITS header.
  char* headerstr=NULL;
  // 1-dimensional image buffer.
  float* image1d=NULL;

  int status=EXIT_SUCCESS;

  // Register HEATOOL
  set_toolname("fudgexp");
  set_toolversion("0.01");

  do { // ERROR handling loop

    // Read parameters by PIL:
    status = fudgexp_getpar(&par);
    CHECK_STATUS_BREAK(status);

    // Initialize the random number generator.
    sixt_init_rng((int)time(NULL), &status);
    CHECK_STATUS_BREAK(status);

    // Open the photon list.
    plf=openPhotonFile(par.PhotonList, READWRITE, &status);
    CHECK_STATUS_BREAK(status);

    // Load the exposure map.
    fits_open_file(&fptr, par.ExposureMap, READONLY, &status);
    CHECK_STATUS_BREAK(status);

    // Read the WCS header keywords.
    // Read the entire header into the string buffer.
    int nkeys;
    fits_hdr2str(fptr, 1, NULL, 0, &headerstr, &nkeys, &status);
    CHECK_STATUS_BREAK(status);
    // Parse the header string and store the data in the wcsprm data
    // structure.
    int nreject, nwcs;
    if (0!=wcspih(headerstr, nkeys, 0, 3, &nreject, &nwcs, &wcs)) {
      SIXT_ERROR("parsing of WCS header failed");
      status=EXIT_FAILURE;
      break;
    }
    if (nreject>0) {
      SIXT_ERROR("parsing of WCS header failed");
      status=EXIT_FAILURE;
      break;
    }

    // Determine the image dimensions.
    int naxis;
    fits_get_img_dim(fptr, &naxis, &status);
    CHECK_STATUS_BREAK(status);
    if (2!=naxis) {
      SIXT_ERROR("specified FITS HDU does not contain a 2-dimensional image");
      status=EXIT_FAILURE;
      break;
    }
    fits_get_img_size(fptr, naxis, naxes, &status);
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the image.
    map = (float**)malloc(naxes[0]*sizeof(float*));
    CHECK_NULL_BREAK(map, status, 
		     "memory allocation for exposure map failed");
    long ii;
    for (ii=0; ii<naxes[0]; ii++) {
      map[ii] = (float*)malloc(naxes[1]*sizeof(float));
      CHECK_NULL_BREAK(map[ii], status, 
		       "memory allocation for exposure map failed");
    }
    CHECK_STATUS_BREAK(status);

    // Allocate memory for the image input buffer.
    image1d=(float*)malloc(naxes[0]*naxes[1]*sizeof(float));
    CHECK_NULL_BREAK(image1d, status, 
		     "memory allocation for exposure map buffer failed");

    // Read the image from the file.
    int anynul;
    float null_value=0.;
    long fpixel[2] = {1, 1};   // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1).
    long lpixel[2] = {naxes[0], naxes[1]}; // Upper right corner.
    long inc[2]    = {1, 1};
    fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value, 
		     image1d, &anynul, &status);
    CHECK_STATUS_BREAK(status);
    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure and generate a probability distribution function,
    // i.e., sum up the pixels.
    for(ii=0; ii<naxes[0]; ii++) {
      long jj;
      for(jj=0; jj<naxes[1]; jj++) {
	map[ii][jj] = image1d[ii+ naxes[0]*jj];
      }
    }
    // END of loading the exposure map.
    
    
    // Loop over all photons in the list.
    long row=1;
    for (ii=0; ii<plf->nrows; ii++) {

      // Get the next photon from the list.
      Photon ph;
      status=PhotonFile_getRow(plf, &ph, row);
      CHECK_STATUS_BREAK(status);

      // Determine the pixel coordinates corresponding to the photon
      // directions of origin.
      double pixcrd[2];
      double imgcrd[2];
      double world[2] = {ph.ra*180./M_PI, ph.dec*180./M_PI};
      double phi, theta;
      int status2=0;
      wcss2p(wcs, 1, 2, world, &phi, &theta, imgcrd, pixcrd, &status2);
      if (3==status2) { 
	// Pixel does not correspond to valid world coordinates.
	continue;
      } else if (0!=status2) {
	SIXT_ERROR("projection failed");
	status=EXIT_FAILURE;
	break;
      }
      
      // Determine whether the photon should be discarded.
      double p=sixt_get_random_number(&status);
      CHECK_STATUS_BREAK(status);
      long xx = (long)(pixcrd[0]-0.5);
      long yy = (long)(pixcrd[1]-0.5);
      if ((xx<0)||(xx>=naxes[0]) || (yy<0)||(yy>=naxes[1]) || 
	  (p > map[xx][yy]/par.ExposureTime)) {
	// Delete the photon from the list.
	fits_delete_rows(plf->fptr, row, 1, &status);
	CHECK_STATUS_BREAK(status);
      } else {
	// Proceed with the next line.
	row++;
      }
    }
    CHECK_STATUS_BREAK(status);
    // END loop over all photons.
        
  } while(0); // End of error handling loop

  // --- Clean Up ---

  if (NULL!=fptr) fits_close_file(fptr, &status);

  if (NULL!=headerstr) free(headerstr);  
  freePhotonFile(&plf, &status);

  if (NULL!=map) {
    long ii;
    for (ii=0; ii<naxes[0]; ii++) {
      free(map[ii]);
    }
    free(map);
  }
  if (NULL!=image1d) {
    free(image1d);
  }

  wcsfree(wcs);

  // Clean up the random number generator.
  sixt_destroy_rng();

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}

