#include "sourceimage.h"


////////////////////////////
SourceImage* get_SourceImage() 
{
  SourceImage* si=NULL;

  // Allocate memory:
  si = (SourceImage*)malloc(sizeof(SourceImage));
  if(si!=NULL) {
    si->naxis1=0;
    si->naxis2=0;
    si->pixel=NULL;
  }

  return(si);
}



/////////////////////////////////////
SourceImage* get_SourceImage_fromFile(char* filename, int* status)
{
  SourceImage* si=NULL;
  fitsfile* fptr=NULL;
  float* input_buffer=NULL;
  char msg[MAXMSG]; 

  do { // Beginning of ERROR handling loop

    // Get an empty SourceImage using the standard Constructor without any arguments:
    si = get_SourceImage();
    if(si==NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the SourceImage!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Open image FITS file
    headas_chat(5, "open extended SourceImage FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;
    
    // Determine the width of the image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    if (naxes[0] != naxes[1]) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: SourcImage must be square!\n"); // TODO: really??
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      si->naxis1 = (int)naxes[0];
      si->naxis2 = (int)naxes[1];
    }

    // Determine the width of one image pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &si->cdelt1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &si->cdelt2, comment, status)) break;

    // Determine the WCS coordinates of the image.
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &si->crpix1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &si->crpix2, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &si->crval1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &si->crval2, comment, status)) break;

    // Convert from [deg] to [rad]
    si->cdelt1 *= M_PI/180.;
    si->cdelt2 *= M_PI/180.;
    si->crval1 *= M_PI/180.;
    si->crval2 *= M_PI/180.;

    // Determine the edges of the covered area:
    si->minra  = si->crval1 - si->cdelt1*(si->crpix1-0.5);
    si->maxra  = si->crval1 + si->cdelt1*(si->naxis1-(si->crpix1-0.5));
    si->mindec = si->crval2 - si->cdelt2*(si->crpix2-0.5);
    si->maxdec = si->crval2 + si->cdelt2*(si->naxis2-(si->crpix2-0.5));
    

    // Allocate memory for the pixels of the image:
    si->pixel = (struct SourceImagePixel**)
      malloc(si->naxis1*sizeof(struct SourceImagePixel*));
    if (si->pixel!=NULL) {
      int count;
      for(count=0; (count<si->naxis1)&&(*status==EXIT_SUCCESS); count++) {
	si->pixel[count] = (struct SourceImagePixel*)
	  malloc(si->naxis2*sizeof(struct SourceImagePixel));
	if(si->pixel[count]==NULL) {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: could not allocate memory for storing the "
		  "extended SourceImage!\n");
	  HD_ERROR_THROW(msg, *status);
	  break;
	}
      }
    } else {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing the "
	      "extended SourceImage!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } 
    
    // Allocate memory for input buffer (1D array):
    input_buffer=(float*)malloc(si->naxis1*si->naxis2*sizeof(float));
    if(input_buffer==NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing the "
	      "extended SourceImage!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    // END of memory allocation
    

    // READ the FITS image:
    int anynul;
    double null_value=0.;
    long fpixel[2] = {1, 1};   // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    long lpixel[2] = {si->naxis1, si->naxis2}; // upper right corner
    long inc[2] = {1, 1};
    if (fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure and apply the correct factors to convert the pixel
    // information (surface brightness) into real count rate for eROSITA
    // (assuming a Raymond-Smith spectrum for the clusters).
    int x, y;
    for(x=0; x<si->naxis1; x++) {
      for(y=0; y<si->naxis2; y++) {
	si->pixel[x][y].rate = input_buffer[x+ si->naxis2*y]; // [photons/s]
	si->pixel[x][y].t_next_photon = 0.;
      }
    }

  } while(0); // END of Error handling loop


  // --- clean up ---

  // free the input buffer
  if(input_buffer) free(input_buffer);

  // close FITS file (if open)
  if(fptr!=NULL) fits_close_file(fptr, status);

  return(si);
}




/////////////////////
void free_SourceImage(SourceImage* si) 
{
  if(si != NULL) {
    if(si->naxis1 > 0) {
      int count;
      for(count=0; count<si->naxis1; count++) {
	if(si->pixel[count] != NULL) free(si->pixel[count]);
      }
      free(si->pixel);
    }
    free(si);
  }
}


//////////////////////////////////////////
SourceImageCatalog* get_SourceImageCatalog() 
{
  SourceImageCatalog* sic = NULL;

  sic = (SourceImageCatalog*) malloc(sizeof(SourceImageCatalog));

  if(sic!=NULL) {
    sic->nimages = 0;
    sic->images = NULL;
  }

  return(sic);
}




////////////////////////////
void free_SourceImageCatalog(SourceImageCatalog* sic) 
{
  if (sic!=NULL) {
    if (sic->nimages > 0) {
      int count;
      for(count=0; count<sic->nimages; count++) {
	free_SourceImage(sic->images[count]);
      }
    }
    if (sic->images!=NULL) {
      free(sic->images);
    }
    free(sic);
  }
}

