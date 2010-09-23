#include "sourceimage.h"


SourceImage* get_SourceImage() 
{
  SourceImage* si=NULL;

  // Allocate memory:
  si = (SourceImage*)malloc(sizeof(SourceImage));
  if(si!=NULL) {
    si->naxis1=0;
    si->naxis2=0;
    si->pixel=NULL;
    si->accumulated=0;
    si->total_rate=0.;
    si->t_last_photon=0.;
  }

  return(si);
}



SourceImage* getEmptySourceImage(struct SourceImageParameters* sip, int* status)
{
  SourceImage* si=NULL;

  // Check if the requested dimensions are reasonable.
  if ((sip->naxis1<0) || (sip->naxis2<0)) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: array dimensions for SourceImage pixel array cannot "
		   "be negative!\n", *status);
    return(si);
  }

  // Obtain a bare SourceImage object from the standard constructor.
  si=get_SourceImage();
  if (NULL==si) return(si);

  // Allocate memory.
  si->pixel=(double**)malloc(sip->naxis1*sizeof(double*));
  if (NULL==si->pixel) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: could not allocate memory to store "
		   "the SourceImage!\n", *status);
    return(si);
  }
  si->naxis1 = sip->naxis1;
  int xcount, ycount;
  for(xcount=0; xcount<si->naxis1; xcount++) {
    si->pixel[xcount]=(double*)malloc(sip->naxis2*sizeof(double));
    if (NULL==si->pixel[xcount]) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SourceImage!\n", *status);
      return(si);
    }
    // Clear the pixels.
    for(ycount=0; ycount<si->naxis2; ycount++) {
      si->pixel[xcount][ycount] = 0.;
    }
  }
  si->naxis2 = sip->naxis2;
  
  // Set the properties of the pixel array.
  si->cdelt1 = sip->cdelt1;
  si->cdelt2 = sip->cdelt2;

  return(si);
}



SourceImage* get_SourceImage_fromFile(char* filename, int* status)
{
  SourceImage* si=NULL;
  fitsfile* fptr=NULL;

  do { // Beginning of ERROR handling loop

    // Open image FITS file
    headas_chat(5, "open extended SourceImage FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;
    
    // Load the image using a different constructor which accepts FITS file pointers.
    si = get_SourceImage_fromHDU(fptr, status);

  } while(0); // END of Error handling loop

  // --- clean up ---

  // close FITS file (if open)
  if(fptr!=NULL) fits_close_file(fptr, status);

  return(si);
}



SourceImage* get_SourceImage_fromHDU(fitsfile* fptr, int* status)
{
  SourceImage* si=NULL;
  double* input_buffer=NULL;
  char msg[MAXMSG];

  do { // Beginning of ERROR handling loop

    headas_chat(5, "load SourceImage from FITS file HDU ...\n");

    // Get an empty SourceImage using the standard Constructor without 
    // any arguments:
    si = get_SourceImage();
    if(si==NULL) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store "
		     "the SourceImage!\n", *status);
      break;
    }

    // Determine the width of the image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    si->naxis1 = (int)naxes[0];
    si->naxis2 = (int)naxes[1];

    sprintf(msg, " NAXIS1: %d, NAXIS2: %d\n", si->naxis1, si->naxis2);
    headas_chat(5, msg);

    // Determine the width of one image pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &si->cdelt1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &si->cdelt2, comment, status)) break;

    sprintf(msg, " CDELT1: %lf, CDELT2: %lf (degree)\n", si->cdelt1, si->cdelt2);
    headas_chat(5, msg);

    // Determine the WCS coordinates of the image.
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &si->crpix1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &si->crpix2, comment, status)) break;

    sprintf(msg, " CRPIX1: %lf, CRPIX2: %lf\n", si->crpix1, si->crpix2);
    headas_chat(5, msg);

    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &si->crval1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &si->crval2, comment, status)) break;

    sprintf(msg, " CRVAL1: %lf, CRVAL2: %lf (degree)\n", si->crval1, si->crval2);
    headas_chat(5, msg);


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
    

    // Load the spectra specified in the FITS header.
    *status = loadSpectra(fptr, &si->spectrumstore);
    if (EXIT_SUCCESS!=*status) break;


    // Allocate memory for the pixels of the image:
    si->pixel = (double**)malloc(si->naxis1*sizeof(double*));
    if (si->pixel!=NULL) {
      int count;
      for(count=0; (count<si->naxis1)&&(*status==EXIT_SUCCESS); count++) {
	si->pixel[count] = (double*)malloc(si->naxis2*sizeof(double));
	if(si->pixel[count]==NULL) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory for storing the "
			 "extended SourceImage!\n", *status);
	  break;
	}
      }
    } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory for storing the "
		     "extended SourceImage!\n", *status);
      break;
    } 
    
    // Allocate memory for input buffer (1D array):
    input_buffer=(double*)malloc(si->naxis1*si->naxis2*sizeof(double));
    if(input_buffer==NULL) {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory for storing the "
		     "extended SourceImage!\n", *status);
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
    if (fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure and generate a probability distribution function,
    // i.e., sum up the pixels.
    si->total_rate = 0.;
    int x, y;
    for(x=0; x<si->naxis1; x++) {
      for(y=0; y<si->naxis2; y++) {
	si->total_rate += input_buffer[x+ si->naxis1*y]; // [photons/s]
	si->pixel[x][y] = si->total_rate;
      }
    }

    // Normalization of the SourceImage such that the sum
    // of all pixel values is 1 and create a probability 
    // distribution.
    for(x=0; x<si->naxis1; x++) {
      for(y=0; y<si->naxis2; y++) {
	si->pixel[x][y] *= 1./si->total_rate;
      }
    }
    // Set accumulation flag.
    si->accumulated = 1;

  } while(0); // END of Error handling loop


  // --- Clean Up ---

  // Free the input buffer.
  if(input_buffer) free(input_buffer);

  return(si);
}



void saveSourceImage(SourceImage* si, char* filename, int* status)
{
  fitsfile *fptr=NULL;
  double *image1d=NULL;

  // Print information to STDOUT.
  char msg[MAXMSG];
  sprintf(msg, "Store SourceImage in file '%s' ...\n", filename);
  headas_chat(5, msg);
 
  do { // ERROR handling loop

    // If the specified file already exists, remove the old version.
    remove(filename);

    // Create a new FITS-file:
    if (fits_create_file(&fptr, filename, status)) break;

    // Allocate memory for the 1-dimensional image buffer (required for
    // output to FITS file).
    image1d = (double*)malloc(si->naxis1*si->naxis2*sizeof(double));
    if (!image1d) {
      *status = EXIT_FAILURE;
      HD_ERROR_THROW("Error allocating memory!\n", *status);
      break;
    }

    if (1==si->accumulated) {
      // TODO Convert the probability distribution in the SourceImage
      // to a probability density.
      printf("Error: You  must update the function saveSourceImage()!\n");
      exit(0);
    }

    // Store the source image in the 1-dimensional buffer to handle it 
    // to the FITS routine.
    int x, y;
    for (x=0; x<si->naxis1; x++) {
      for (y=0; y<si->naxis2; y++) {
	image1d[(x+ si->naxis1*y)] = si->pixel[x][y];
      }
    }
    
    // Create an image in the FITS-file (primary HDU):
    long naxes[2] = {(long)(si->naxis1), (long)(si->naxis2)};
    if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, status)) break;
    //                                   |-> naxis
    //    int hdutype;
    if (fits_movabs_hdu(fptr, 1, NULL, status)) break;


    // Write the header keywords for the SourceImage.
    if (fits_write_key(fptr, TSTRING, "HDUCLAS1", "Image", "", status)) break;
    if (fits_write_key(fptr, TSTRING, "CTYPE1", "RA---TAN",
		       "sky coordinate system", status)) break;
    if (fits_write_key(fptr, TSTRING, "CTYPE2", "DEC--TAN",
		       "sky coordinate system", status)) break;
      
    if (fits_write_key(fptr, TSTRING, "CUNIT1", "degree", "", status)) break;
    if (fits_write_key(fptr, TSTRING, "CUNIT2", "degree", "", status)) break;
    if (fits_write_key(fptr, TDOUBLE, "CRPIX1", &si->crpix1, 
                       "X axis reference pixel", status)) break;
    if (fits_write_key(fptr, TDOUBLE, "CRPIX2", &si->crpix2, 
    		       "Y axis reference pixel", status)) break;

    double dbuffer;
    dbuffer = si->crval1 *180./M_PI;
    if (fits_write_key(fptr, TDOUBLE, "CRVAL1", &dbuffer, 
		       "coord of X ref pixel", status)) break;
    dbuffer = si->crval2 *180./M_PI;
    if (fits_write_key(fptr, TDOUBLE, "CRVAL2", &dbuffer, 
    		       "coord of Y ref pixel", status)) break;
    dbuffer = si->cdelt1 *180./M_PI;
    if (fits_write_key(fptr, TDOUBLE, "CDELT1", &dbuffer, 
		       "X axis increment", status)) break;
    dbuffer = si->cdelt2 *180./M_PI;
    if (fits_write_key(fptr, TDOUBLE, "CDELT2", &dbuffer, 
		       "Y axis increment", status)) break;
    
    HDpar_stamp(fptr, 1, status);
    if (EXIT_SUCCESS!=*status) break;
    // END of writing header keywords.


    // Write the image to the file:
    long fpixel[2] = {1, 1};  // Lower left corner.
    //                |--|--> FITS coordinates start at (1,1)
    // Upper right corner.
    long lpixel[2] = {si->naxis1, si->naxis2}; 
    fits_write_subset(fptr, TDOUBLE, fpixel, lpixel, image1d, status);

  } while (0); // END of ERROR handling loop

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);

  if (NULL!=image1d) free(image1d);
}



void free_SourceImage(SourceImage* si) 
{
  if(si != NULL) {
    if((si->naxis1 > 0)&&(NULL!=si->pixel)) {
      int count;
      for(count=0; count<si->naxis1; count++) {
	if(si->pixel[count] != NULL) free(si->pixel[count]);
      }
      free(si->pixel);
    }
    free(si);
  }
}



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



int addSourceImage2Catalog(SourceImageCatalog* sic, fitsfile* fptr) 
{
  int status=EXIT_SUCCESS;

  // Check if the SourceImageCatalog is empty.
  // If yes, we have to use malloc to get the memory, otherwise we
  // use realloc the resize the formerly allocated memory.
  if (0==sic->nimages) {
    // Allocate memory.
    sic->images = (SourceImage**)malloc(sizeof(SourceImage*));
    if (NULL==sic->images) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for ClusterImageCatalog failed!\n", status);
      return(status);
    }
    sic->nimages=1; // Initial value.
	
  } else { // SourceImageCatalog already contains SourceImage objects.
    // Resize the formerly allocated memory.
    sic->images = (SourceImage**)realloc(sic->images, (sic->nimages+1)*sizeof(SourceImage*));
    if (NULL==sic->images) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for ClusterImageCatalog failed!\n", status);
      return(status);
    }
    sic->nimages++;
  } // END of memory allocation.

  // Load the cluster image in the current HDU.
  sic->images[sic->nimages-1] = get_SourceImage_fromHDU(fptr, &status);

  return(status);
}



void getRandomSourceImagePixel(SourceImage* si, int* x, int* y) 
{
  double rnd = sixt_get_random_number();

  // Perform a binary search to obtain the x-coordinate.
  int high = si->naxis1-1;
  int low = 0;
  int mid;
  int ymax = si->naxis2-1;
  while (high > low) {
    mid = (low+high)/2;
    if (si->pixel[mid][ymax] < rnd) {
      low = mid+1;
    } else {
      high = mid;
    }
  }
  *x = low;

  // Search for the y coordinate:
  high = si->naxis2-1;
  low = 0;
  while (high > low) {
    mid = (low+high)/2;
    if (si->pixel[*x][mid] < rnd) {
      low = mid+1;
    } else {
      high = mid;
    }
  }
  *y = low;
  // Now x and y have pixel positions [integer pixel].
}

