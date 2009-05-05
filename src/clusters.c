#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>


#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "sixt.h"



struct ClusterPixel {
  float rate;           // count rate in this pixel
  //  double t_last_photon; // last photon emitted from this direction
};

typedef struct {
  struct ClusterPixel **pixel;
  int naxis1, naxis2;    // width of the image [pixel]
  double cdelt1, cdelt2; // width of one pixel [rad]
  double crpix1, crpix2; // [pixel]
  double crval1, crval2; // [rad]

  double minra, maxra;   // minimum right ascension covered by the image [rad]
  double mindec, maxdec; // maximum    -"-
} ClusterImage;

typedef struct {
  int nimages;
  ClusterImage* images; /* nimages */
} ClusterImageCatalog;




// Constructor for Cluster Images.
ClusterImage* get_ClusterImage() 
{
  ClusterImage* ci=NULL;

  // Allocate memory:
  ci = (ClusterImage*)malloc(sizeof(ClusterImage));
  if(ci!=NULL) {
    ci->naxis1=0;
    ci->naxis2=0;
    ci->pixel=NULL;
  }

  return(ci);
}



// Constructor: Reads a cluster image from a FITS file and stores it in the 
// corresponding data structure.
ClusterImage* get_ClusterImage_fromFile(char* filename, int* status)
{
  ClusterImage* ci=NULL;
  fitsfile* fptr=NULL;
  float* input_buffer=NULL;
  char msg[MAXMSG]; 


  do { // Beginning of ERROR handling loop

    // Get an empty ClusterImage using the Constructor:
    ci = get_ClusterImage();
    if(ci==NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing "
	      "the cluster image!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Open PSF FITS file
    headas_chat(5, "open cluster image FITS file '%s' ...\n", filename);
    if (fits_open_image(&fptr, filename, READONLY, status)) break;
    
    // Determine the width of the PSF image.
    long naxes[2];
    if (fits_get_img_size(fptr, 2, naxes, status)) break;
    if (naxes[0] != naxes[1]) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: cluster image must be square!\n"); // TODO: really??
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      ci->naxis1 = (int)naxes[0];
      ci->naxis2 = (int)naxes[1];
    }

    // Determine the width of one detector pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &ci->cdelt1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CDELT2", &ci->cdelt2, comment, status)) break;
    // Rescale from [deg] to [rad]:
    ci->cdelt1 *= M_PI/180.;
    ci->cdelt2 *= M_PI/180.;

    // From the header keywords determine the minimum and maximum
    // right ascension and declination covered by the image:
    if (fits_read_key(fptr, TDOUBLE, "CRPIX1", &ci->crpix1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRPIX2", &ci->crpix2, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL1", &ci->crval1, comment, status)) break;
    if (fits_read_key(fptr, TDOUBLE, "CRVAL2", &ci->crval2, comment, status)) break;
    // Rescale from [deg] to [rad]:
    ci->crval1 *= M_PI/180.;
    ci->crval2 *= M_PI/180.;

    // Determine the edges of the covered area:
    ci->minra  = ci->crval1 - ci->cdelt1* ci->crpix1;
    ci->maxra  = ci->crval1 + ci->cdelt1*(ci->naxis1-ci->crpix1);
    ci->mindec = ci->crval2 - ci->cdelt2* ci->crpix2;
    ci->maxdec = ci->crval2 + ci->cdelt2*(ci->naxis2-ci->crpix2);
    

    // Allocate memory for the pixels of the image:
    ci->pixel = (struct ClusterPixel**)malloc(ci->naxis1*sizeof(struct ClusterPixel*));
    if (ci!=NULL) {
      int count;
      for(count=0; (count<ci->naxis1)&&(*status==EXIT_SUCCESS); count++) {
	ci->pixel[count] = (struct ClusterPixel*)
	  malloc(ci->naxis2*sizeof(struct ClusterPixel));
	if(ci->pixel[count]==NULL) {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: could not allocate memory for storing the "
		  "cluster image!\n");
	  HD_ERROR_THROW(msg, *status);
	  break;
	}
      }
    } else {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing the "
	      "cluster image!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } 
    
    // Allocate memory for input buffer (1D array):
    input_buffer=(float*)malloc(ci->naxis1*ci->naxis2*sizeof(float));
    if(input_buffer==NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing the "
	      "cluster image!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    // END of memory allocation
    

    // READ the FITS image:
    int anynul;
    double null_value=0.;
    long fpixel[2] = {1, 1};   // lower left corner
    //                |--|--> FITS coordinates start at (1,1)
    long lpixel[2] = {ci->naxis1, ci->naxis2}; // upper right corner
    long inc[2] = {1, 1};
    if (fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure.
    int x, y;
    for(x=0; x<ci->naxis1; x++) {
      for(y=0; y<ci->naxis2; y++) {
	ci->pixel[x][y].rate = input_buffer[x+ ci->naxis2*y];
      }
    }

  } while(0); // END of Error handling loop


  // --- clean up ---

  // free the input buffer
  if(input_buffer) free(input_buffer);

  // close FITS file (if open)
  if(fptr!=NULL) fits_close_file(fptr, status);

  return(ci);
}



// Destructor: release allocated memory.
void free_ClusterImage(ClusterImage* ci) 
{
  if(ci != NULL) {
    if(ci->naxis1 > 0) {
      int count;
      for(count=0; count<ci->naxis1; count++) {
	if(ci->pixel[count] != NULL) free(ci->pixel[count]);
      }
      free(ci->pixel);
    }
    free(ci);
  }
}


// Constructor for the ClusterImageCatalog:
ClusterImageCatalog* get_ClusterImageCatalog() 
{
  ClusterImageCatalog* cic = NULL;

  cic = (ClusterImageCatalog*) malloc(sizeof(ClusterImageCatalog));

  if(cic!=NULL) {
    cic->nimages = 0;
    cic->images = NULL;
  }

  return(cic);
}


// Desctructor for the ClusterImageCatalog:
void free_ClusterImageCatalog(ClusterImageCatalog* cic) 
{
  if (cic!=NULL) {
    if (cic->nimages > 0) {
      int count;
      for(count=0; count<cic->nimages; count++) {
	//free_ClusterImage(&(cic->images[count])); // TODO
      }
    }
    free(cic);
  }
}

