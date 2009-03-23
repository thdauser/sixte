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
  int width;         // width of the image [pixel]
  double pixelwidth; // width of one pixel [rad]
  double minra, maxra;   // minimum right ascension covered by the image [rad]
  double mindec, maxcec; // maximum    -"-
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
    ci->width=0;
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

    // Allocate memory:
    ci = get_ClusterImage();
    if(ci==NULL) {
      *status=EXIT_FAILURE;
      sprintf(msg, "Error: could not allocate memory for storing the cluster image!\n");
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
      sprintf(msg, "Error: cluster image must be square!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    } else {
      ci->width = (int)naxes[0];
    }

    // Determine the width of one detector pixel.
    char comment[MAXMSG]; // buffer
    if (fits_read_key(fptr, TDOUBLE, "CDELT1", &ci->pixelwidth, comment, status)) break;

    ci->pixelwidth = 3.32/3600. * M_PI/180.; // [rad]

    // Allocate memory for the pixels of the image:
    ci->pixel = (struct ClusterPixel**)malloc(ci->width*sizeof(struct ClusterPixel*));
    if(ci!=NULL) {
      int count;
      for(count=0; (count<ci->width)&&(*status==EXIT_SUCCESS); count++) {
	ci->pixel[count] = (struct ClusterPixel*)
	  malloc(ci->width*sizeof(struct ClusterPixel));
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
    input_buffer=(float*)malloc(ci->width*ci->width*sizeof(float));
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
    long lpixel[2] = {ci->width, ci->width}; // upper right corner
    long inc[2] = {1, 1};
    if (fits_read_subset(fptr, TFLOAT, fpixel, lpixel, inc, &null_value, 
			 input_buffer, &anynul, status)) break;

    
    // Transfer the image from the 1D input buffer to the 2D pixel array in
    // the data structure.
    int x, y;
    for(x=0; x<ci->width; x++) {
      for(y=0; y<ci->width; y++) {
	ci->pixel[x][ci->width-1-y].rate = input_buffer[x*ci->width + y];
	//	ci->pixel[x][ci->width-1-y].t_last_photon = 0.;
      }
    }


  } while(0); // END of Error handling loop


  // --- clean up ---

  if(input_buffer) free(input_buffer);

  // close FITS file (if open)
  if(fptr!=NULL) fits_close_file(fptr, status);

  return(ci);
}



// Destructor: release allocated memory.
void free_ClusterImage(ClusterImage* ci) 
{
  if(ci != NULL) {
    if(ci->width > 0) {
      int count;
      for(count=0; count<ci->width; count++) {
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

