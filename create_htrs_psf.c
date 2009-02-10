#include <stdio.h>
#include <stdlib.h>

#include "global_constants.h"
#include "detector.h"
#include "psf.h"


int main()
{
  struct Detector detector;
  struct PSF_Store store;

  int count, count2;

  char msg[MAXMSG];
  int status = EXIT_SUCCESS;



  do { // Error handling loop 

    // Detector Setup:
    detector.type = HTRS;
    detector.width = 7;
    detector.pixelwidth = 3.0e-3;    // in [m]

    status = htrs_get_detector(&detector);



    // PSF Setup:
    store.N_elements = 1;
    store.width = 700;
    store.pixelwidth = (detector.width*detector.pixelwidth)/store.width;

    // Get memory for the PSF.
    store.psf = (struct PSF *) malloc(store.N_elements * sizeof(struct PSF));
    if (store.psf) {   // memory was allocated successfully
      for (count=0; count<store.N_elements; count++) {
	store.psf[count].data = (double **) malloc(store.width * sizeof(double *));
	if (store.psf[count].data) {
	  for (count2=0; count2<store.width; count2++) {
	    store.psf[count].data[count2] = (double *) 
	      malloc(store.width * sizeof(double));
	    if (!store.psf[count].data[count2]) {
	      status = EXIT_FAILURE;
	    }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    // Check if all memory was allocated successfully
    if (status != EXIT_SUCCESS) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg, status);  
      break;
    }


    // Use the PSF photon partition data from
    // Oosterbroek: "HTRS: sizing and positioning for IXO CDF design".
    double psf_parts[7][7];
    for(count=0; count<7; count++) {
      for(count2=0; count2<7; count2++) {
	psf_parts[count][count2] = 0.;
      }
    }
    psf_parts[0][1] = 0.010;
    psf_parts[0][2] = 0.024;
    psf_parts[0][3] = 0.024;
    psf_parts[0][4] = 0.010;

    psf_parts[1][1] = 0.024;
    psf_parts[1][2] = 0.035;
    psf_parts[1][3] = 0.035;
    psf_parts[1][4] = 0.035;
    psf_parts[1][5] = 0.024;

    psf_parts[2][0] = 0.024;
    psf_parts[2][1] = 0.035;
    psf_parts[2][2] = 0.034;
    psf_parts[2][3] = 0.034;
    psf_parts[2][4] = 0.035;
    psf_parts[2][5] = 0.024;

    psf_parts[3][0] = 0.010;
    psf_parts[3][1] = 0.035;
    psf_parts[3][2] = 0.034;
    psf_parts[3][3] = 0.023;
    psf_parts[3][4] = 0.034;
    psf_parts[3][5] = 0.035;
    psf_parts[3][6] = 0.010;

    psf_parts[4][0] = 0.024;
    psf_parts[4][1] = 0.035;
    psf_parts[4][2] = 0.034;
    psf_parts[4][3] = 0.034;
    psf_parts[4][4] = 0.035;
    psf_parts[4][5] = 0.024;

    psf_parts[5][1] = 0.024;
    psf_parts[5][2] = 0.035;
    psf_parts[5][3] = 0.035;
    psf_parts[5][4] = 0.035;
    psf_parts[5][5] = 0.024;

    psf_parts[6][1] = 0.010;
    psf_parts[6][2] = 0.024;
    psf_parts[6][3] = 0.024;
    psf_parts[6][4] = 0.010;

    // END of PSF data setup


    
    // Calculate the square PSF
    struct Point2d position;
    int x[2], y[2], normalization=0;
    double fraction[2];

    // Determine normalization
    for(count=0; count<store.width; count++) {
      for(count2=0; count2<store.width; count2++) {
	position.x = (count-store.width/2+0.5)*store.pixelwidth;
	position.y = (count2-store.width/2+0.5)*store.pixelwidth;
	htrs_get_pixel(detector, position, x, y, fraction);
	
	if ((x[0]==3)&&(y[0]==3)) {
	  normalization++;
	}
      }
    }

    for(count=0; count<store.width; count++) {
      for(count2=0; count2<store.width; count2++) {
	position.x = (count-store.width/2+0.5)*store.pixelwidth;
	position.y = (count2-store.width/2+0.5)*store.pixelwidth;
	htrs_get_pixel(detector, position, x, y, fraction);
	
	if (x[0] != INVALID_PIXEL) {
	  store.psf[0].data[count][count2] = psf_parts[x[0]][y[0]]/normalization;//6468.;
	  // Normalization (number of PSF pixels                   <-|
	  // per HTRS pixel)
	}
      }
    }

    // Store the PSF in a FITS file.
    remove("psf.htrs.image.fits");
    status = save_psf_image(store, "psf.htrs.image.fits", &status);


  } while (0); // END of error handling loop


  return(status);
}
