#include <stdio.h>
#include <stdlib.h>

#include "sixt.h"
#include "detectors.h"
#include "psf.h"


int main()
{
  Detector* detector;
  PSF psf;

  int count, count2;

  char msg[MAXMSG];
  int status = EXIT_SUCCESS;


  do { // Error handling loop 

    detector = get_Detector(&status);

    // Detector Setup:
    detector->type = HTRS;
    detector->width = 7;
    detector->pixelwidth = 4.0e-3;    // in [m]

    detector = htrs_get_Detector(detector, &status);


    // PSF Setup:
    psf.N_elements = 1;
    psf.width = 700;
    psf.pixelwidth = (detector->width*detector->pixelwidth)/psf.width;

    // Get memory for the PSF.
    psf.item = (PSF_Item *) malloc(psf.N_elements * sizeof(PSF_Item));
    if (psf.item != NULL) {   // memory was allocated successfully
      for (count=0; count<psf.N_elements; count++) {
	psf.item[count].energy = 1.0;  // Default values.
	psf.item[count].angle = 0.0;
	psf.item[count].data = (double **) malloc(psf.width * sizeof(double *));
	if (psf.item[count].data != NULL) {
	  for (count2=0; count2<psf.width; count2++) {
	    psf.item[count].data[count2] = (double *) 
	      malloc(psf.width * sizeof(double));
	    if (psf.item[count].data[count2] == NULL) {
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
    for(count=0; count<psf.width; count++) {
      for(count2=0; count2<psf.width; count2++) {
	position.x = (count-psf.width/2+0.5)*psf.pixelwidth;
	position.y = (count2-psf.width/2+0.5)*psf.pixelwidth;
	htrs_get_pixel(detector, position, x, y, fraction);
	
	if ((x[0]==3)&&(y[0]==3)) {
	  normalization++;
	}
      }
    }

    for(count=0; count<psf.width; count++) {
      for(count2=0; count2<psf.width; count2++) {
	position.x = (count-psf.width/2+0.5)*psf.pixelwidth;
	position.y = (count2-psf.width/2+0.5)*psf.pixelwidth;
	htrs_get_pixel(detector, position, x, y, fraction);
	
	if (x[0] != INVALID_PIXEL) {
	  psf.item[0].data[count][count2] = psf_parts[x[0]][y[0]]/normalization;
	  // Normalization (number of PSF pixels                   <-|
	  // per HTRS pixel)
	}
      }
    }

    // Store the PSF in a FITS file.
    remove("htrs.psf.hexagons.fits");
    status = save_psf_image(&psf, "htrs.psf.hexagons.fits", &status);

  } while (0); // END of error handling loop


  return(status);
}
