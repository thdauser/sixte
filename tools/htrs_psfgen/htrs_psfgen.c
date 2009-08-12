#include "sixt.h"
#include "htrsdetector.h"
#include "psf.h"


int main()
{
  HexagonalPixels hexagonalPixels;
  PSF psf;

  int count, count2;
  int status = EXIT_SUCCESS;


  do { // Error handling loop 

    // Initialization of HexagonalPixels data structure.
    struct HexagonalPixelsParameters hpparameters = {
      .npixels = 37,
      .pixelwidth = 4.e-3  // in [m]
    };
    if(EXIT_SUCCESS!=(status=initHexagonalPixels(&hexagonalPixels, &hpparameters))) break;

    // PSF Setup:
    psf.N_elements = 1;

    // Get memory for the PSF.
    psf.item = (PSF_Item *) malloc(psf.N_elements * sizeof(PSF_Item));
    if (psf.item != NULL) {   // Memory was allocated successfully!
      for (count=0; count<psf.N_elements; count++) {
	psf.item[count].naxis1 = 700;
	psf.item[count].naxis2 = 700;
	psf.item[count].cdelt1 = (14*hexagonalPixels.h)/psf.item[count].naxis1;
	psf.item[count].cdelt2 = (14*hexagonalPixels.h)/psf.item[count].naxis2;
	
	psf.item[count].energy = 1.0;  // Default values.
	psf.item[count].angle = 0.0;
	psf.item[count].data = (double **) 
	  malloc(psf.item[count].naxis1 * sizeof(double *));
	if (psf.item[count].data != NULL) {
	  for (count2=0; count2<psf.item[count].naxis1; count2++) {
	    psf.item[count].data[count2] = (double *) 
	      malloc(psf.item[count].naxis2 * sizeof(double));
	    if (psf.item[count].data[count2] == NULL) {
	      status = EXIT_FAILURE;
	    }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    // Check if all memory was allocated successfully
    if (EXIT_SUCCESS!=status) {
      HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", status);  
      break;
    }


    // Use the PSF photon partition data from
    // Oosterbroek: "HTRS: sizing and positioning for IXO CDF design".
    double psf_parts[37];
    for(count=0; count<37; count++) {
      psf_parts[count] = 0.;
    }
    psf_parts[0] = 0.023;

    psf_parts[1] = 0.034;
    psf_parts[2] = 0.034;
    psf_parts[3] = 0.034;
    psf_parts[4] = 0.034;
    psf_parts[5] = 0.034;
    psf_parts[6] = 0.034;

    psf_parts[7] = 0.035;
    psf_parts[8] = 0.035;
    psf_parts[9] = 0.035;
    psf_parts[10] = 0.035;
    psf_parts[11] = 0.035;
    psf_parts[12] = 0.035;
    psf_parts[13] = 0.035;
    psf_parts[14] = 0.035;
    psf_parts[15] = 0.035;
    psf_parts[16] = 0.035;
    psf_parts[17] = 0.035;
    psf_parts[18] = 0.035;

    psf_parts[19] = 0.024;
    psf_parts[20] = 0.024;
    psf_parts[21] = 0.010;
    psf_parts[22] = 0.024;
    psf_parts[23] = 0.024;
    psf_parts[24] = 0.010;
    psf_parts[25] = 0.024;
    psf_parts[26] = 0.024;
    psf_parts[27] = 0.010;
    psf_parts[28] = 0.024;
    psf_parts[29] = 0.024;
    psf_parts[30] = 0.010;
    psf_parts[31] = 0.024;
    psf_parts[32] = 0.024;
    psf_parts[33] = 0.010;
    psf_parts[34] = 0.024;
    psf_parts[35] = 0.024;
    psf_parts[36] = 0.010;
    // END of PSF data setup

    
    // Calculate the square PSF
    struct Point2d position;
    int pixel, normalization=0;

    // Determine normalization
    for(count=0; count<psf.item[0].naxis1; count++) {
      for(count2=0; count2<psf.item[0].naxis2; count2++) {
	position.x = (count-psf.item[0].naxis1/2+0.5)*psf.item[0].cdelt1;
	position.y = (count2-psf.item[0].naxis2/2+0.5)*psf.item[0].cdelt2;
	getHexagonalPixel(&hexagonalPixels, position, &pixel);
	
	if (0==pixel) {
	  normalization++;
	}
      }
    }

    for(count=0; count<psf.item[0].naxis1; count++) {
      for(count2=0; count2<psf.item[0].naxis2; count2++) {
	position.x = (count-psf.item[0].naxis1/2+0.5)*psf.item[0].cdelt1;
	position.y = (count2-psf.item[0].naxis2/2+0.5)*psf.item[0].cdelt2;
	getHexagonalPixel(&hexagonalPixels, position, &pixel);
	
	if (INVALID_PIXEL != pixel) {
	  //                 |-------|--> x and y are exchanged to obtain the
	  //                 |       |    correct orientation in the FITS image
	  psf.item[0].data[count2][count] = psf_parts[pixel]/normalization;
	  // Normalization (number of PSF pixels per HTRS pixel) <-|
	}
      }
    }

    // Store the PSF in a FITS file.
    remove("htrs.psf.hexagons.fits");
    status = save_psf_image(&psf, "htrs.psf.hexagons.fits", &status);

  } while (0); // END of error handling loop


  // --- Clean up ---
  cleanupHexagonalPixels(&hexagonalPixels);

  return(status);
}
