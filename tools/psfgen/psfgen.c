/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "sixt.h"

// GNU scientific library: error function (calculation of gaussian integral)
#include <gsl/gsl_sf_erf.h>

#include "psf.h"

#define TOOLSUB psfgen_main
#include "headas_main.c"



/////////////////////////////////////////////////
// Main routine.
int psfgen_main()
{
  PSF* psf=NULL;

  // Output file (FITS file).
  char psffile[MAXFILENAME];

  // Focal length of the telescope [m].
  double focallength; 

  // Number of pixels along one edge of the (square) PSF image.
  int width;

  // Size of the image pixels [m].
  double pixelwidth;

  // HEW of on-axis Gaussian [deg].
  double hew; 
         
  // Parameters for the King profile.
  double rc;
  double alpha;

  // Error status.
  int status=EXIT_SUCCESS;  


  // HEATOOLs: register program
  set_toolname("psfgen");
  set_toolversion("0.03");

  do { // Beginning of outer ERROR handling loop

    // Allocate memory:
    psf=(PSF*)malloc(sizeof(PSF));
    if (NULL==psf) {
      status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for PSF data structure failed");
      break;
    }
    psf->data    =NULL;
    psf->energies=NULL;
    psf->thetas  =NULL;
    psf->phis    =NULL;


    if ((status=PILGetInt("Width", &width))) {
      SIXT_ERROR("failed reading the width of the PSF image");
      break;
    }
    if ((status=PILGetReal("PixelWidth", &pixelwidth))) {
      SIXT_ERROR("failed reading the width of the PSF pixels");
      break;
    }
    if ((status=PILGetReal("FocalLength", &focallength))) {
      SIXT_ERROR("failed reading the focal length of the telescope");
      break;
    }
    if ((status=PILGetFname("PSFFile", psffile))) {
      SIXT_ERROR("failed reading the name of PSF the output file");
      break;
    }
    int type;
    if ((status=PILGetInt("type", &type))) {
      SIXT_ERROR("failed reading the PSF type");
      break;
    }
    if (1==type) {  // Simple Gaussian PSF
      if ((status=PILGetReal("HEW", &hew))) {
	SIXT_ERROR("failed reading the HEW");
	break;
      }
      hew = hew/(3600.); // Rescaling from [arc sec] -> [deg]
      psf->nenergies = 1;
      psf->nthetas   = 1;
      psf->nphis     = 1;

    } // END of Simple Gauss PSF (type==1)
    
    else if (2==type) { // King profile.
      if ((status=PILGetReal("rc_a", &rc))) {
	SIXT_ERROR("failed reading the King parameter rc");
	break;
      }
      if ((status=PILGetReal("alpha_x", &alpha))) {
	SIXT_ERROR("failed reading the King parameter alpha");
	break;
      }
      psf->nenergies = 1;
      psf->nthetas   = 1;
      psf->nphis     = 1;
      
    } // End of type 2, King profile

    else if (3==type) { // King profile for pn Camera
      double rc_a, rc_b, rc_c, rc_d = 0.;
      double alpha_x, alpha_y, alpha_z, alpha_w = 0.;
      
      if ((status=PILGetReal("rc_a", &rc_a))) {
	SIXT_ERROR("failed reading the King parameter rc_a");
	break;
      }
      if ((status=PILGetReal("rc_b", &rc_b))) {
	SIXT_ERROR("failed reading the King parameter rc_b");
	break;
      }
      if ((status=PILGetReal("rc_c", &rc_c))) {
	SIXT_ERROR("failed reading the King parameter rc_c");
	break;
      }
      if ((status=PILGetReal("rc_d", &rc_d))) {
	SIXT_ERROR("failed reading the King parameter rc_d");
	break;
      }
      if ((status=PILGetReal("alpha_x", &alpha_x))) {
	SIXT_ERROR("failed reading the King parameter alpha_x");
	break;
      }
      if ((status=PILGetReal("alpha_y", &alpha_y))) {
	SIXT_ERROR("failed reading the King parameter alpha_y");
	break;
      }
      if ((status=PILGetReal("alpha_z", &alpha_z))) {
	SIXT_ERROR("failed reading the King parameter alpha_z");
	break;
      }
      if ((status=PILGetReal("alpha_w", &alpha_w))) {
	SIXT_ERROR("failed reading the King parameter alpha_w");
	break;
      }
      
      double Energy = 1.5;
      double Theta = 0.;
      
      rc = rc_a + rc_b*Energy + rc_c*Theta + rc_d*Energy*Theta;
      alpha = alpha_x + alpha_y*Energy + alpha_z*Theta + alpha_w*Energy*Theta; 
      printf("rc is %f and alpha is %f\n", rc, alpha);
      
      psf->nenergies = 1;
      psf->nthetas   = 1;
      psf->nphis     = 1;
      
    } // End of type 3, King profile
			
    else { // invalid PSF type
      status=EXIT_FAILURE;
      SIXT_ERROR("invalid PSF type");
      break;
    }

    // --- END of PIL parameter input ---


    // --- PSF initialization ---
    
    int count1, count2, count3;
    int xcount;

    // Get memory for the PSF data.
    psf->data = (PSF_Item ***)malloc(psf->nenergies*sizeof(PSF_Item**));
    if (NULL!=psf->data) {   // Memory was allocated successfully.
      for (count1=0; count1<psf->nenergies; count1++) {
	psf->data[count1] = (PSF_Item **)malloc(psf->nthetas*sizeof(PSF_Item*));
	if (NULL!=psf->data[count1]) {   // Memory was allocated successfully.
	  for (count2=0; count2<psf->nthetas; count2++) {
	    psf->data[count1][count2] = (PSF_Item *)malloc(psf->nphis*sizeof(PSF_Item));
	    if (NULL!=psf->data[count1][count2]) {   // Memory was allocated successfully.
	      for (count3=0; count3<psf->nphis; count3++) {
		psf->data[count1][count2][count3].naxis1 = width;
		psf->data[count1][count2][count3].naxis2 = width;
		psf->data[count1][count2][count3].cdelt1 = pixelwidth;
		psf->data[count1][count2][count3].cdelt2 = pixelwidth;

		psf->data[count1][count2][count3].data = (double**) 
		  malloc(psf->data[count1][count2][count3].naxis1 * sizeof(double**));
		if (psf->data[count1][count2][count3].data) {
		  for (xcount=0; xcount<psf->data[count1][count2][count3].naxis1; xcount++) {
		    psf->data[count1][count2][count3].data[xcount] = 
		      (double*)malloc(psf->data[count1][count2][count3].naxis2 * sizeof(double));
		    if (!psf->data[count1][count2][count3].data[xcount]) {
		      status=EXIT_FAILURE;
		      break;
		    }
		  }
		  CHECK_STATUS_BREAK(status);
		} else { status=EXIT_FAILURE; }
	      }
	      CHECK_STATUS_BREAK(status);
	    } else { status=EXIT_FAILURE; }
	  }
	  CHECK_STATUS_BREAK(status);
	} else { status=EXIT_FAILURE; }
      }
      CHECK_STATUS_BREAK(status);
    } else { status=EXIT_FAILURE; }
    // Check if all necessary memory has been allocated successfully:
    if (status!=EXIT_SUCCESS) {
      SIXT_ERROR("not enough memory to store PSF data");  
      break;
    }

    // --- END of PSF initialization ---


    // --- PSF data generation ---

    if (1==type) {
      // Create a simple Gaussian PSF with a given HEW.
      psf->energies = (double*)malloc(sizeof(double));
      psf->thetas   = (double*)malloc(sizeof(double));
      psf->phis     = (double*)malloc(sizeof(double));
      if ((NULL==psf->energies) || 
	  (NULL==psf->thetas) || 
	  (NULL==psf->phis)) {
	SIXT_ERROR("not enough memory to store PSF data");  
	break;
      }
      psf->energies[0] = 1.;
      psf->thetas[0]   = 0.;
      psf->phis[0]     = 0.;

      // Fill the PSF array with a 2D Gaussian distribution.
      // Determine sigma in unit of [detector pixels].
      double sigma = 
	hew*M_PI/180.          // [rad]
	/(2.*sqrt(2.*log(2.))) // HEW -> sigma
	/atan(psf->data[0][0][0].cdelt1/focallength); // [rad] -> [detector pixels]
      headas_chat(5, "PSF Sigma: %.2lf pixel\n", sigma);
      double x, y;
      for (count1=0; count1<psf->data[0][0][0].naxis1; count1++) {
	for (count2=0; count2<psf->data[0][0][0].naxis2; count2++) {
	  x = (double)(count1-psf->data[0][0][0].naxis1/2);
	  y = (double)(count2-psf->data[0][0][0].naxis2/2);
	  psf->data[0][0][0].data[count1][count2] = 
	    (gsl_sf_erf_Q(x/sigma) - gsl_sf_erf_Q((x+1.)/sigma))* 
	    (gsl_sf_erf_Q(y/sigma) - gsl_sf_erf_Q((y+1.)/sigma));
	}
      }
    }

    if (2==type) {
      // Create a PSF with a King profile.
      psf->energies = (double*)malloc(sizeof(double));
      psf->thetas   = (double*)malloc(sizeof(double));
      psf->phis     = (double*)malloc(sizeof(double));
      if ((NULL==psf->energies) || 
	  (NULL==psf->thetas) || 
	  (NULL==psf->phis)) {
	SIXT_ERROR("not enough memory to store PSF data");  
	break;
      }
      psf->energies[0] = 1.;
      psf->thetas[0]   = 0.;
      psf->phis[0]     = 0.;

      // Calculate the King profile.
      //double A = (alpha-1.)/(pow(rc,2.) * M_PI)
      //*psf->data[0][0][0].cdelt1*psf->data[0][0][0].cdelt2;
      double norm=0.;
      for (count1=0; count1<psf->data[0][0][0].naxis1; count1++) {
	for (count2=0; count2<psf->data[0][0][0].naxis2; count2++) {
	  double x = 
	    (count1-psf->data[0][0][0].naxis1/2+0.5)*psf->data[0][0][0].cdelt1;
	  double y = 
	    (count2-psf->data[0][0][0].naxis2/2+0.5)*psf->data[0][0][0].cdelt2;;
	  psf->data[0][0][0].data[count1][count2] = 
	    1./(pow(1.+pow(sqrt(pow(x,2.)+pow(y,2.))/rc,2.),alpha));
	  norm += psf->data[0][0][0].data[count1][count2];
	}
      }

      // Normalization.
      for (count1=0; count1<psf->data[0][0][0].naxis1; count1++) {
	for (count2=0; count2<psf->data[0][0][0].naxis2; count2++) {
	  psf->data[0][0][0].data[count1][count2] *= 1./norm;
	}
      }
    }

    else if (3==type) {
      // Create a PSF with the King Profile with rc, 
      // the King core radius and alpha, the King slope
      psf->energies = (double*)malloc(sizeof(double));
      psf->thetas   = (double*)malloc(sizeof(double));
      psf->phis     = (double*)malloc(sizeof(double));
      if ((NULL==psf->energies) || 
	  (NULL==psf->thetas) || 
	  (NULL==psf->phis)) {
	SIXT_ERROR("not enough memory to store PSF data");  
	break;
      }
      psf->energies[0] = 1.;
      psf->thetas[0]   = 0.;
      psf->phis[0]     = 0.;
      
      double x, y = 0.;
      double norm = 0.;
      for (count1=0; count1<psf->data[0][0][0].naxis1; count1++) {
	for (count2=0; count2<psf->data[0][0][0].naxis2; count2++) {
	  x = (double)(count1-psf->data[0][0][0].naxis1/2);
	  y = (double)(count2-psf->data[0][0][0].naxis2/2);
	  psf->data[0][0][0].data[count1][count2] =
	    1/(pow(1+pow(sqrt(pow(x,2)+pow(y,2))/(rc/4.1),2),alpha));//*1/(pow(1+pow(y/rc,2),alpha));
	  norm = norm + psf->data[0][0][0].data[count1][count2];
	  
	}
      }
      printf("The norm is %f\n", norm);
      for (count1=0; count1<psf->data[0][0][0].naxis1; count1++) {
	for (count2=0; count2<psf->data[0][0][0].naxis2; count2++) {
	  psf->data[0][0][0].data[count1][count2] /= norm;
	}
      }
    }
    // END of which PSF type.

    // --- END of PSF data generation ---


    // Create FITS file and store the PSF data (for 
    // all off-axis angles and energies in the same file).
    remove(psffile);
    savePSFImage(psf, psffile, &status);

  } while (0); // END of outer Error handling loop


  // --- Clean up ---

  destroyPSF(&psf);

  if (EXIT_SUCCESS==status) {
    headas_chat(3, "finished successfully!\n\n");
    return(EXIT_SUCCESS);
  } else {
    return(EXIT_FAILURE);
  }
}
