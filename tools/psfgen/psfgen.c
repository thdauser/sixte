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
  // Output file (FITS file)
  char psf_filename[FILENAME_LENGTH];
  double focal_length; // [m]
  // HEW of on-axis Gaussian [deg].
  double hew;          
  // Number of pixels along one edge of the (square) PSF image.
  int psf_width;
  double psf_pixelwidth;
  double rc;
  double alpha;

  // Buffer for error output messages.
  char msg[MAXMSG];        
  // Error status.
  int status=EXIT_SUCCESS;  


  // HEATOOLs: register program
  set_toolname("psfgen");
  set_toolversion("0.01");

  do { // Beginning of outer ERROR handling loop

    // Allocate memory:
    psf=(PSF*)malloc(sizeof(PSF));
    if (NULL==psf) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: memory allocation for PSF data structure failed!", status);
      break;
    }
    psf->data    =NULL;
    psf->energies=NULL;
    psf->thetas  =NULL;
    psf->phis    =NULL;


    if ((status = PILGetInt("psf_width", &psf_width))) {
      sprintf(msg, "Error reading the width of the PSF!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if ((status = PILGetReal("psf_pixelwidth", &psf_pixelwidth))) {
      sprintf(msg, "Error reading the width of the PSF pixels!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if ((status = PILGetReal("focal_length", &focal_length))) {
      sprintf(msg, "Error reading the focal length of the telescope!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if ((status = PILGetFname("psf_filename", psf_filename))) {
      sprintf(msg, "Error reading the name of PSF the output file!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    int type;
    if ((status = PILGetInt("type", &type))) {
      sprintf(msg, "Error reading the PSF type!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    if (type == 1) {  // Simple Gaussian PSF
      if ((status = PILGetReal("hew", &hew))) {
	sprintf(msg, "Error reading the HEW!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      hew = hew/(3600.);  // Rescaling from [arc sec] -> [deg]
      psf->nenergies = 1;
      psf->nthetas   = 1;
      psf->nphis     = 1;

    } // END of Simple Gauss PSF (type==1)
    
    else if (type == 3) { // King profile for pn Camera
      double rc_a, rc_b, rc_c, rc_d = 0.;
      double alpha_x, alpha_y, alpha_z, alpha_w = 0.;
      
      if ((status = PILGetReal("rc_a", &rc_a))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("rc_b", &rc_b))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("rc_c", &rc_c))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("rc_d", &rc_d))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("alpha_x", &alpha_x))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("alpha_y", &alpha_y))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("alpha_z", &alpha_z))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
	break;
      }
      if ((status = PILGetReal("alpha_w", &alpha_w))) {
	sprintf(msg, "Error reading the King parameter file!\n");
	HD_ERROR_THROW(msg,status);
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
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Invalid PSF type!\n");
      HD_ERROR_THROW(msg,status);
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
		psf->data[count1][count2][count3].naxis1 = psf_width;
		psf->data[count1][count2][count3].naxis2 = psf_width;
		psf->data[count1][count2][count3].cdelt1 = psf_pixelwidth;
		psf->data[count1][count2][count3].cdelt2 = psf_pixelwidth;

		psf->data[count1][count2][count3].data = (double**) 
		  malloc(psf->data[count1][count2][count3].naxis1 * sizeof(double**));
		if (psf->data[count1][count2][count3].data) {
		  for (xcount=0; xcount<psf->data[count1][count2][count3].naxis1; xcount++) {
		    psf->data[count1][count2][count3].data[xcount] = 
		      (double*)malloc(psf->data[count1][count2][count3].naxis2 * sizeof(double));
		    if (!psf->data[count1][count2][count3].data[xcount]) {
		      status = EXIT_FAILURE;
		      break;
		    }
		  }
		  if (EXIT_SUCCESS!=status) break;
		} else { status = EXIT_FAILURE; }
	      }
	      if (EXIT_SUCCESS!=status) break;
	    } else { status = EXIT_FAILURE; }
	  }
	  if (EXIT_SUCCESS!=status) break;
	} else { status = EXIT_FAILURE; }
      }
      if (EXIT_SUCCESS!=status) break;
    } else { status = EXIT_FAILURE; }
    // Check if all necessary memory has been allocated successfully:
    if (status != EXIT_SUCCESS) {
      HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", status);  
      break;
    }

    // --- END of PSF initialization ---


    // --- PSF data generation ---

    if (type == 1) {
      // Create a simple Gaussian PSF with a given HEW.
      psf->energies = (double*)malloc(sizeof(double));
      psf->thetas   = (double*)malloc(sizeof(double));
      psf->phis     = (double*)malloc(sizeof(double));
      if ((NULL==psf->energies) || (NULL==psf->thetas) || (NULL==psf->phis)) {
	HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", status);  
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
	/atan(psf->data[0][0][0].cdelt1/focal_length); // [rad] -> [detector pixels]
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
    else if (type == 3) {
      // Create a PSF with the King Profile with rc, the King core radius and alpha, the King slope
      psf->energies = (double*)malloc(sizeof(double));
      psf->thetas   = (double*)malloc(sizeof(double));
      psf->phis     = (double*)malloc(sizeof(double));
      if ((NULL==psf->energies) || (NULL==psf->thetas) || (NULL==psf->phis)) {
	HD_ERROR_THROW("Error: not enough memory to store PSF data!\n", status);  
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


    // Create FITS file and store the PSF data (for all off-axis angles and 
    // energies in the same file).
    remove(psf_filename);
    savePSFImage(psf, psf_filename, &status);

  } while (0); // END of outer Error handling loop


  // --- Clean up ---

  destroyPSF(&psf);

  if(status==EXIT_SUCCESS) headas_chat(5, "finished successfully\n\n");
  
  return(status);
}



///////////////////////////////////////////////////////////////////////
// This function determines the number of lines in an ASCII file.
int linecount(const char *filename) {
  int lines = 0;
  FILE *file = fopen(filename, "r"); // open the file

  if (file) {
    int ch, prev = '\n' /* so empty files have no lines */;
    
    while ( (ch = fgetc(file)) != EOF ) { /* Read all chars in the file. */
      if ( ch == '\n' ) {
	++lines; /* Bump the counter for every newline. */
      }

      prev = ch; /* Keep a copy to later test whether... */
    }

    fclose(file);

    if ( prev != '\n' ) { /* ...the last line did not end in a newline. */
      ++lines; /* If so, add one more to the total. */
    }
  } else {
    lines = -1;
  }

  return(lines);
}



