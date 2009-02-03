/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and calculates the PSF
// according to a model with parameters given in an ASCII file.
//
/////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @date     2008/11
// @param    inputfile - filename of the ASCII input file
// @param    outputfile - filename of the FITS output file
//
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

// GNU scientific library: error function (calculation of gaussian integral)
#include <gsl/gsl_sf_erf.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB create_psf_main
#include "headas_main.c"

#include "psf.h"

#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG 256          // maximum length of an output/error message
//#define SIMPLE_GAUSS_PSF 1  // produce a simple on-axis Gauss PSF


// reads all parameters of 'create_psf' using PIL
int create_psf_getpar(int *psf_width, double *psf_pixelwidth, double *focal_length,
		      char inputfile[], char outputfile[]);

// does the acutal work
int create_psf_work(const int psf_width, const double psf_pixelwidth, const double,
		    const char inputfile[], const char outputfile[]);

// Determines the number of lines in an ASCII file.
int linecount(const char *filename);



/////////////////////////////////////////////////
// Main routine.
int create_psf_main()
{
  // Input parameters:
  int psf_width;
  double psf_pixelwidth;           // in [\mu m]
  double focal_length;             // in [m]
  // name of the input file containing the fit parameters in ASCII format
  char inputfile[FILENAME_LENGTH];     
  char outputfile[FILENAME_LENGTH];// name of output file (FITS file)

  int status=0;                    // error report status


  // HEATOOLs: register program
  set_toolname("create_psf");
  set_toolversion("0.01");


  // read parameters using the PIL
  status = create_psf_getpar(&psf_width, &psf_pixelwidth, &focal_length,
			     inputfile, outputfile);


  if(!status) {
    // perform the actual work: open and create files respectively, 
    // transfer source data etc.
    status = create_psf_work(psf_width, psf_pixelwidth, focal_length, 
			     inputfile, outputfile);

    headas_chat(5, "finished\n");
  }

  return(status);
}



///////////////////////////////////////////////////////////////////////
// This routine reads all parameters of 'create_psf' using the PIL.
int create_psf_getpar(
		      int *psf_width,         // width of the PSF in pixels
		      double *psf_pixelwidth, // width of a PSF pixel in [\mu m]
		      double *focal_length,   // [m]
		      // filename of the input file (ASCII) containing the 
		      // parameters of the PSF fits
		      char inputfile[],  
		      char outputfile[]       // filename of the FITS output file
		      ) 
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages


  if ((status = PILGetInt("psf_width", psf_width))) {
    sprintf(msg, "Error reading the width of the PSF!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("psf_pixelwidth", psf_pixelwidth))) {
    sprintf(msg, "Error reading the width of the PSF pixels!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("focal_length", focal_length))) {
    sprintf(msg, "Error reading the focal length of the telescope!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("inputfile", inputfile))) {
    sprintf(msg, "Error reading the name of the input file!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the name of the output file!\n");
    HD_ERROR_THROW(msg,status);
  }


  return(status);
}



///////////////////////////////////////////////////////////////////////////////////
// This routine does the acutal work: create a FITS file with a table, 
// open ASCII files, create PSF array, store it to FITS file and finally clean up.
int create_psf_work(
		    const int psf_width,
		    const double psf_pixelwidth,   // in [\mu m]
		    const double focal_length,     // in [m]
		    // filename of the input file containing fit parameters 
		    // in ASCII format
		    const char inputfile[],  
		    const char outputfile[]  // filename of the FITS output file
		    )
{
  FILE *input_fptr=NULL;   // ASCII file pointer to input file
  struct PSF_Store store;  // data structure, containing several PSFs for the 
                           // individual off-axis angles and energies
  int count1, count2, count3;

  int status=EXIT_SUCCESS;    // error status
  char msg[MAXMSG];           // buffer for error output messages


  do {  // error handling loop (is only run once)

    // Initialization:
    // Set the parameters of the PSF.
    store.width = psf_width;
    store.pixelwidth = psf_pixelwidth*1.e-6;  // transformation from [\mu m] -> [m]

    float X[store.width], Y[store.width];
    for (count1=0; count1<store.width; count1++) {
      X[count1] = (count1-store.width/2) * store.pixelwidth;
      Y[count1] = X[count1];
    }

    // Get the number of lines in the ASCII file.
#ifdef SIMPLE_GAUSS_PSF
    store.N_elements = 1;
#else
    store.N_elements = linecount(inputfile);
#endif
    if (store.N_elements <= 0) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Input file '%s' contains no data!\n", inputfile);
      HD_ERROR_THROW(msg,status);  
      return(status);
    }

    // Get memory for the PSF data.
    store.psf = (struct PSF *) malloc(store.N_elements * sizeof(struct PSF));
    if (store.psf) {   // memory was allocated successfully
      for (count1=0; count1<store.N_elements; count1++) {
	store.psf[count1].data = (double **) malloc(store.width * sizeof(double **));
	if (store.psf[count1].data) {
	  for (count2=0; count2<store.width; count2++) {
	    store.psf[count1].data[count2] = 
	      (double *) malloc(store.width * sizeof(double));
	    if (!store.psf[count1].data[count2]) {
	      status = EXIT_FAILURE;
	      break;
	    }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    
    // check if all necessary memory was allocated successfully
    if (status == EXIT_FAILURE) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg,status);  
      return(status);
    }



#ifdef SIMPLE_GAUSS_PSF
    // Create simple Gaussian PSF with a given HEW
    store.psf[0].angle = 0.;
    store.psf[0].energy = 1.;

    // Fill the PSF array with a 2D Gaussian distribution.
    double hew =  // HEW in detector pixels
      tan(5 /3600. /180. * M_PI)*focal_length  /store.pixelwidth;
    //    tan (HEW)             *focal length  /pixelwidth
    double sigma = hew /(2.*sqrt(2.*log(2.)));  // sigma in detector pixels
    headas_chat(5, "HEW: %lf pixel\n", hew);
    headas_chat(5, "PSF Sigma: %lf pixel\n", sigma);
    double x, y;
    for (count2=0; count2<store.width; count2++) {
      for (count3=0; count3<store.width; count3++) {
	x = (double)(count2-store.width/2);
	y = (double)(count3-store.width/2);
	store.psf[0].data[count2][count3] = // 0.38 * 1./196;
	  (gsl_sf_erf_Q(x/sigma) - gsl_sf_erf_Q((x+1.)/sigma))* 
	  (gsl_sf_erf_Q(y/sigma) - gsl_sf_erf_Q((y+1.)/sigma));
      }
    }

    
#else
    ////////////////////////////////////////////////////////////////////////////////
    // Load the input ASCII file and store the contained fit parameters in an array 
    // of a structure. Then calculate the corresponding PSF from these data.

    headas_chat(5, "import ASCII file '%s' with PSF parameters and calculate "
		"PSF ...\n", inputfile);

    // Structure containing the PSF fit parameters for the individual off-axis angles 
    // and photon energies.
    struct PSF_Parameter {
      int offaxis_angle;    // [arc min]
      int energy;           // [keV]
      float chi;
      
      float A1, x01, y01, sx1, sy1, prefactor1, c1, phase1;
      float A2, x02, y02, sx2, sy2, prefactor2, c2, phase2;
      float A3, x03, y03, sx3, sy3, prefactor3, c3, phase3;
    };

    // Open the ASCII file containing the PSF parameters.
    if (!(input_fptr=fopen(inputfile, "r+"))) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error accessing the ASCII file '%s'!\n", inputfile);
      HD_ERROR_THROW(msg,status);
      break;
    }

    struct PSF_Parameter psf_parameter;
    count1=0;
    while (!(feof(input_fptr))) {
      // read the PSF parameters from the ASCII file
      fscanf(input_fptr, "%d %d %e ", &psf_parameter.offaxis_angle, 
	     &psf_parameter.energy, &psf_parameter.chi);

      fscanf(input_fptr, "%e %e %e %e %e ", &psf_parameter.A1, &psf_parameter.x01, 
	     &psf_parameter.y01, &psf_parameter.sx1, &psf_parameter.sy1);
      fscanf(input_fptr, "%e %e %e ", &psf_parameter.prefactor1, &psf_parameter.c1, 
	     &psf_parameter.phase1);

      fscanf(input_fptr, "%e %e %e %e %e ", &psf_parameter.A2, &psf_parameter.x02, 
	     &psf_parameter.y02, &psf_parameter.sx2, &psf_parameter.sy2);
      fscanf(input_fptr, "%e %e %e ", &psf_parameter.prefactor2, &psf_parameter.c2, 
	     &psf_parameter.phase2);
	
      fscanf(input_fptr, "%e %e %e %e %e ", &psf_parameter.A3, &psf_parameter.x03, 
	     &psf_parameter.y03, &psf_parameter.sx3, &psf_parameter.sy3);
      fscanf(input_fptr, "%e %e %e\n", &psf_parameter.prefactor3, &psf_parameter.c3, 
	     &psf_parameter.phase3);
	

      // clear the PSF data array
      for (count2=0; count2<store.width; count2++) {
	for (count3=0; count3<store.width; count3++) {
	  store.psf[count1].data[count2][count3] = 0.;
	}
      }

      // Store the off-axis angle and the energy:
      //                                          [arc min] -> [degree]
      store.psf[count1].angle = (double)psf_parameter.offaxis_angle/60.; 
      store.psf[count1].energy = (double)psf_parameter.energy;

      // Calculate the PSF from the fit model and the given parameters.
      for (count2=0; count2<store.width; count2++) {
	for (count3=0; count3<store.width; count3++) {
	  float phi = atan2(Y[count3]-psf_parameter.y01, X[count2]-psf_parameter.x01);

	  store.psf[count1].data[count2][count3] =
	      psf_parameter.A1/(2.*M_PI*psf_parameter.sx1*psf_parameter.sy1)
	      *exp(-0.5*(pow((X[count2]-psf_parameter.x01)/psf_parameter.sx1,2.) + 
			 pow((Y[count3]-psf_parameter.y01)/psf_parameter.sy1,2.) ))
	      *(1 + psf_parameter.prefactor1*
		cos(psf_parameter.c1*phi+psf_parameter.phase1))

              +psf_parameter.A2/(2*M_PI*psf_parameter.sx2*psf_parameter.sy2)
	      *exp(-0.5*(pow((X[count2]-psf_parameter.x02)/psf_parameter.sx2,2.) + 
			 pow((Y[count3]-psf_parameter.y02)/psf_parameter.sy2,2.) ))
	      *(1 + psf_parameter.prefactor2*
		cos(psf_parameter.c2*phi+psf_parameter.phase2))

	      +psf_parameter.A3/(2*M_PI*psf_parameter.sx3*psf_parameter.sy3)
	      *exp(-0.5*(pow((X[count2]-psf_parameter.x03)/psf_parameter.sx3,2.) + 
			 pow((Y[count3]-psf_parameter.y03)/psf_parameter.sy3,2.) ))
	      *(1 + psf_parameter.prefactor3*
		cos(psf_parameter.c3*phi+psf_parameter.phase3));
	}
      }
      
      count1++;  // next line in the ASCII file
    }   // end of the loop over the lines in the input file

    // close the ASCII input file
    if(input_fptr) fclose(input_fptr);
    input_fptr=NULL;	

#endif


    // Create FITS file and store the PSF data (for all off-axis angles and 
    // energies in the same file).
    //save_psf_to_fits(store, outputfile, &status);

    // Test: output to FITS image:
    remove("test.psf.image.fits");
    save_psf_image(store, "test.psf.image.fits", &status);


  } while(0);  // END of error loop


  //---------------
  // clean up:
  headas_chat(5, "closing files ...\n");

  // close the ASCII input file
  if(input_fptr) fclose(input_fptr);

  // free PSF memory
  free_psf_store(store);

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




