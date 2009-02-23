/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and converts 
// the PSF ASCII data files determined by a simulation of Peter 
// Friedrich (MPE) to a fite in FITS format.
//
/////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @data     2008/06
// @param    outputfile - filename of the FITS output file
//
/////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "strftcpy.h"

#define TOOLSUB conv_psf2fits_main
#include "headas_main.c"

#define LINELENGTH 1024     // maximum linelength in the ASCII file
#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG 256          // maximum length of an output/error message

#include "fits_psf.h"



// reads all parameters of 'conv_rosat2fits' using PIL
int conv_psf2fits_getpar(double *hew, char outputfile[]);

// does the acutal work: create FITS file with table, open ASCII file, 
// transfer sources, close files, clean up
int conv_psf2fits_work(const double sigma, const char outputfile[]);




///////////////////////////////////////////////////////////////////////
// main procedure
int conv_psf2fits_main()
{
  char outputfile[FILENAME_LENGTH]; // name of output file (FITS file)  
  int status=0;                     // error report status
  double hew;                       // on-axis HEW
             // The data read from the files are smeared by a convolution with an
             // additional Gauss function in order to include the on-axis spread of
             // the PSF, which is not included in the data.)


  // HEATOOLs: register program
  set_toolname("conv_psf2fits");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = conv_psf2fits_getpar(&hew,outputfile);

  // sigma of the Gaussian on-axis blurring of the PSF
  double sigma = hew/2. /(61.2*60.)*384. /1.177410022;   
  //                      |         |     |-> sqrt(2.*log(2.))
  //                      |         |-> width of the detector 
  //                      |             (number of pixels, eROSITA specific)
  //                      |-> transformation from arcsec in pixel coordinates 
  //                          (eROSITA specific)
  // Sigma is given in [pixels].


  if(!status) {
    // perform the actual work: open and create files respectively, 
    // transfer source data etc.
    status = conv_psf2fits_work(sigma, outputfile);

    headas_chat(5, "finished\n");
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'conv_rosat2fits' using PIL
int conv_psf2fits_getpar(
			 double *hew,       //
			 char outputfile[]  // filename of the FITS output file
			 ) 
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages

  if ((status = PILGetReal("onaxis_HEW", hew))) {
    sprintf(msg, "Error reading the on-axis HEW!\n");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the name of the outputfile!\n");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}



///////////////////////////////////////////////////////////////////////////////
// does the acutal work: create FITS file with table, open ASCII files,
//  create PSF array, store it to FITS file, and finally clean up
int conv_psf2fits_work(
		       const double sigma,      // sigma of on-axis PSF blurring
		       const char outputfile[]  // filename of the FITS output file
		       )
{
  // The following data are specific for the particular PSF simulation!

  const int det_width = 384;

  // contains the filenames of the PSF data files (ASCII files with 
  // photon event lists) for each off-axis angle and each energy
  const char filenames[PSF_N_ANGLES][PSF_N_ENERGIES][FILENAME_LENGTH] = 
    {{"psf/sim_pet_00offaxis_1keV.dat", "psf/sim_pet_00offaxis_4keV.dat", "psf/sim_pet_00offaxis_7keV.dat"},
     {"psf/sim_pet_05offaxis_1keV.dat", "psf/sim_pet_05offaxis_4keV.dat", "psf/sim_pet_05offaxis_7keV.dat"},
     {"psf/sim_pet_10offaxis_1keV.dat", "psf/sim_pet_10offaxis_4keV.dat", "psf/sim_pet_10offaxis_7keV.dat"},
     {"psf/sim_pet_15offaxis_1keV.dat", "psf/sim_pet_15offaxis_4keV.dat", "psf/sim_pet_15offaxis_7keV.dat"},
     {"psf/sim_pet_20offaxis_1keV.dat", "psf/sim_pet_20offaxis_4keV.dat", "psf/sim_pet_20offaxis_7keV.dat"},
     {"psf/sim_pet_25offaxis_1keV.dat", "psf/sim_pet_25offaxis_4keV.dat", "psf/sim_pet_25offaxis_7keV.dat"},
     {"psf/sim_pet_30offaxis_1keV.dat", "psf/sim_pet_30offaxis_4keV.dat", "psf/sim_pet_30offaxis_7keV.dat"}};

  // store the information about the index-energy and index-angle relations in the PSF structure data,
  // in order to be able look up, which energy/angle corresponds to which index:
  const double angles[PSF_N_ANGLES] = {0.,                  //  [floating point pixel]
				       5./60.,
				       10./60.,
				       15./60.,
				       20./60.,
				       25./60.,
				       30./60. };

  const double energies[PSF_N_ENERGIES] = {1., 4., 7.};     // energies in [keV]
  const double N_TOTAL_PHOTONS = 1000000.;  // total number of photons that have been simulated for the PSF event list
  //  const double PSF_THRESHOLD = 3./N_TOTAL_PHOTONS;



  fitsfile *output_fptr=NULL; // fitsfile pointer to output file
  FILE *input_fptr=NULL;      // ASCII file pointer to input file
  char line[LINELENGTH];      // input buffer for ASCII file
  char *cbuffer[1];           // string buffer
  double ****psf;             // PSF data array
  int count1, count2, count3, count4;
  int pos1, pos2;
  double normalization_factor = 0.;  // normalizing the normal distribution

  int status=EXIT_SUCCESS;    // error status
  char msg[MAXMSG];           // buffer for error output messages


  do {  // error handling loop (is only run once)

    for (count3=-4; count3<=4; count3++) {
      for (count4=-4; count4<=4; count4++) {
	normalization_factor += 
	  1./(sigma*sqrt(2.*M_PI))*exp(-0.5*(pow((double)count3/sigma,2.)+pow((double)count4/sigma,2.)));
      }
    }

    // get memory for the PSF
    psf = (double ****) malloc(PSF_N_ANGLES * sizeof(double ***));
    if (psf) {   // memory was allocated successfully
      for (count1=0; count1<PSF_N_ANGLES; count1++) {
	psf[count1] = (double ***) malloc(PSF_N_ENERGIES * sizeof(double **));
	if (psf[count1]) {
	  for (count2=0; count2<PSF_N_ENERGIES; count2++) {
	    psf[count1][count2] = (double **) malloc(det_width * sizeof(double *));
	    if (psf[count1][count2]) {
	      for (count3=0; count3<det_width; count3++) {
		psf[count1][count2][count3] = (double *) malloc(det_width * sizeof(double));
		if (!psf[count1][count2][count3]) {
		  status = EXIT_FAILURE;
		}
	      }
	    } else { status = EXIT_FAILURE; }
	  }
	} else { status = EXIT_FAILURE; }
      }
    } else { status = EXIT_FAILURE; }
    
    // check if all necessary memory was allocated successfully
    if (status == EXIT_FAILURE) {
      sprintf(msg, "Error: not enough memory to store PSF data!\n");
      HD_ERROR_THROW(msg,status);  
      // quit function, because can't store psf data, if memory allocation failed
      return(status);
    }

    // get memory for char input buffer
    cbuffer[0]=malloc(10*sizeof(char));
    if (!cbuffer[0]) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error allocating memory for input buffer!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }


    /////////////////////////////////////////////////////////
    // fill the psf.data array with data from the ASCII files
    headas_chat(5, "load PSF data ...\n");

    long row=0;      // row in the FITS table

    // loop over all different off-axis angles (break, if an error has occured)
    for (count1=0; (count1<PSF_N_ANGLES)&&(status==EXIT_SUCCESS); count1++) {
      // loop over the different energies
      for (count2=0; (count2<PSF_N_ENERGIES)&&(status==EXIT_SUCCESS); count2++) {
	
	// clear the detector array
	for (count3=0; count3<det_width; count3++) {
	  for (count4=0; count4<det_width; count4++) {
	    psf[count1][count2][count3][count4] = 0.;
	  }
	}

	// open PSF data file (ASCII)
	headas_chat(5,"import PSF file '%s' ...\n", filenames[count1][count2]);
	if (!(input_fptr=fopen(filenames[count1][count2], "r+"))) {
	  status = EXIT_FAILURE;
	  sprintf(msg, "Error opening the PSF data file '%s' (ASCII)!\n", filenames[count1][count2]);
	  HD_ERROR_THROW(msg,status);
	  break;
	}


	// loop over all events in the list
	headas_chat(5, "processing photon events ...\n");
	while (fgets(line, LINELENGTH, input_fptr)) {
	  // parse ASCII line
	  strftcpy(cbuffer[0],line,24,3);
	  pos1 = atoi(cbuffer[0]) - 1;    // convert to integer number
	  
	  strftcpy(cbuffer[0],line,31,3);
	  pos2 = atoi(cbuffer[0]) - 1;
	  
	  // add event to detector array considering the smearing effect by the optics
	  // (convolution with Gauss, HEW = 15" bzw. 1.6 pixel => sigma)
	  for (count3=-4; count3<=4; count3++) {
	    for (count4=-4; count4<=4; count4++) {
	      if ((pos1+count3>=0)&&(pos1+count3<det_width) && (pos2+count4>=0)&&(pos2+count4<det_width)) {
		psf[count1][count2][pos1][pos2] +=
		  (1./N_TOTAL_PHOTONS) * 1./(sigma*sqrt(2.*M_PI))*
		  exp(-0.5* (pow((double)count3/sigma,2.)+pow((double)count4/sigma,2.) )) / normalization_factor;
		//          normlization_factor -> correction factor, as we deal with a step function and not with a 
		//                              continuous function. The integral of the probability distribution has 
		//                              to be 1.
	      }
	    }
	  }

	} // end of loop over all events
	
	// close psf file
	if(input_fptr) fclose(input_fptr);
	input_fptr=NULL;	
	
      }  // end of loop over energies
    }  // end of loop over off-axis angles
    

    // determine size of PSF sub-rectangles
    // (don't save entire PSF but only the relevant region around the central peak, which has a 
    // probability greater than 0)
    int n = 0;     // width and
    int m = 0;     // height of sub-rectangle
 /*
    int xl, xh, yl, yh;

    for (count1=0; count1<PSF_N_ANGLES; count1++) {
      for (count2=0; count2<PSF_N_ENERGIES; count2++) {

	// get the width of the sub-rectangle
	for (xl=0; xl<det_width; xl++) {
	  for (count4=0; count4<det_width; count4++) {
	    if (psf[count1][count2][xl][count4] > PSF_THRESHOLD) {
	      break;
	    }
	  }
	  if (count4<det_width) {
	    break;
	  }
	}

	for (xh=det_width; xh>=0; xh--) {
	  for (count4=0; count4<det_width; count4++) {
	    if (psf[count1][count2][xl][count4] > PSF_THRESHOLD) {
	      break;
	    }
	  }
	  if (count4<det_width) {
	    break;
	  }
	}

	// find maximum necessary sub-rectangle
	if ((xh-xl>=0)&&(xh-xl>n)) {
	  n = xh-xl;
	}

      }
    }
*/	
    n = det_width;
    m = det_width;



    ///////////////////////////
    // create FITS file

    // delete old FITS output file
    remove(outputfile);

    // create new FITS output file
    if (fits_create_file(&output_fptr, outputfile, &status)) break;
    headas_chat(5, "FITS file '%s' created ...\n", outputfile);


    // variables for the creation of the fits-table (field-type and -format of columns)
    char *ftype[PSF_NFIELDS];
    char *fform[PSF_NFIELDS];
    char *funit[PSF_NFIELDS];
    // create a binary table in the FITS file
    psf_create_tbl_parameter(ftype, fform, funit, det_width);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, PSF_NFIELDS, ftype, fform, funit, "PSF" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status)*/

    // write headers
    fits_write_key (output_fptr, TSTRING, "TELESCOP", "eROSITA", "name of the telescope", &status);
    fits_write_key (output_fptr, TSTRING, "COMMENT", "DESCRIPT", "simulated PSF for the eROSITA Wolter telescope",&status);
    fits_write_key (output_fptr, TINT, "n", &n, "width of PSF sub-rectangle",&status);
    fits_write_key (output_fptr, TINT, "m", &m, "height of PSF sub-rectangle",&status);

    // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, &status);

    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written into FITS file '%s' ...\n", outputfile);


    // create relevant PSF sub-rectangles
    double *sub_psf=0;     // PSF sub-rectangle
    
    sub_psf = (double *) malloc((long)n*(long)m*sizeof(double));
    if (!sub_psf) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error allocating memory!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
	
    // copy sub-rectangle to 1d-array for each energy and off-axis angle
    for (count1=0; count1<PSF_N_ANGLES; count1++) {
      for (count2=0; count2<PSF_N_ENERGIES; count2++) {
	int x=0, y=0;         // coordinates of upper left corner of sub-rectangle

	for (count3=x; count3<x+n; count3++) {
	  for (count4=y; count4<y+m; count4++) {
	    sub_psf[((count3-x)*n+count4-y)] = psf[count1][count2][count3][count4];
	  }
	}
    
	// write data row to FITS file                                                    
	if ((status=insert_psf_fitsrow(angles[count1], energies[count2], x, y, sub_psf, 
				       n*m, output_fptr, ++row))) break;
      }
    }

    if (!sub_psf) free (sub_psf);

  
  } while(0);  // end of error loop


  // clean up:
  headas_chat(5, "closing files ...\n");
  // close the ASCII input file
  if(input_fptr) fclose(input_fptr);
  // close the FITS file
  if(output_fptr) fits_close_file(output_fptr, &status);

  // free PSF memory
  if (psf) {
    for (count1=0; count1<PSF_N_ANGLES; count1++) {
      if (psf[count1]) {
	for (count2=0; count2<PSF_N_ENERGIES; count2++) {
	  if (psf[count1][count2]) {
	    for (count3=0; count3<det_width; count3++) {
	      if (psf[count1][count2][count3]) {
		free (psf[count1][count2][count3]);
	      }
	    }
	    free(psf[count1][count2]);
	  }
	}
	free(psf[count1]);
      }
    }
    free(psf);
    psf=NULL;
  }

  // free memory of input buffer
  if(cbuffer[0]) free(cbuffer[0]);  

  return(status);

}



