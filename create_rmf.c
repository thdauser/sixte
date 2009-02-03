/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and creates a FITS
// file containing the detector response matrix (probabilities for each
// detector pixel and PHA channel).
//
/////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @date     2008/06
// @param    outputfile - FITS file containing detector response matrix
//
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "fits_pha.h"

#define TOOLSUB create_rmf_main
#include "headas_main.c"

#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG 256          // maximum length of an output/error message
#define N_PHA_CHANNELS 3    // number of detector PHA channels
#define NFIELDS 2          // number of columns in FITS binary table


// reads all parameters using the PIL
int create_rmf_getpar(char outputfile[]);

// does the acutal work
int create_rmf_work(const char outputfile[]);



////////////////////////////////////////////////////////////////////////////////////////
// main procedure
int create_rmf_main()
{
  char outputfile[FILENAME_LENGTH];     // name of output file (FITS file)

  int status=0;             // error report status


  // HEATOOLs: register program
  set_toolname("create_rmf");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = create_rmf_getpar(outputfile);


  if(!status) {
    // perform the actual work: open and create files respectively, transfer source data etc.
    status = create_rmf_work(outputfile);

    headas_chat(5, "finished\n");
  }

  return(status);
}



///////////////////////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'conv_rosat2fits' using PIL
int create_rmf_getpar(
			  char outputfile[]  // filename of the FITS output file
			  ) 
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages

  if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the 'outputfile' parameter");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
int create_rmf_work(
		    const char outputfile[]  // filename of the FITS output file
		    )
{
  const int det_width = 384;    // width of detector [pixels]

  fitsfile *output_fptr=NULL;   // fitsfile pointer to output file
  float ***rmf=NULL;  // detector redistribution matrix
  int count1, count2, count3;

  int status=0;               // error status
  char msg[MAXMSG];           // buffer for error output messages

  do {  // error handling loop (is only run once)

    headas_chat(5, "creating detector response matrix ...\n");
    // get memory for detector response matrix
    rmf = (float***) malloc(N_PHA_CHANNELS * sizeof(float**));
    if (rmf) {
      for (count1=0; count1<N_PHA_CHANNELS; count1++) {
	rmf[count1] = (float**) malloc(det_width * sizeof(float*));
	if (rmf[count1]) {
	  for (count2=0; count2<det_width; count2++) {
	    rmf[count1][count2] = (float*) malloc(det_width * sizeof(float));
	    if (!rmf[count1][count2]) {
	      status = EXIT_FAILURE;
	      sprintf(msg, "Error: Not enough memory available to store the detector redistribution matrix!\n");
	      HD_ERROR_THROW(msg,status);
	      break;
	    }
	  }
	} else {
	  status = EXIT_FAILURE;
	  sprintf(msg, "Error: Not enough memory available to store the detector redistribution matrix!\n");
	  HD_ERROR_THROW(msg,status);
	  break;
	}
      }
    } else {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: Not enough memory available to store the detector redistribution matrix!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }


    // fill RMF with entries
    for (count1=0; count1<N_PHA_CHANNELS; count1++) {
      for (count2=0; count2<det_width; count2++) {
	for (count3=0; count3<det_width; count3++) {
	  rmf[count1][count2][count3] = 1.;
	}
      }
    }



    // delete old FITS output file
    remove(outputfile);
  
    // create new FITS output file
    if (fits_create_file(&output_fptr, outputfile, &status)) break;
    headas_chat(5, "FITS file '%s' created ...\n", outputfile);

    // variables for the creation of the fits-table (field-type and -format of columns)
    char *ftype[NFIELDS];
    char *fform[NFIELDS];
    char *funit[NFIELDS];
    // create a binary table in the FITS file
    rmf_create_tbl_parameter(ftype, fform, funit);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, NFIELDS, ftype, fform, funit, "MATRIX" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status)*/

    // write headers
    fits_write_key (output_fptr, TSTRING, "TELESCOP", "eROSITA", "name of the telescope", &status);
    fits_write_key (output_fptr, TSTRING, "MISSION", "SpecXG", "spectrumXgamma", &status);
    fits_write_key (output_fptr, TSTRING, "COMMENT", "DESCRIPT", "detector response matrix", &status);

    // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, &status);

    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written into FITS file '%s' ...\n", outputfile);
    

    // write table with energy bin information of PHA channels
    headas_chat(5, "transferring data ...\n");
    long row = 0;
    for (count1=0; count1<N_PHA_CHANNELS; count1++) {
      /*
      // transfer detector pixel array to 1dimensional array
      float matrix1d[det_width*det_width];
      for (count2=0; count2<det_width; count2++) {
	for (count3=0; count3<det_width; count3++) {
	  matrix1d[count2*det_width+count3] = detrsp_matrix[count1][count2][count3];
	}
      }
      */

      // insert new row with source data from ASCII file to FITS table
      if ((status=insert_rmf_fitsrow(count1+1, (float)count1, (float)(count1+1), output_fptr, ++row))) break;
    }
    headas_chat(3, "%ld PHA channels transferred ...\n", row);
    
  } while(0);  // end of error loop


  // clean up:
  headas_chat(5, "closing files ...\n");

  if (rmf) {
    for (count1=0; count1<N_PHA_CHANNELS; count1++) {
      if (rmf[count1]) {
	for (count2=0; count2<det_width; count2++) {
	  if (rmf[count1][count2]) {
	    free(rmf[count1][count2]);
	  }
	}
	free(rmf[count1]);
      }
    }
    free(rmf);
  }

  // close the FITS file
  if(output_fptr) fits_close_file(output_fptr, &status);

  return(status);
}

