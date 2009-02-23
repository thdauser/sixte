/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and creates a FITS
// file containing a spectrum, which can be used as input for the
// measurement simulation.
//
/////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @date     2008/08
// @param    outputfile - FITS file containing the spectrum
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
#include "detector.h"
#include "global_constants.h"

#define TOOLSUB create_spectrum_main
#include "headas_main.c"

#define NFIELDS 2              // number of columns in FITS binary table for spectrum


// reads all parameters using PIL
int create_spectrum_getpar(char outputfile[], char rmffile[]);

// does the acutal work: create FITS file with table and store spectrum
int create_spectrum_work(const char outputfile[], const char rmffile[]);



////////////////////////////////////////////////////////////////////////////////////////
// main procedure
int create_spectrum_main()
{
  char outputfile[FILENAME_LENGTH];     // name of output file (FITS file)
  char rmffile[FILENAME_LENGTH];        // name of the RMF file 
                                        // (containing the energy bins of PHA channels in EBOUNDS table)

  int status=0;             // error report status


  // HEATOOLs: register program
  set_toolname("create_spectrum");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = create_spectrum_getpar(outputfile, rmffile);


  if(!status) {
    // perform the actual work: open and create files respectively, transfer source data etc.
    status = create_spectrum_work(outputfile, rmffile);

    headas_chat(5, "finished\n");
  }

  return(status);
}



///////////////////////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'conv_rosat2fits' using PIL
int create_spectrum_getpar(
			   char outputfile[],  // filename of the FITS output file
			   char rmffile[]      // filename of the RMF file 
			                       // (containing the energy bins in the EBOUNDS table)
			   ) 
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages

  if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the 'outputfile' parameter");
    HD_ERROR_THROW(msg,status);
  } 

  else if ((status = PILGetFname("rmffile", rmffile))) {
    sprintf(msg, "Error reading the name of the RMF file");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// does the acutal work: create FITS file with table, open ASCII file, transfer sources, close files, clean up
int create_spectrum_work(
			 const char outputfile[],  // filename of the FITS output file
			 const char rmffile[]
			 )
{
  fitsfile *output_fptr=NULL;      // fitsfile pointer to output file
  struct Ebounds ebounds;          // energy bins of PHA channels
  int N_channels;                  // number of PHA channels
  int count;

  int status=EXIT_SUCCESS;         // error status
  //  char msg[MAXMSG];                // buffer for error output messages


  do {  // error handling loop (is only run once)
    headas_chat(5, "reading RMF file ...\n");

    // get the energy bins of the PHA channels
    if ((status=get_ebounds(&ebounds, &N_channels, rmffile))!=EXIT_SUCCESS) break;


    headas_chat(5, "creating spectrum ...\n");

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
    spectrum_create_tbl_parameter(ftype, fform, funit);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, NFIELDS, ftype, fform, funit, "SPECTRUM" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status)*/

    // write headers
    fits_write_key (output_fptr, TSTRING, "TELESCOP", "eROSITA", "name of the telescope", &status);
    fits_write_key (output_fptr, TSTRING, "MISSION", "SpecXG", "spectrumXgamma", &status);
    fits_write_key (output_fptr, TSTRING, "COMMENT", "DESCRIPT", "sample source spectrum", &status);

    // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, &status);

    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written into FITS file '%s' ...\n", outputfile);
    

    // write table with energy bin information of PHA channels
    headas_chat(5, "creating data ...\n");
    long row = 0;
    for (count=0; count<N_channels; count++) {
      // insert new row into spectrum
      float factor = 1.;
      
      /*
      if (count < 500) {
	factor = 1.;
      } else if ((count >= 500) && (count < 1100)) {
	factor = 4.188;
      } else if ((count >= 1100) && (count < 2000)) {
	factor = 10.275;
      } else {
	factor = 0.;
      }
      */
      //factor = factor * 1./pow((double)(count+1), 2.);   // power law shape
      
      float amplitude = 1. * (ebounds.row[count].E_max -ebounds.row[count].E_min) * factor;
      if ((status=insert_spectrum_fitsrow(count, amplitude, output_fptr, ++row))) break;
    }
    headas_chat(3, "%ld PHA channels created ...\n", row);
    
  } while(0);  // end of error loop


  // clean up:
  headas_chat(5, "closing files ...\n");

  free_ebounds(ebounds);

  // close the FITS file
  if(output_fptr) fits_close_file(output_fptr, &status);

  return(status);
}

