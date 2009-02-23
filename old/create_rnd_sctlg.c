/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and creates 
// a specified number of randomly distributed sources on the sky,
// using an exponential logN-logS distribution with a specified powerlaw index.
// The sources a written to an FITS source file.
//
/////////////////////////////////////////////////////////////////////////
//
// @author  Christian Schmid
// @data    2008/04
// @param   n_sources - number of sources to be created
// @param   output-file - filename of the FITS output file
// @param   powerlaw_index - index of the powerlaw source distribution (logN-logS)
//
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "create_rnd_source.h"
#include "fits_ctlg.h"
#include "strftcpy.h"

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#define TOOLSUB create_rnd_sctlg_main
#include "headas_main.c"

#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG  256         // maximum length of an output/error message



// reads all parameters of 'create_rnd_sctlg' using PIL
int create_rnd_sctlg_getpar(int *n_sources, double *powerlaw_index, double *lower_threshold, double *upper_threshold, char outputfile[]);

// does the actual work: create the requested number of random sources and 
// add them to FITS source catalog
int create_rnd_sctlg_work(const int n_sources, const double powerlaw_index, const double lower_threshold, const double upper_threshold, const char outputfile[]);



////////////////////////////////////////////////////////////////////////////
// main procedure
int create_rnd_sctlg_main() 
{
  int n_sources;                        // number of sources, the program should create
  double powerlaw_index;                // powerlaw index for logN-logS distribution
  double l_thres, u_thres;              // lower and upper threshold for created sources [erg cm^-2 s^-1]
  char outputfile[FILENAME_LENGTH];     // name of the outputfile

  int status=0;                         // error status


  // HEATOOLs: register program
  set_toolname("create_rnd_sctlg");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = create_rnd_sctlg_getpar(&n_sources, &powerlaw_index, &l_thres, &u_thres, outputfile);


  if(!status) {
    // perform the actual work: create the requested number of random sources
    // and insert them to the FITS source catalog
    status = create_rnd_sctlg_work(n_sources, powerlaw_index, l_thres, u_thres, outputfile);

    headas_chat(5, "finished\n");
  }



  return(status);
}




///////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'create_rnd_sctlg' using PIL
int create_rnd_sctlg_getpar(
			    int *n_sources,         // number of sources to be created
			    double *powerlaw_index, // index of the powerlaw of the source distribution (logN-logS)
			    double *l_thres,        // lower and upper threshold for created sources
			    double *u_thres,
			    char outputfile[]       // name of the FITS output file
			    )
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages

  if ((status = PILGetInt("n_sources", n_sources))) {
    sprintf(msg, "Error reading the 'number of sources' parameter");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("powerlaw_index", powerlaw_index))) {
    sprintf(msg, "Error reading the 'powerlaw_index' parameter");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("lower_threshold", l_thres))) {
    sprintf(msg, "Error reading the 'lower threshold' parameter");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetReal("upper_threshold", u_thres))) {
    sprintf(msg, "Error reading the 'upper threshold' parameter");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the 'outputfile' parameter");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}




////////////////////////////////////////////////////////////////////////////
// does the actual work: create the requested number of random sources and 
// add them to FITS source catalog
int create_rnd_sctlg_work(
			  const int n_sources,         // number of sources to be created
			  const double powerlaw_index, // index of the powerlaw of the source distribution (logN-logS)
			  const double lower_threshold, // lower and
			  const double upper_threshold, // upper threshold for created sources [erg cm^-2 s^-1]
			  const char outputfile[]      // filename of the FITS output file
			  ) 
{
  fitsfile *output_fptr=NULL;      // FITS file pointer to the output file
  double rasc, dec;                // right ascension and declination of a new source
  float countrate, flux;           // countrate of a new source
  long rows=0;                     // row counter for the FITS file

  int status=0;                    // error status


  do {   // error handling loop (is only run once)

    // delete old FITS output file
    remove(outputfile);

    // create new FITS file
    if (fits_create_file(&output_fptr, outputfile, &status)) break;
    headas_chat(5, "FITS file '%s' created ...\n", outputfile);
    
    // format-data to create the new fits table
    char *ftype[N_SOURCE_FIELDS];
    char *fform[N_SOURCE_FIELDS];
    char *funit[N_SOURCE_FIELDS];
    // create a binary table in the FITS file
    create_srctbl_parameters(ftype, fform, funit);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, N_SOURCE_FIELDS, ftype, fform, funit, "RND_SOURCES" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status) */


    // write description data into the header of the fits-file
    if (fits_write_key(output_fptr, TSTRING, "COMMENT", "SOURCE CATALOGUE","source catalogue for eROSITA simulation", &status)) break;
    if (fits_write_key(output_fptr, TSTRING, "COMMENT", "random sources","randomly distributed sources", &status)) break;
    if (fits_write_key(output_fptr, TSTRING, "COMMENT", "r.a., dec, count rate", "right ascension, declination, photon cps", &status)) break;

  // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, &status);

    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written into FITS '%s' file ...\n", outputfile);


    // initialize random number generator
    HDmtInit(1);

    // create sources that are distributed randomly on the surface of the unit sphere,
    // with powerlaw intensity distribution
    headas_chat(5, "creating random sources ...\n");

    int counter;
    // conversion factor from count rate to energy flux for photon index 2:
    const float conversion_factor = 1.e7 *   log(2./0.5) /      (1./0.5 - 1./10.) * 1.6022e-16 / 2471;
    //                              J->erg   reference energy band                  keV -> J     eROSITA collecting area

    for (counter=0; counter<n_sources; ) {
      // create random source with intensity (countrate) from powerlaw distribution
      flux = get_rnd_intensity(powerlaw_index) * lower_threshold;

      //      printf("%e\n", flux);
      // check if count rate is outside desired range
      if ((flux<lower_threshold) || (flux>upper_threshold)) continue;

      // convert energy flux to count rate
      countrate = flux/conversion_factor;

      // create random position on unit sphere
      create_rnd_source_position(&rasc, &dec);
      rasc = rasc*180./M_PI;
      dec = dec*180./M_PI;

      // insert it into FITS table
      add_srctbl_row(output_fptr, rows++, rasc, dec, countrate, &status);

      // increase counter after source creation (If no source is created, the counter remains as it is.)
      counter++;
    }

    headas_chat(1, "%d sources created in the flux range from %e to %e ...\n", counter, lower_threshold,upper_threshold);

  } while(0);  // end of error handling loop


  // clean up:
  headas_chat(5, "closing files ...\n");

  // release random number generator
  HDmtFree();

  // close the FITS file
  if (output_fptr) fits_close_file(output_fptr, &status);

  return(status);
}


